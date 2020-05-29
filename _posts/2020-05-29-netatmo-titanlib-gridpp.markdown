---
layout: post
title:  "Using Netatmo measurements in weather forecasting"
date:   2020-05-29 15:58:36 +0200
author: Thomas Nipen (thomasn@met.no)
tags: netatmo optimal_interpolation quality_control
---

Networks of personal weather stations have grown immensely and now represent a major source of weather
observations. MET-Norway has since 2018 used [Netatmo](https://netatmo.com)'s network of observations in
correcting temperature forecasts on Yr (https://www.yr.no). The weather forecasts on Yr, are based on a
sophisticated weather model. However, in many cases, the model output deviates significantly from the observed
weather. Netatmo observations are used to correct these forecasts.

In this post, we will shown how to do this using two open source software packages titanlib and gridpp.
Titanlib is used to remove observations that are likely to be invalid, whereas gridpp creates a gridded
temperature field representing the current conditions. Both libraries are written in C++, but they have python
interfaces, which will be used here.

## Getting started
First install [titanlib](https://github.com/metno/titanlib) and [gridpp](https://github.com/metno/gridpp),
which if you already have their dependencies could be as easy as:

{% highlight bash %}
pip3 install titanlib
pip3 install gridpp
{% endhighlight %}

You will also need to download two test files in NetCDF format: a [file containing
observations](https://thredds.met.no//thredds/fileServer/metusers/thomasn/gridpp/obs.nc); and a [file
containing weather model
output](https://thredds.met.no//thredds/fileServer/metusers/thomasn/gridpp/analysis.nc).

## Quality control of observations

Most of Netatmo's weather stations have reliable measurements, however occasionally stations are located in
directy sunlight or too close to buildings. These must be filtered out. Let's start by reading observations
and their station's metadata from file:

{% highlight python %}
import titanlib
import netCDF4

with netCDF4.Dataset('obs.nc', 'r') as file:
    obs_lats = file.variables['latitude'][:]
    obs_lons = file.variables['longitude'][:]
    obs_elevs = file.variables['altitude'][:]
    obs = file.variables['air_temperature_2m'][:]
{% endhighlight %}

Titanlib supports a variety of quality control methods, but we will focus on the spatial consistency test
(SCT) here. The SCT compares each observations to what is expected given the other observations in the nearby
area. If the deviation is large, the observation is removed. The SCT has several parameters, which are
described in detail on the [titanlib wiki](https://github.com/metno/titanlib/wiki/Spatial-consistency-test).
For a given observation, all other observations within the `outer_radius` (in meters) will be used to compute
a cross-validation value for the observation. For computational efficiency, the cross-validation results will
be reused for all observations within the `inner_radius` such that the proceedure doesn't have to be
performed independently for every observation. To further reduce computation time, only the `num_max` closest
observations are used, even if there are more than this many observations within the outer radius. The test
is only performed if there are at least `num_min` observations within the outer radius.

Let's set up the parameters and run the test:

{% highlight python %}
inner_radius = 50000
outer_radius = 100000
num_min = 5
num_max = 100
num_iterations = 1
num_min_prof = 20
dzmin = 100
dhmin = 10000
dz = 200
t2pos = np.full(len(obs), 4)
t2neg = np.full(len(obs), 4)
eps2 = np.full(len(obs), 0.5)

flags, sct, rep = titanlib.sct(obs_lats, obs_lons, obs_elevs, obs, num_min, num_max, inner_radius, outer_radius, num_iterations, num_min_prof, dzmin, dhmin , dz, t2pos, t2neg, eps2)
{% endhighlight %}

The `flags` array now has the same length as the input observations and uses a value of `0` denoting
observations that passed the test and a value of `1` for those that are flagged. We will only keep the
non-flagged observations in the next step:

{% highlight python %}
index_valid_obs = np.where(flags == 0)[0]
index_invalid_obs = np.where(flags != 0)[0]
{% endhighlight %}

Let's plot the observations, marking flagged ones with a black edge:

{% highlight python %}
import matplotlib.pylab as mpl
mpl.scatter(obs_lons[index_valid_obs], obs_lats[index_valid_obs],
        c=obs[index_valid_obs] - 273.15, s=20, cmap="RdBu_r")
mpl.scatter(obs_lons[index_invalid_obs], obs_lats[index_invalid_obs],
        c=obs[index_invalid_obs] - 273.15, s=20, edgecolors="k", cmap="RdBu_r")
mpl.xlim(0, 35)
mpl.ylim(55, 75)
mpl.gca().set_aspect(2)
cb = mpl.colorbar()
cb.set_label(r"Temperature ($\degree C$)")
mpl.show()
{% endhighlight %}

The SCT has identified one fault observation in the southern part of the domain:

![Result of the quality control]({{ site.url }}/assets/img/titan_sct.png)

## Creating the analysis

Once we have a set of trustworth observations, we can assimilate those into the NWP background. Gridpp
provides a function that performs optimal interpolation (OI), which merges a background field and a set of
observations, based on their relative accuracies. We will first use the deterministic OI scheme, which takes
a single ensemble member. We will take the control member, which has index 0 in the background file.

{% highlight python %}
import gridpp
import numpy as np

with netCDF4.Dataset('analysis.nc', 'r') as file:
    index_control = 0
    blats = file.variables['latitude'][:]
    blons = file.variables['longitude'][:]
    belevs = file.variables['surface_geopotential'][0, 0, index_control, :, :] / 9.81
    bgrid = gridpp.Grid(blats, blons, belevs)
    background = file.variables['air_temperature_2m'][0, 0, index_control, :, :]

points = gridpp.points(obs_lats[index_valid_obs], obs_lons[index_valid_obs], obs_elevs[index_valid_obs])
{% endhighlight %}

OI does not require you to specify the error variance of the background and the observations. Instead, only
their ratios is needed. Since we trust observations more than the background, we set the ratio to 0.1. OI
also requires the background values at the observation points. In many cases a simple method like nearest
neighbour or bilinear is sufficient, however the gridpp functions allows you to calculate any values if
desired.

{% highlight python %}
variance_ratios = 0.1 * np.ones(points.size())
pbackground = gridpp.bilinear(grid, points, background)
{% endhighlight %}

Optimal interpolation needs a structure function. We will use a Barnes function (Barnes 1973), with a horizontal
decorrelation length of 100 km and a vertical decorrelation length of 200 m. This sets the limits to how far an
observation will affect the adjustment of the background.
{% highlight python %}
h = 100000
v = 200
structure = gridpp.BarnesStructure(h, v)
{% endhighlight %}

When OI processes a gridpoint, it needs to invert a matrix containing all the observations. In areas where
the observation network is dense, this matrix can be come large and will result in large computation times.
However, in most cases reducing the number of observations doesn't negatively impact the accuracy of the
analysis. We can set this using the `max_points` argument:

{% highlight python %}
max_points = 10
{% endhighlight %}

Now we have defined all the required inputs and we can create a gridded analysis:
{% highlight python %}
analysis = gridpp.optimal_interpolation(bgrid, background, points,
        obs[index_valid_obs], variance_ratios, pbackground, structure, max_points)
{% endhighlight %}

## Plotting the result

Finally, we can plot the analysis increments:
{% highlight python %}
import matplotlib.pylab as mpl

diff = analysis - background
mpl.pcolormesh(blons, blats, diff, cmap="RdBu_r", vmin=-2, vmax=2)
mpl.xlim(0, 35)
mpl.ylim(55, 75)
mpl.gca().set_aspect(2)
mpl.show()
{% endhighlight %}

![Analysis increment map]({{ site.url }}/assets/img/analysis_increment.png)

## Ensemble mode

Gridpp also supports an Ensemble-based Statistical Interpolation (EnSI; Lussana et al. 2019) scheme that takes
uses spatial structure information from an ensemble of NWP model runs. This is the method used in the
operational temperature analyses used on Yr.no, as described by Nipen et al. 2020. To test this, you need to
load the full ensemble and ensure that the ensemble dimension is at the end. Also note that EnSI requires the
specification of the observation variance (`psigmas`).

{% highlight python %}
with netCDF4.Dataset('analysis.nc', 'r') as file:
    background_ens = file.variables['air_temperature_2m'][0, 0, :, :, :]
    background_ens = np.moveaxis(background_ens, 0, 2)

num_members = background_ens.shape[2]
pbackground_ens = np.zeros([points.size(), num_members])
for e in range(num_members):
    pbackground_ens[:, e] = gridpp.bilinear(bgrid, points, background_ens[:, :, e])

psigmas = 0.5 * np.ones(points.size())

analysis_ens = gridpp.optimal_interpolation_ensi(bgrid, background_ens, points,
        obs[index_valid_obs], psigmas, pbackground_ens, structure, max_points)
diff = (analysis_ens - background_ens)[:, :, index_control]
{% endhighlight %}

![Analysis increment map for ensemble]({{ site.url }}/assets/img/analysis_increment_ens.png)

## References

Barnes, S. L., 1973: Mesoscale objective map analysis using weighted time-series observations. NOAA Tech.
Memo. ERL NSSL-62 [NTIS COM-73-10781], March, National Severe Storms Laboratory, Norman, Oklahoma, 60 pp.

Lussana, C, Seierstad, IA, Nipen, TN, Cantarello, L. Spatial interpolation of two‐metre temperature over
Norway based on the combination of numerical weather prediction ensembles and in situ observations. Q J R
Meteorol Soc. 2019; 145: 3626– 3643 ([link](https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/qj.3646))

Nipen, T.N., I.A. Seierstad, C. Lussana, J. Kristiansen, and Ø. Hov, 2020: Adopting Citizen Observations in
Operational Weather Prediction. Bull. Amer. Meteor. Soc., 101, E43–E57
([link](https://journals.ametsoc.org/doi/full/10.1175/BAMS-D-18-0237.1))

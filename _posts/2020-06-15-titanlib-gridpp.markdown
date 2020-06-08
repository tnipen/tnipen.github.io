---
layout: post
title:  "Surface analyses with titanlib and gridpp"
date:   2020-06-15 10:00:00 +0200
author: Thomas Nipen (thomasn@met.no)
tags: optimal_interpolation quality_control
---

MET Norway produces gridded surface analyses every hour, which are used to show current condtions on
MET Norway's weather site Yr ([https://www.yr.no](https://www.yr.no)); These analyses are created by
combining gridded output from a high resolution numerical weather prediction (NWP) model with observations
from in-situ weather stations, using optimal interpolation (OI). Currently 2 m temperature and hourly
precipitation are produced.

The analyses are produced using two open-source software packages, developed at MET Norway. **Titanlib** is
a quality control library that flags suspicious observations and **gridpp** is a post-processing library that
can, among other things, assimilate observations into a gridded background field from NWP models. Both were
designed around the need to integrate larger amounts of observations, with the rapid emergence of
crowdsourced weather observations, such as personal weather stations from [Netatmo](https://netatmo.com).

Both libraries are written in C++, but they also have python interfaces (and in the future R interfaces).
This post gives a step-by-step guide for creating surface analyses of temperature using titanlib and gridpp
in python.

## Getting started
First install [titanlib](https://github.com/metno/titanlib) and [gridpp](https://github.com/metno/gridpp),
which could be as easy as running:

{% highlight bash %}
pip3 install titanlib
pip3 install gridpp
{% endhighlight %}

... if you already have their non-python dependencies. If not, check out their respective webpages for
installation instructions. You will also need to download the following four files:

Type | Filename | Size
--- | --- |---
Observations | [obs.nc](https://thredds.met.no//thredds/fileServer/metusers/thomasn/gridpp/obs.nc) | 9.3 KB
Ensemble background | [analysis.nc](https://thredds.met.no//thredds/fileServer/metusers/thomasn/gridpp/analysis.nc) | 110 MB
Output template | [template.nc](https://thredds.met.no//thredds/fileServer/metusers/thomasn/gridpp/template.nc) | 45 MB
Python code | [surface_analysis.py]({{ site.url }}/assets/scripts/surface_analysis.py) | 4.3 KB

The python code is a collection of all code shown on this page.

## Quality control of observations

The first step is to ensure that the observations we use are trustworthy. Let's start by reading observations
and their station's metadata from file into arrays:

{% highlight python %}
import titanlib
import netCDF4
import numpy as np

with netCDF4.Dataset('obs.nc', 'r') as file:
    obs_lats = file.variables['latitude'][:]
    obs_lons = file.variables['longitude'][:]
    obs_elevs = file.variables['altitude'][:]
    obs = file.variables['air_temperature_2m'][:, 0]
{% endhighlight %}

### Spatial constistency test

Titanlib supports a variety of quality control methods, but we will focus on the spatial consistency test
(SCT) here. The SCT compares each observations to what is expected given the other observations in the nearby
area. If the deviation is large, the observation is removed.

The SCT has several parameters, which are described in detail on the [titanlib
wiki](https://github.com/metno/titanlib/wiki/Spatial-consistency-test). For a given observation, all other
observations within the `outer_radius` (in meters) will be used to compute a cross-validation value for the
observation. For computational efficiency, the cross-validation results will be reused for all observations
within the `inner_radius` such that the proceedure doesn't have to be performed independently for every
observation. To further reduce computation time, only the `num_max` closest observations are used, even if
there are more than this many observations within the outer radius. The test is only performed if there are at
least `num_min` observations within the outer radius.

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

flags, sct, rep = titanlib.sct(obs_lats, obs_lons, obs_elevs, obs,
        num_min, num_max, inner_radius, outer_radius, num_iterations,
        num_min_prof, dzmin, dhmin , dz, t2pos, t2neg, eps2)
{% endhighlight %}

The `flags` array now has the same length as the input observations and uses a value of `0` denoting
observations that passed the test and a value of `1` for those that are flagged. We will only keep the
non-flagged observations in the next step:

{% highlight python %}
index_valid_obs = np.where(flags == 0)[0]
index_invalid_obs = np.where(flags != 0)[0]
{% endhighlight %}

### Plotting the QC flags

Let's plot the observations, marking flagged ones with a black edge:

{% highlight python %}
import matplotlib.pylab as mpl
mpl.scatter(obs_lons[index_valid_obs], obs_lats[index_valid_obs],
        c=obs[index_valid_obs] - 273.15, s=50, linewidths=0, cmap="RdBu_r")
mpl.scatter(obs_lons[index_invalid_obs], obs_lats[index_invalid_obs],
        c=obs[index_invalid_obs] - 273.15, s=50, edgecolors="r", linewidths=1.5, cmap="RdBu_r")
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
observations, based on their relative accuracies.

Let's start by reading in the data from the NWP data file, and the 1 km output template. We will first use
the deterministic OI scheme, which takes a single ensemble member and for this we will take the control
member, which has index 0 in the background file.

{% highlight python %}
import gridpp
import numpy as np

with netCDF4.Dataset('analysis.nc', 'r') as file:
    index_control = 0
    blats_2500m = file.variables['latitude'][:]
    blons_2500m = file.variables['longitude'][:]
    belevs_2500m = file.variables['surface_geopotential'][0, 0, index_control, :, :] / 9.81
    bgrid_2500m = gridpp.Grid(blats_2500m, blons_2500m, belevs_2500m)
    background_2500m = file.variables['air_temperature_2m'][0, 0, index_control, :, :]

with netCDF4.Dataset('template.nc', 'r') as file:
    blats = file.variables['latitude'][:]
    blons = file.variables['longitude'][:]
    belevs = file.variables['altiude'][:]
    bgrid = gridpp.Grid(blats, blons, belevs)
{% endhighlight %}

`bgrid_2500m` and `bgrid` are objects that encapsulate the grid definitions.

### Downscaling to 1 km

The background has 2.5 km grid spacing. To downscale to 1 km, we do a height correction. Operationally, we
use a dynamic elevation gradient, but here we will use a constant -0.0065&deg;C/m elevation gradient. Gridpp
provides a method that does this downscaling, by taking the grid definitions, the original background field,
and the desired elevation gradient as arguments:

{% highlight python %}
gradient = -0.0065
background = gridpp.simple_gradient(bgrid_2500m, bgrid, background_2500m, gradient)
{% endhighlight %}

This yields a 1 km gridded field stored in `background`.

### Observation operator

OI also requires the background values at the observation points, also called the observation operator. In
many cases a simple method like nearest neighbour or bilinear is sufficient, however the background could also
be computed using any method. We will use gridpp's bilinear interpolator:

{% highlight python %}
points = gridpp.Points(obs_lats[index_valid_obs], obs_lons[index_valid_obs], obs_elevs[index_valid_obs])
pbackground = gridpp.bilinear(grid, points, background)
{% endhighlight %}

`gridpp.Points` creates an object that stores that encapsulates the point metadata (just like `gridpp.Grid`
does for a grid).

OI does not require you to specify the error variance of the background and the observations, only
their ratios is needed. Since we trust observations more than the background, we set the ratio to 0.1:

{% highlight python %}
variance_ratios = np.full(points.size(), 0.1)
{% endhighlight %}

`variance_ratios` is a vector, one value for each observation, which means different observations can be
assigned different variance ratios if there is reason to believe that the observations have different error
statistics.

### Specifying the structure function

Optimal interpolation needs a structure function. We will use a Barnes function (Barnes 1973), with a horizontal
decorrelation length of 100 km and a vertical decorrelation length of 200 m. This sets the limits to how far an
observation will affect the adjustment of the background.
{% highlight python %}
h = 100000
v = 200
structure = gridpp.BarnesStructure(h, v)
{% endhighlight %}

This structure function truncates the correlations past horizontal distances that are 3.64 greater than `h`.
Note, that this behvariour can be altered by specifing a truncation distance `hmax` in the structure
function.

Currently, gridpp also supports a Cressman structure function, which uses a linear decorrelation function
instead of an exponential. We aim to include further structure functions in the future.

### Running the OI

When OI processes a gridpoint, it needs to invert a matrix containing all the observations within the
influence radius. In areas where the observation network is dense, this matrix can be come large and will
result in large computation times. However, in most cases reducing the number of observations doesn't
negatively impact the accuracy of the analysis. We can set this using the `max_points` argument:

{% highlight python %}
max_points = 50
{% endhighlight %}

Note that this setting isn't important in this example with only 75 observations, but when including Netatmo
observations, this is critical.

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

Gridpp also supports an Ensemble-based Statistical Interpolation (EnSI; Lussana et al. 2019) scheme that
uses spatial structure information from an ensemble of NWP model runs. This is the method used in the
operational temperature analyses used on Yr, as described by Nipen et al. 2020. To test this, you need to
load the full ensemble and ensure that the ensemble dimension is at the end. Also note that EnSI requires the
standard error of the observation (`psigmas`), which we set to 0.5&deg;C. `psigmas` is a vector, one value for
each observation, which means different observations can be assigned different standard errors if there is
reason to believe that the observations have different error statistics.

{% highlight python %}
with netCDF4.Dataset('analysis.nc', 'r') as file:
    background_ens_2500m = file.variables['air_temperature_2m'][0, 0, :, :, :]
    background_ens_2500m = np.moveaxis(background_ens_2500m, 0, 2)

num_members = background_ens_2500m.shape[2]
pbackground_ens = np.zeros([points.size(), num_members])
for e in range(num_members):
    background_ens = gridpp.simple_gradient(bgrid_2500m, bgrid, background_ens_2500m[:, :, e], gradient)
    pbackground_ens[:, e] = gridpp.bilinear(bgrid, points, background_ens)

psigmas = np.full(points.size(), 0.5)

analysis_ens = gridpp.optimal_interpolation_ensi(bgrid, background_ens, points,
        obs[index_valid_obs], psigmas, pbackground_ens, structure, max_points)
diff = (analysis_ens - background_ens)[:, :, index_control]
{% endhighlight %}

![Analysis increment map for ensemble]({{ site.url }}/assets/img/analysis_increment_ens.png)

## Further work

We have serveral features we want to implement. First of all, the OI presented here assumes Gaussian error
characteristics. For non-Gaussian variables, such as wind and precipitation, we will support the
transformation of the variable.

We also hope to support other structure functions. Gridpp can support structure functions that can provide a
correlation function like this (in C++):

{% highlight c++ %}
float corr(const Point& p1, const Point& p2)
{% endhighlight %}

where `p1` and `p2` are points described by latitude, longitude, altitude, and land area fraction, and that
can use any arguments passed in the initialization of the class (such as `h`, `v`,`hmax` for `BarnesStructure`).

Finally, and improvement to the ensemble OI scheme presented here (EnSI) will be included that better handles
variables such as precipitation, where the ensemble variance often is 0, when all members have no
precipitation.

## References

Barnes, S. L., 1973: Mesoscale objective map analysis using weighted time-series observations. NOAA Tech.
Memo. ERL NSSL-62 [NTIS COM-73-10781], March, National Severe Storms Laboratory, Norman, Oklahoma, 60 pp.

Lussana, C, Seierstad, IA, Nipen, TN, Cantarello, L. Spatial interpolation of two‐metre temperature over
Norway based on the combination of numerical weather prediction ensembles and in situ observations. Q J R
Meteorol Soc. 2019; 145: 3626– 3643 ([link](https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/qj.3646))

Nipen, T.N., I.A. Seierstad, C. Lussana, J. Kristiansen, and Ø. Hov, 2020: Adopting Citizen Observations in
Operational Weather Prediction. Bull. Amer. Meteor. Soc., 101, E43–E57
([link](https://journals.ametsoc.org/doi/full/10.1175/BAMS-D-18-0237.1))

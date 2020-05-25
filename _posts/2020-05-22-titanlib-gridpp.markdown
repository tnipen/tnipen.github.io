---
layout: post
title:  "Surface analyses with titanlib and gridpp"
date:   2020-05-22 15:58:36 +0200
tags: optimal_interpolation quality_control
---

MET-Norway constructs a gridded analysis every hour based on an NWP background and all observations we can
get a hold on. This includes conventional weather stations, but also an emerging source of observations:
citizen weather stations. Specifically, we use observations from [Netatmo](https://netatmo.com)'s network of weather stations.

To help create gridded analyses, MET-Norway has built two software packages: titanlib and gridpp. Titanlib is
a quality control package that can flag suspicious observations and gridpp is a post-processing package that
can, among other things, assimilate observations into a gridded background field from NWP models. This
tutorial shows how to integrate observations with an NWP background field using python.

# Set up
First install [titanlib](https://github.com/metno/titanlib) and [gridpp](https://github.com/metno/gridpp):

{% highlight bash %}
pip3 install titanlib
pip3 install gridpp
{% endhighlight %}

You will also need two test files, an [observation file](https://thredds.met.no//thredds/fileServer/metusers/thomasn/gridpp/obs.nc) and a [gridded background file](https://thredds.met.no//thredds/fileServer/metusers/thomasn/gridpp/analysis.nc), which is a single timestep from an NWP ensemble run.


# Quality control of observations

The first step is to remove faulty observations. We will use the titan library for this. We will use the
spatial consistency test.

{% highlight python %}
import titanlib

# Load observations and metadata from file
with netCDF4.Dataset('obs.nc', 'r') as file:
    obs_lats = file.variables['latitude'][:]
    obs_lons = file.variables['longitude'][:]
    obs_elevs = file.variables['altitude'][:]
    obs = file.variables['air_temperature_2m'][:]

flags, sct, rep = titanlib.sct(obs_lats, obs_lons, obs_elevs, obs, minnumobs, maxnumobs, inner_radius, outer_radius, niterations, nminprof, dzmin, dhmin , dz, t2pos, t2neg, eps2)
{% endhighlight %}

The `flags` array has the same length as the observations and has a `0` denoting observations that passed the
QC and `1` for those that are flagged. We will only keep the non-flagged observations in the next step:

{% highlight python %}
index_valid_obs = np.where(flags == 0)[0]
obs_lats = obs_lats[index_valid_obs]
obs_lons = obs_lons[index_valid_obs]
obs_elevs = obs_elevs[index_valid_obs]
obs = obs[index_valid_obs]
{% endhighlight %}

# Surface analysis

This tutorial will only run the deterministic optimal interpolation scheme, which needs a single ensemble
member. We will take the control member, which has index 0 in the background file.

{% highlight python %}
import gridpp
import netCDF4
import numpy as np

with netCDF4.Dataset('analysis.nc', 'r') as file:
    index_control = 0
    blats = file.variables['latitude'][:]
    blons = file.variables['longitude'][:]
    belevs = file.variables['surface_geopotential'][0, 0, index_control, :, :] / 9.81
    bgrid = gridpp.Grid(blats, blons, belevs)
    background = file.variables['air_temperature_2m'][0, 0, index_control, :, :]

points = gridpp.points(obs_lats, obs_lons, obs_elevs)
variance_ratios = 0.1 * np.ones(points.size())
pbackground = gridpp.bilinear(grid, points, background)
{% endhighlight %}

Optimal interpolation needs a structure function. We will use a Barnes function, with a horizontal
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
        obs, variance_ratios, pbackground, structure, max_points)
{% endhighlight %}

# Plotting the result

Finally, we can plot the analysis increments:
{% highlight python %}
import matplotlib.pylab as mpl

diff = analysis - background
mpl.pcolormesh(blons, blats, diff, cmap="RdBu_r", vmin=-2, vmax=2)
mpl.plot(obs_lons, obs_lats, 'ko', mfc="w", ms=4)
mpl.xlim(0, 35)
mpl.ylim(55, 75)
mpl.gca().set_aspect(2)
mpl.show()
{% endhighlight %}

![Analysis increment map]({{ site.url }}/assets/img/analysis_increment.png)

# Ensemble mode

Gridpp also supports an Ensemble-based Statistical Interpolation (EnSI; Lussana et al 2019) scheme that takes
uses spatial structure information from an ensemble of NWP model runs. To test this, you need to load the
full ensemble:

{% highlight python %}
with netCDF4.Dataset('analysis.nc', 'r') as file:
    background_ens = np.moveaxis(file.variables['air_temperature_2m'][0, 0, :, :, :], 0, 2)

num_members = background_ens.shape[2]
pbackground_ens = np.zeros([points.size(), num_members])
for e in range(num_members):
    pbackground_ens[:, e] = gridpp.bilinear(bgrid, points, background_ens[:, :, e])

psigmas = 0.5 * np.ones(points.size())

analysis_ens = gridpp.optimal_interpolation_ensi(bgrid, background_ens,
        points, obs, psigmas, pbackground_ens, structure, max_points)
{% endhighlight %}

![Analysis increment map for ensemble]({{ site.url }}/assets/img/analysis_increment_ens.png)

# References

Lussana, C, Seierstad, IA, Nipen, TN, Cantarello, L. Spatial interpolation of two‐metre temperature over
Norway based on the combination of numerical weather prediction ensembles and in situ observations. Q J R
Meteorol Soc. 2019; 145: 3626– 3643 ([link](https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/qj.3646))

Nipen, T.N., I.A. Seierstad, C. Lussana, J. Kristiansen, and Ø. Hov, 2020: Adopting Citizen Observations in
Operational Weather Prediction. Bull. Amer. Meteor. Soc., 101, E43–E57
([link](https://journals.ametsoc.org/doi/full/10.1175/BAMS-D-18-0237.1))

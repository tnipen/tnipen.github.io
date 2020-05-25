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

Next, [download](https://thredds.met.no//thredds/fileServer/metusers/thomasn/gridpp/obs.nc) the observations file.

# Quality control

The first step is to remove faulty observations. We will use the titan library for this. We will use the
spatial consistency test.

{% highlight python %}
import titanlib

# Observations
with netCDF4.Dataset('obs.nc', 'r') as file:
    obs_lats = file.variables['latitude'][:]
    obs_lons = file.variables['longitude'][:]
    obs_elevs = file.variables['altitude'][:]
    obs = file.variables['air_temperature_2m'][:]

flags, sct, rep = titanlib.sct(obs_lats, obs_lons, obs_elevs, obs, minnumobs, maxnumobs, inner_radius, outer_radius, niterations, nminprof, dzmin, dhmin , dz, t2pos, t2neg, eps2)

index_valid_obs = np.where(flags == 0)[0]
obs_lats = obs_lats[index_valid_obs]
obs_lons = obs_lons[index_valid_obs]
obs_elevs = obs_elevs[index_valid_obs]
obs = obs[index_valid_obs]
{% endhighlight %}

# Surface analysis

For this you need NWP output. [Download](https://thredds.met.no//thredds/fileServer/metusers/thomasn/gridpp/analysis.nc) the ensemble analysis file.

{% highlight python %}
import gridpp
import netCDF4
import numpy as np

with netCDF4.Dataset('analysis.nc', 'r') as file:
    blats = file.variables['latitude'][:]
    blons = file.variables['longitude'][:]
    belevs = file.variables['surface_geopotential'][0, 0, 0, :, :] / 9.81
    bgrid = gridpp.Grid(blats, blons, belevs)
    background = file.variables['air_temperature_2m'][0, 0, 0, :, :]

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
output = gridpp.optimal_interpolation(bgrid, background, points, obs, variance_ratios, pbackground, structure, max_points)
{% endhighlight %}

## Plotting the result

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

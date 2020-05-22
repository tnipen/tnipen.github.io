---
layout: post
title:  "Surface analyses with titanlib and gridpp"
date:   2020-05-22 15:58:36 +0200
tags: optimal_interpolation quality_control
---

This tutorial shows how to integrate observations with an NWP background field using python. First install
[titanlib](https://github.com/metno/titanlib) and [gridpp](https://github.com/metno/gridpp):

{% highlight bash %}
pip3 install titanlib
pip3 install gridpp
{% endhighlight %}

# Quality control

Next, perform some quality control (QC) checks:

{% highlight python %}
import titanlib

titanlib.sct()
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
    bgrid = gridpp.Grid(lats, lons)
    background = np.moveaxis(np.squeeze(file.variables['air_temperature_2m'][:]), 0, 2)

points = gridpp.points(obs_lats, obs_lons, obs_elevs)
variance_ratios = 0.1 * np.ones(points.size())
pbackground = gridpp.bilinear(grid, points, background)
h = 10000
v = 200
structure = gridpp.BarnesStructure(h, v)
max_points = 10
output = gridpp.optimal_interpolation(grid, ivalues, points, obs, variance_ratios, pbackground, max_points)
{% endhighlight %}


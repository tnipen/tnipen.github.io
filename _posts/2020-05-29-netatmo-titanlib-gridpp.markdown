---
layout: post
title:  "Using Netatmo in weather forecasting"
date:   2020-05-29 15:58:36 +0200
author: Thomas Nipen (thomasn@met.no)
tags: netatmo optimal_interpolation quality_control
---

Networks of personal weather stations are growning immensely in size and now represent a major source of
weather information. MET-Norway has since 2018 used [Netatmo](https://netatmo.com)'s network of observations
in correcting temperature forecasts on Yr ([www.yr.no](https://www.yr.no)). Although the weather forecasts on
Yr are based on a sophisticated weather models, there are many cases where the model deviates significantly
from the observed weather. Netatmo observations are important to correct these errors. And since there is
bound to be a Netatmo station nearby, the forecasts can be corrected almost anywhere.

In this post, we will shown how to use temperature measurements from Netatmos network by using two open
source software packages: titanlib and gridpp. **Titanlib** is used to remove observations that are likely to be
invalid, whereas **gridpp** creates a gridded temperature field representing the current conditions. Both
libraries are written in C++, but they have python interfaces, which will be used here.

## Getting started
First, install [titanlib](https://github.com/metno/titanlib) and [gridpp](https://github.com/metno/gridpp),
which could be as easy as running:

{% highlight bash %}
pip3 install titanlib
pip3 install gridpp
{% endhighlight %}

... if you already have their non-python dependencies. If not, check out their respective webpages for
installation instructions.

For this tutorial, you also need to download [weather model output in NetCDF
format](https://thredds.met.no//thredds/fileServer/metusers/thomasn/gridpp/model.nc) and [Netatmo
observations in JSON format]({{ site.utl }}/assets/img/netatmo.json). The Netatmo observations were retrived
from their [open API](https://dev.netatmo.com/).

The python code below can also be [downloaded](TOOD).

## Quality control using Titanlib

Most of Netatmo's weather stations make reliable measurements, in fact, more than 80% are used in
MET-Norway's operational temperature forecasts. However, occasionally stations that are exposed to direct
sunlight or are located too close to buildings giving readings that are not representative of the area they
are in. These must be filtered out. Titanlib was developed specifically for dense networks of amateur weather
stations.

Let's start by reading the observations and their station's metadata from the downloaded file:

{% highlight python %}
import titanlib
import netCDF4
lats = lons = elevs = values = np.array(0)
with open(filename) as fid:
    text = fid.read()
    data = json.loads(text)
    for item in data["body"]:
        lon, lat = item["place"]["location"]
        elev = item["place"]["altitude"]
        for key, measure in item["measures"].items():
            if "type" in measure:
                types = measure["type"]
                if "temperature" in types:
                    unixtime = list(measure["res"].keys())[0]
                    I = types.index("temperature")
                    value = measure["res"][unixtime][I]
                    lats = np.append(lats, lat)
                    lons = np.append(lons, lon)
                    elevs = np.append(elevs, elev)
                    values = np.append(values, value)
{% endhighlight %}

The `obs` variable now contains the observations [in Kelvin] in a 1-D array, with `obs_lats`, `obs_lons`, and `obs_elevs`
being corresponding 1-D arrays of the stations' latitudes [degrees], longitudes [degrees], and elevations
[m], respectively.

Titanlib supports a variety of quality control methods, but we will use the **buddy check**. the **spatial
consistency test (SCT)** and the **isolation check**. The easiest way to run multiple checks in titanlib is
to create a dataset object:

{% highlight python %}
dataset = gridpp.Dataset(obs_lats, obs_lons, obs_elevs, obs)
{% endhighlight %}

### The buddy check

The buddy check can then we run on this dataset as follows:

{% highlight python %}
radius = 10000
min_buddies = 10
threshold = 2
dataset.buddy_check(radius, min_buddies, threshold)
{% endhighlight %}

This check will flag observations in the dataset that deviate by more than 2 standard deviations compared to
the average values of all observations within a 10 km radius. The test is only run if there are at least 10
stations within this radius.

# The Spatial Consistency Test

The next check we will apply is the SCT, which compares each observations to what is expected given the other
observations in the nearby area. If the deviation is large, the observation is removed. The SCT has several
parameters, which are described in detail on the [titanlib
wiki](https://github.com/metno/titanlib/wiki/Spatial-consistency-test).  For a given observation, all other
observations within the `outer_radius` (in meters) will be used to compute a cross-validation value for the
observation. For computational efficiency, the cross-validation results will be reused for all observations
within the `inner_radius` such that the proceedure doesn't have to be performed independently for every
observation. To further reduce computation time, only the `num_max` closest observations are used, even if
there are more than this many observations within the outer radius. The test is only performed if there are
at least `num_min` observations within the outer radius.

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

dataset.sct(num_min, num_max, inner_radius, outer_radius, num_iterations, num_min_prof, dzmin, dhmin, dz, t2pos, t2neg, eps2)
{% endhighlight %}

Finally, we use an isolation check to remove remote stations. We do this because when there are no nearby
stations, there is no information to corroborate the measurement and therefore we cannot guarantee that the
station has a valid measurement.

{% highlight python %}
radius = 15000
num min= 5
dz = 200
dataset.isolation_check(num_min, radius, dz)
{% endhighlight %}

The final information about which observations passed are available in `dataset.flags`. This array
has the same length as the input observations and uses `0` to denote observations that passed
the test and `1` for those that are flagged as suspicious. We will only keep the non-flagged observations in
the next step:

{% highlight python %}
index_valid_obs = np.where(flags == 0)[0]
index_invalid_obs = np.where(flags != 0)[0]
{% endhighlight %}

## Plotting the quality control flags

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

![Result of the quality control]({{ site.url }}/assets/img/netatmo_titan.png)

## Correcting the weather model

Once we have a set of trustworth observations, we can correct the output from the weather model. The model
file contains a gridded forecasts for the current conditions on a 2.5 km by 2.5 km grid. Gridpp contains,
among other things, functions for merging observations into gridded forecasts using optimal interpolation
(OI). This method merges the two data sources, based on their relative accuracies: The observations are
generally more accurate than the model output, and will therefore be weighted more heavily.

Let's continue by reading the gridded model data from file:

{% highlight python %}
import gridpp
import numpy as np

file = netCDF4.Dataset('model.nc', 'r')
blats = file.variables['latitude'][:]
blons = file.variables['longitude'][:]
belevs = file.variables['surface_geopotential'][0, :, :] / 9.81
bgrid = gridpp.Grid(blats, blons, belevs)
background = file.variables['air_temperature_2m'][0, :, :]
file.close()

points = gridpp.points(obs_lats[index_valid_obs], obs_lons[index_valid_obs], obs_elevs[index_valid_obs])
{% endhighlight %}

OI has a number of tuning parameters. The first one is the ratio of the error variance of the observations to
the model. Since we trust observations more than the background, we set the ratio to 0.1. OI
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

![Analysis increment map]({{ site.url }}/assets/img/netatmo_analysis.png)

## References

Barnes, S. L., 1973: Mesoscale objective map analysis using weighted time-series observations. NOAA Tech.
Memo. ERL NSSL-62 [NTIS COM-73-10781], March, National Severe Storms Laboratory, Norman, Oklahoma, 60 pp.

Lussana, C, Seierstad, IA, Nipen, TN, Cantarello, L. Spatial interpolation of two‐metre temperature over
Norway based on the combination of numerical weather prediction ensembles and in situ observations. Q J R
Meteorol Soc. 2019; 145: 3626– 3643 ([link](https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/qj.3646))

Nipen, T.N., I.A. Seierstad, C. Lussana, J. Kristiansen, and Ø. Hov, 2020: Adopting Citizen Observations in
Operational Weather Prediction. Bull. Amer. Meteor. Soc., 101, E43–E57
([link](https://journals.ametsoc.org/doi/full/10.1175/BAMS-D-18-0237.1))

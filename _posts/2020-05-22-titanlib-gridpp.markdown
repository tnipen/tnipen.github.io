---
layout: post
title:  "Surface analyses with titanlib and gridpp"
date:   2020-05-22 15:58:36 +0200
tags: optimal_interpolation quality_control
---

This tutorial shows how to integrate observations with an NWP background field using python. First install
[titanlib]: https://github.com/metno/titanlib and [gridpp]: https://github.com/metno/gridpp:

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

{% highlight python %}
import gridpp
{% endhighlight %}


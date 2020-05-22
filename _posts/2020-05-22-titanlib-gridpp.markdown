---
layout: post
title:  "Using titanlib and gridpp to produce surface analyses"
date:   2020-05-22 15:58:36 +0200
categories: jekyll update
---

This tutorial shows how to integrate observations with an NWP background field using python. First install
[titanlib][https://github.com/metno/titanlib] and [gridpp][https://github.com/metno/gridpp]:

{% highlight bash %}
pip3 install titanlib
pip3 install gridpp
{% endhighlight %}

Next, perform some quality control (QC) checks:

{% highlight python %}
import titanlib
import gridpp


titanlib.sct()
{% endhighlight %}

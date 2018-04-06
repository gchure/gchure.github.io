---
layout: post
title: Photography
cover_url: 'https://s3.amazonaws.com/gdc-photos/DSC_0482.jpg' 
permalink: photography
---

Sometimes I place myself behind a series of lenses for fun. These are some of
the results from those experiments. I occasionally post other work to my
[instagram page](https://www.instagram.com/gchure/) if you want to see more.

{% for photo in site.data.photos %}

<div style="text-align: left;">

<b>{{ photo.title }}</b> | {{photo.subtitle}}

<img style="width:auto; align: center; margin: auto; display: block; max-height:600px;" src="{{ site.data.settings.s3_url }}{{ photo.file }}">
<br />
</div>
{% endfor %}


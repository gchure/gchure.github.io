---
layout: post
cover_url: assets/beckman.jpg
title: CV
permalink: cv
---

[**Download PDF <i class="icon-pdf"></i>**]({{ site.baseurl}}/assets/GriffinChure_CV_201903.pdf)

## Education
**A.Sc. General Studies**, Utah State University, *2009*<br/>
**B.Sc. Biology** - Cell & Molecular Emphasis (Honors), University of Utah, *2013*<br/>
**B.Sc. Chemistry** - Biological Emphasis & Minor Physics, University of Utah, *2013*<br/>
**Ph.D. Biochemistry & Molecular Biophysics**, California Institute of Technology, *expected 2020* 
<div style="padding-left: 1em; margin-top:-2em;">
<i>Thesis topic:</i> Predictive versatility of the Monod-Wyman-Changeux model of Allostery in the context of transcriptional regulation.<br/> 
<i>Thesis adviser:</i> Professor Rob Phillips

</div>

## Publications
<sup>**☭** </sup> indicates equal contribution

{% for pub in site.data.pubs %}
1. {{ pub.authors}} *{{pub.title}}*, {{ pub.journal }} **{{ pub.volume }}** ({{pub.year}}). DOI: {{ pub.DOI }}. {% for link in pub.links %}   [**\[{{link[0]}}\]**]({{link[1]}}){% endfor %}{% if pub.other %}
<ul style='margin-top: -1.5em;'>
    <li style="margin-left:2em;list-style-type: circle;"> <i> {{pub.other}} </i></li></ul>
    {% endif %}
{% endfor %}
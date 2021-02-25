---
layout: page
title: Gallery
subtitle: From the featured folder
permalink: /gallery/
gallery_path: "assets/img/featured"
tags: [Gallery, Photo]
---

This is a photo gallery made from the static files in the `assets/img/featured` folder. 
I wanted to create automatically a simple gallery from a folder without having to create a markdown page as you would for the portfolio.

{% include gallery.html gallery_path=page.gallery_path %}

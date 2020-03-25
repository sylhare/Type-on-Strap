---
layout: post
title: Feature images
feature-img: "assets/img/pexels/desk-messy.jpeg"
thumbnail: "assets/img/thumbnails/desk-messy.jpeg"
image: "assets/img/thumbnails/desk-messy.jpg" #seo tag
tags: [Test, Lorem]
---

This is an example of a post which includes a feature image specified in the front matter of the post. 
The feature image spans the full-width of the page, and is shown with the title on permalink pages:

```yaml
feature-img: "assets/img/pexels/desk-messy.jpeg"
```

>  - And now it is working

You can also add images aligned in your post using:

{% include aligner.html images="pexels/book-glass.jpeg,pexels/desk-messy.jpeg" %}

```html
{% raw %}
{% include aligner.html images="pexels/book-glass.jpeg,triangle.png" %}
{% endraw %}
```

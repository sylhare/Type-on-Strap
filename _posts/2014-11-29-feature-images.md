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

You can also add images aligned in your post using the `aligner` include.
Make sure to separate all of the image path from in a string separated with `,`.
It by default look into `assets/img/` so give the path from there, example:

{% highlight ruby %}
{% raw %}
{% include aligner.html images="pexels/book-glass.jpeg,triangle.png" %}
{% endraw %}
{% endhighlight %}

{% include aligner.html images="pexels/book-glass.jpeg,pexels/desk-messy.jpeg" %}


Here you have two images side by side, but you can set more and set the amount per columns 
(by specifying the number of columns or let it be automatic using `"auto"`):

{% highlight ruby %}
{% raw %}
{% include aligner.html images="portfolio/cabin.png,portfolio/cake.png,portfolio/circus.png" column=3 %}
{% endraw %}
{% endhighlight %}

{% include aligner.html images="portfolio/cabin.png,portfolio/cake.png,portfolio/circus.png" column=3 %}

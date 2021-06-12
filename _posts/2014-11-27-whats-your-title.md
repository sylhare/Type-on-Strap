---
layout: post
title: What's your title
hide_title: true
feature-img: assets/img/feature-img/story.jpeg
author: mhagnumdw
tags: [Test, Image]
---

This is an example of a post which includes a feature image that has a
text, where you don't want to redisplay the title.
Mind your image size in order for the text to be displayed where you want it to.
The only limit is your imagination.

Here is how the yaml looks inside the post:

```yml
title: What's your title
hide_title: true
feature-img: assets/img/feature-img/story.jpeg
author: mhagnumdw
tags: [Test, Lorem]
```

You may wonder, why is there a title when you are not actually displaying it. <br>
Well that's due to some jekyll limitation:

> You **can't** set the **title** to the **empty string**

The title is used elsewhere than inside the post, for example in the blog page that list this post.
An empty title would break those pages and possibly prevents jekyll to render your blog. 

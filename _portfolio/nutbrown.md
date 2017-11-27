---
layout: post
title: Nutbrown.io
img: "assets/img/portfolio/nutbrown-io.png"
date: September 2014
tags: [Nutbrown, Stephen Nutbrown, Blog]
order: 6
---

![image]({{ site.baseurl }}/{{ page.img }})

This is my blog! You can read a bit more about the technical decisions for how it is hosted [here](https://nutbrown.io/2017/11/19/static-website.html).
Tech stack:
- S3 bucket for static content
- Codebuild for building static content from Jekyll
- Cloudfront distribution for performance and SSL.
- Route53 used for domain with Alias records to cloudfront
- [Github for the codebase](https://github.com/NutterzUK/nutbrown) (Feel free to make a PR!).


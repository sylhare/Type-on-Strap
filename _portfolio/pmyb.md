---
layout: post
title: Promote My Brand (PmyB)
feature-img: "assets/img/portfolio/submarine.png"
img: "assets/img/portfolio/pmyb.png"
date: 2017-09-03
order: 3
---

PromoteMyBrand (https://pmyb.co.uk) was deployed for a startup, who are now doing extremely well at influencer marketing. The business has employees familiar with WordPress. WordPress wouldn't be my first choice, but on this occasion it made sense as it is a blog-like page and the company wanted several writers to be able to create content for the site. It is hosted on an AWS EC2 instance and the static files are stored in S3 and distributed by Amazons CDN (Cloudfront) for performance. 

Tech Stack:
- RDS MySQL database; able to restore to the second if required.
- EC2 instance with automated backups for hosting the site
- Wordpress offload for each S3/Cloudfront integration: https://deliciousbrains.com/wp-offload-s3/


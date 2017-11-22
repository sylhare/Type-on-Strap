---
layout: post
title: Social Wage
feature-img: "assets/img/portfolio/social-wage.png"
img: "assets/img/portfolio/social-wage.png"
date: 27 September 2015
tags: [Lorem, Ipsum, portfolio]
order: 1
---

[Socialwage](https://socialwage.co.uk) is an influencer marketing agency. The idea is simple, influencers (mostly people with large social media accounts) sign up to be paid for promoting content. A sister site, pmyb.co.uk (the same company) allows companies to run advertising campaigns using the accounts who signed up through socialwage. Social wage allows users to log in and verify their various social media accounts to prove ownership of them and handles PayPal payments. 

This post considers the architecture and technology choices, then looks at some difficult lessons learnt. I keep this list as much for myself to allow some time to reflect and improve, as for others looking at my portfolio for possible recruitment.

Tech stack:
- Jetty webserver.
- Java 8
- Tapestry 5
- MySql database
- Maven
- Jenkins for builds
- Bitbucket for git.
- Hosted on an OVH cloud VPS.
- SDKs for Twitter, Facebook, PayPal etc.

When I look at this stack a few years later, I think back to the decisions behind each component. 

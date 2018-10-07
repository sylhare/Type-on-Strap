---
layout: post
title: AWS static website hosting.
tags: [AWS, Static website hosting, Jekyll, S3, Route53, Cloudfront]
---
I wanted to set up a website for posting tech-related content. My requirements were that it's secure, easy to edit, costs almost nothing to run and I don't need to worry about scaling. Additionally, i'd identified that it would be a bonus if it can be serverless, so I don't need to maintain/patch/update anything. This post goes through the decisions around infrastructure and software to set this up, and in particular why I avoided the the usual php frameworks (e.g Joomla, Wordpress) for my personal blog, instead going with Jekyll and AWS. The end result looks like this:
<br/>![image]({{ site.baseurl }}/assets/img/posts/aws-static-website.png)

Everything on this site is open source, so you can also browse the codebase [here](https://github.com/NutterzUK/nutbrown "Github project").

Choosing the right technology for your project is as important as the idea itself. I come from a Java background, but Java isn't the right tool for this usecase. It's important to avoid a particular language or environment becoming your '[golden hammer](https://www.exceptionnotfound.net/the-golden-hammer-anti-pattern-primers/ "golden hammer")', and I'm keen to discuss why. Let's say I went with Java and I picked a familiar framework like [Apache's Tapestry5](http://tapestry.apache.org/ "Apache Tapestry"). It would work, but it comes with two problems:
- Firstly, I now need a server running Java which I need to pay for, patch, backup and maintain.
- Secondly, when I want to re-deploy it, unless I start considering deployment strategies (blue/green etc) then it's going to require a little bit of downtime whilst it re-deploys. 

These limitations are unnecessary, given that the requirement is just to be able to serve HTML/CSS/Javascript to a user's browser so they can view my blog. Similarly, using Wordpress or Joomla then ties me in to needing a server which can process PHP and I will even end up in having to maintain a database! These are both overkill, why should a small blog require a database, or more generally, why should a small blog require any back-end processing at all? Joomla/Wordpress/Tapestry whilst familiar and useful, aren't a good decision for this simple use-case.

Based on this, I knew I wanted to produce static HTML/CSS/Javascript which needs hosting somewhere. What I didn't want to have to do is write raw HTML/CSS/Javascript and worry about various links breaking and adding to a navigation every time I add new pages. I set about looking at reviews for static website generators; there is actually a whole [website dedicated to static website generators](https://www.staticgen.com/ "Static website generators"). These frameworks generally give you some functionality you may expect from server-side frameworks, such as letting you write common layouts and components which you can then 'include', and handling navigation/lists of links, but the difference is that these are 'built' into just normal HTML/CSS/Javascript before you host it (it's not processed on a server, it is processed before you host it). Number 1 on that list, and with the most github stars is [Jekyll](https://jekyllrb.com/).

[Jekyll](https://jekyllrb.com/ "Jekyll"):
Jekyll is a simple, blog-aware, static site generator perfect for personal, project, or organization sites. Think of it like a file-based CMS, without all the complexity. Jekyll takes your content, renders Markdown and Liquid templates, and spits out a complete, static website ready to be served by Apache, Nginx or another web server.

This is perfect as it means I can produce static webpages easily and use Jekyll to handle the complexities of putting together those static files. It's also handy that there are lots of Jekyll themes available to choose from, saving a lot of the time and effort putting together a bespoke template (I have no requirement for my blog to look wildly different from your average website, as it's mainly informational). You can find several [free Jekyll templates here](http://jekyllthemes.org/ "Free Jekyll Templates"), and you can find the extremely familiar looking ['type-on-strap' template here](http://jekyllthemes.org/themes/Type-on-Strap/ "Type-On-Strap Jekll Template").

Another blogger has created a blog on [7 reasons NOT to use a static site generator](https://www.sitepoint.com/7-reasons-not-use-static-site-generator/). To summarise, the main points are:
1. You are on your own: You have to understand markdown and perhaps a version control system like git. As this is a technical blog and I'm a full time software engineer, using git is just a day-to-day norm.
2.  There are too many choices: I find this to be a moot point, particularly as that is true for CMS systems too, but also because I don't think choice is a bad thing.
3. The initial setup time: Yes, it takes some time to set up, but no longer than it does to set up wordpress. Even with one of the out of the box services which set up hosting etc with wordpress installed for you, you still need to go and find a theme and spend a while setting up navigation and theme properties.
4. No admin interface - This I find to also be a moot point, as any decent editor or even notepad++ will show you what you are producing whilst you type. Also, making/maintaining a website is not non-technical work (As much as wix etc would like to tell you otherwise). You will at some point have to write some code, or install some plugins, or debug why it's not doing what you expect.
5. Website consistency - this point talks about users putting scripts and undesired items on your site. This is bizarre; with a static site there is no notion of 'user accounts' for editing. Nobody else is putting anything into your page. This is also a moot point given the number of [wordpress exploits over the years](https://premium.wpmudev.org/blog/wordpress-security-exploits/).
6. Managing large sites - If it's a large site and it wants to scale, throwing unnecessary databases and PHP into the mix isn't the answer.
7. Sites with server-side functionality. Yes, if you want server-side functionality you may need to have something running server-side, but that can be replaced nowadays with the likes of AWS lambda functions which will only charge you for what you use. Furthermore, if you aren't using server-side functionality, having it there is just another thing to be attacked.

For my simple blog, it falls into the final comment Craig wrote, which is that CMS is often overkill. However, if I had a client or a non-technical person editing the site, I may want a nice dashboard for them to look at.

So, i've decided on Jekyll, i've decided on my template, and I need some hosting for static files. AWS will do this for free under the [free tier usage](https://aws.amazon.com/free/). I can store up to 5GB of free storage in Amazon S3 (Their storage solution) and they allow that to be made available via normal requests sent from a browser. In addition, if you use cloudfront you can secure your site with a free SSL certificate. Amazon CloudFront is a content delivery network (CDN). Content delivery networks provide a globally-distributed network of proxy servers which cache content, such as web videos or other bulky media, more locally to consumers, thus improving access speed for downloading the content. This means I can put my files in S3 and when a user requests the file, it'll go through a cloudfront 'edge' location (A data center near the user for lower latency). The edge location will cache that file and next time someone nearby requests it, it'll come from the edge location rather than S3. This gives performance benefits as well as allowing me to install an SSL certificate for free! So what does it cost?

Requirement   | Amazon Service | Price |
------------- | -------------- | ----- |
Domain registration   | Route53   | $39 per year   |
Hosting  | Amazon S3   | Free (Up to 5GB)   |
CDN  | CloudFront   | $0.08 per GB transferred ($0.08 per 680 page loads   |
SSL Certificates  | ACM   | Free   |
Automated deployment  | Codebuild   | Free (up to 100 minutes per month)   |

Total cost: Around $3.25 a month, equal to £2.46 at current rates.

What do I get for that ~£2.46 a month?
- My website hosted with up to 6 replicas of it (AWS promise 99.999999999% durability - the 11 9's they call it)
- My website cached at edge locations around the world.
- SSL certificates which can be expensive elsewhere.
- As much scaling as necessary; AWS will handle that without me touching anything.
- No need to patch/update/backup anything!

The hosting is practically free - it's the .io domain which costs a fair bit. It's £4.99 for [1&1's basic package](https://www.1and1.co.uk/web-hosting) which will run wordpress, but slowly, and certainly doesn't come with a CDN and unlimited scaling, and doesn't even include the .io domain which will cost another £39.99 on top (That's in GBP, it's the same number in USD at AWS). Including the domain, at 1&1 that comes to £8.52. The £2.46 at AWS is looking good! Of course if I get another few 50,000 requests per month it'll be the same price, but I'm fairly confident the basic package at 1&1 won't be up to handling that many requests, and i'm fairly confident I won't be getting that level of traffic.

All in all, if you are wanting to host a small blog and want to have SSL certificates, scale, have a custom domain, code versioning, and no worry of down time or maintenance, then this looks like the solution to go for.

I will post another blog post about how to go about setting up this infrastructure in the same way as soon as I get time.

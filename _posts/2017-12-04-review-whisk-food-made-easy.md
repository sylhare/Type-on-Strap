---
layout: post
title: Review - Whisk.com, really easy home-cooked food!
tags: [Cooking, Food, Review, Whisk.com]
published: true
---
I recently went to a [Coghlans School Indian cookery class]("https://www.cookingexpert.co.uk/collections/daytime-short-cookery-courses") and learnt something I wasn't expecting to. 
Cooking is really and doesn't have to take very long. For me though, the cooking course was easy as everyone had a step-by-step guide with all the right ingredients in front of them. I don't have this at home, as cooking requires planning, a task which should be, but isn't, at the front of my mind doing the weekly shop (Usually i'm thinking "Just gimme some food, I just want to go home and relax"). Whisk.com tries to fix this problem.

I did consider getting one of those [Hello-Fresh](https://www.hellofresh.co.uk) boxes delivered every week. This would make it easy, but I'd be paying £11.66 for a meal for two, a similar price to getting a take-away delivered directly to my door, which hardly seems worth it.

Then I found whisk.com. The site recommends meals based on meals you have previously rated and reviewed. You can then get it to add the ingredients for those meals to your basket at Tesco, Asda, Ocado or Waitrose.
All you have to do is pop on, see something you like the look of and add it to your basket. The ingredients then turn up at your house without costing any more than just buying them yourself from the supermarket. 

![image]({{ site.baseurl }}/assets/img/posts/whisk/whisk.png)

When it comes to cooking, you have the recipes and instructions all on the website, so you can just take your phone or tablet into the kitchen and follow along.

I bought the ingredients for a [really easy casserole](https://whisk.com/recipes/from-partners?url=http%3A%2F%2Fwww.deliciousmagazine.co.uk%2Frecipes%2Feasy-sausage-casserole%2F&viewSource=Shopping%20List) to give it a try. The ingredients turned up with my weekly shopping delivery and I had some fun cooking it. If you prefer to shop at a supermarket in person rather than online, Whisk will email you a shopping list to print off. I that I went for a very simple recipe, it was just to try it out, but everything went exactly to plan.

From a tech point of view there are a few pros and cons:

Pros:
- The site is really easy to use and stores dishes you've ordered so that when you get around to cooking you have all the recipes in one place.
- The website claims to use AI to recommend food dishes to you and to select the correct ingredients from the supermarket, we will see how good this is as I continue to use the site.
- You can review your basket to remove any duplicate items before the items are added to your basket (If you have recipes which both require the same ingredients you may end up with duplicates).
- The website gives you a really clear comparison of your basket price across the supermarkets it supports.

The big down-side is:
- As the supermarkets don't have an 'add-to-basket' API, Whisk relies on the end user inputting their username/password for the supermarket into Whisk.com's website. Whisk logs in on behalf of the user and adds the items, and hopefully doesn't use those credentials to buy themselves a few beers whilst there. Those credentials are saved by Whisk for future use, which means they must be stored in plain text (or encrypted with a key they hold to decrypt when they use it) - not good! From a security perspective, I was upset by this and can draw parallels to other screen scraping sites like mint.com which require usernames/passwords. Screen-scraping is an old and insecure method for accessing accounts. On this occasion, and for the first time I have ever done it, I shared my credentials with Whisk.com, mainly as I knew my supermarket (Asda) ask for additional security details whenever placing orders using my card, so I have some confidence that whisk can't buy themselves food using my account. 

Overall, I had a great experience with Whisk and i'll continue to use it. I'm looking forwards to trying out some more dishes and I may post the occasional recipe review.



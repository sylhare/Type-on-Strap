---
layout: post
title: Hack24 Winner - TrumpBot
feature-img: "assets/img/portfolio/trumpbot.png"
img: "assets/img/portfolio/trumpbot.png"
date: September 2014
order: 1
---

Hack24 was a fantastic event. I joined with a couple of colleagues, Lewis Foster and Tom Welch to take part in a challenge set by Accelerate Places. Our team name was of course TLS (Tom-Lewis-Steve, rather than transport layer security). The challenge was: "We want our hackers to think about helping people understand the ecosystem of news. How can we help people understand the broader eco-system and make more considered decision about how they consume and understand news? Make it inclusive, be neutral, tell a story and have fun!" 

If you have 5 mins and want to watch the winning video, here it is:
<iframe width="560" height="315" src="https://www.youtube.com/embed/px7ZnlLCVao" frameborder="0" allowfullscreen></iframe>

Trumbot was a hilarious bot with an attitude and some serious tech backing it up. Here is roughly how it worked:
- Trumbot follows you on twitter.
- In the background, a server on AWS collects data from Snopes.com and others, populating a database of known 'fake news'
- A Google chrome extension also lets the wider community feed in to this by marking website or pages as fake news. An algorithm weights websites based on the number of votes on each site.
- You post on twitter. Trumbot is watching...
- It pulls the text from your post. If the post contains a link to a 'known-bad' URL (From the community), or it has similar words to any of the fake news stories (Measured using Lesk's algorithm (1985 and WordNet similarity for Java), it'll know it's fake.
- If it's fake... Trumbot posts back some potentially horrendous abuse in response to your fakenews post. Of course, the world wouldn't be complete without a whole API dedicated to delivering Donald Trump based insults [this can be found here!](http://2017.compciv.org/syllabus/assignments/homework/serials/trump-tweets-json.html) The really funny thing was that 'potentially horrendous' was indeed pretty bad, Donald has said some bad things in his time. We had to actually implement a filter on the API responses to avoid the amount of racism and various death threats that came from it. However, in the end, we had a toned down version of trump delivering slightly less horrendous insults. 

It ended up looking like this:
![image]({{ site.baseurl }}/assets/img/portfolio/trumpbot/trumpbot.png)

Here's the hack24 page showing the winning team! [https://www.hack24.co.uk/2017-winners/](https://www.hack24.co.uk/2017-winners/)

It really was great to work with such a fantastic team and it was a really well organised event; thanks to everyone involved!

The prize is also awesome; we won crystal maze tickets each and a [jumping parrot minidrone](http://www.argos.co.uk/product/3791550)! 

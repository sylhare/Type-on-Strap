---
layout: post
title: "In Favor of Developing Ethical Best Practices in AI Research"
short-summary: "What are best practices AI researchers can follow to avoid unintended consequences of their research?"
summary: "What are the ethical responsibilities of AI researchers? Or to put it in more pragmatic terms, what are best practices AI
researchers can follow to avoid unintended consequences of their research?"
thumbnail: "assets/img/posts/2019-02-21-ethical_best_practices/icon.png"
opinion: true
hide: true
author: <a href='http://web.stanford.edu/~shushman/'>Shushman Choudhury</a>, <a href='https://twitter.com/michellearning'>Michelle Lee</a>, and <a href='https://twitter.com/andrey_kurenkov'>Andrey Kurenkov</a>
tags: [opinion, ethics]
---

*Disclaimer: this is an opinion piece that represents the views of its authors, and not all of SAIL. We would like to thank SAIL professors Fei-Fei Li, Chris Manning, and Dorsa Sadigh for reviewing and approving this post. We would also like to thank Margaret Mitchell, Timnit Gebru, Joy Buolamwini, Rachel Thomas, and many others for inspiring us to write this.*

What are the ethical responsibilities of AI researchers? Or to put it in more pragmatic terms, what are best practices AI
researchers can follow to avoid unintended consequences of their research?

<figure>
    <img class="postimagehalf" style="width:69%;" src="{{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image13.png"/> 
    <img class="postimagehalf" style="width:30%;" src="{{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image2.png"/> 
    <figcaption> 
      Left: Microsoft’s <a href="https://qz.com/653084/microsofts-disastrous-tay-experiment-shows-the-hidden-dangers-of-ai/&sa=D&ust=1550739471096000">infamous chatbot</a> Tay (now discontinued). Right: https://xkcd.com/1390/
    </figcaption>
</figure>
<br>
Despite the meteoric rise of AI research over the
past decade, our research community still does not
regularly and openly discuss questions of ethics and responsibility.
Every researcher learns a set of best practices for doing
impactful research -- research published at conferences and
in journals -- but not all of us are asked or expected to
develop best practices for ethical research that prevents potential
misuse. This is already or increasingly the norm for other influential
disciplines such as [law](https://engagedscholarship.csuohio.edu/cgi/viewcontent.cgi?article%3D2157%26context%3Dclevstlrev&sa=D&ust=1550739471097000),
[medicine](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2598142/&sa=D&ust=1550739471097000),
and [engineering](http://theinstitute.ieee.org/ieee-roundup/blogs/blog/why-schools-are-getting-more-serious-about-teaching-engineering-students-about-ethics&sa=D&ust=1550739471098000), and we believe it should be normalized in the education and
research of AI as well.

{% figure caption:'Evidence over growing concerns about ethics and AI. Source: [cbiinsights.com](https://www.cbinsights.com/research/google-amazon-facebook-apple-hiring-techlash/)' %}
[<img class="postimage_75" src="{{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image16.png"/>]({{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image16.png)
{% endfigure %}

Yes, the topic has received growing attention, particularly within
the subfield of research focused on Fairness,
Accountability, and Transparency for AI (see  [FAT ML](http://www.fatml.org/&sa=D&ust=1550739471099000) and
[FAT*](https://fatconference.org/&sa=D&ust=1550739471099000)).
However, as PhD students working in AI, we have never been explicitly
instructed or even incentivized to speak openly and enthusiastically
about this topic. This, despite a documented
growing concern about AI and ethics in the general public (see
above).

This needs to change. Academic AI
researchers are now routinely having an impact in industry,
journalists are now increasingly turning to researchers for quotes, and
people in general can agree that AI as a field has never had as much
influence on society as it does today. At the very least,
all AI researchers and engineers should be aware of the
sorts of ethical hypotheticals and contingencies they may encounter in
their work and how to respond to them. 

In this piece, we intend to promote several best practices that we
think AI researchers and engineers should be aware of. We focus
primarily on what researchers can do to avoid unintended negative
consequences of their work, and do not go into depth what those negative
consequences can be or how to deal with intentionally bad actors. The
ideas we promote are not new (in fact, most of them are ideas suggested
by prominent researchers who we shall credit), nor are they
comprehensive, but they are practices we would like to highlight as a
beginning to a larger discussion on this topic. Our hope is to promote
awareness of these concepts and to inspire other researchers to join
this discussion.

### Education 
If you have read this far into the article, you
hopefully agree that it is reasonable to expect researchers to think
about the broader implications of their research. However, perhaps you
do not know where to begin to embark on this seemingly daunting task.
Fortunately, the first step you can take, as all of us should who care
about this issue, is rather straightforward - to become more informed
about the ethical concerns of AI, at least in your subfield.

**Practice: Familiarize yourself with the basics of AI ethics**

The legal and policy communities have thought about the
concerns of AI as extensively as the technical community has thought
about its development. Even a cursory search on the web will yield some
thought-provoking and well-researched works:


-   Machine learning researchers and practitioners may find it
    insightful to ponder on the [opacity of machine learning algorithms](http://journals.sagepub.com/doi/abs/10.1177/2053951715622512&sa=D&ust=1550739471102000) (even more relevant in the present age of learning with
    deep neural networks).
    
<figure>
    <img class="postimagehalf" src="{{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image5.png"/> 
    <img class="postimagehalf" src="{{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image6.png"/> 
    <figcaption> 
      Left: Redlining (systemic denial of services to residents of certain districts, often racially defined) in the 1940s. Right: Worrying patterns replicated by algorithms (<a href="https://www.skynettoday.com/editorials/biased-ai&sa=D&ust=1550739471102000">source</a>)
    </figcaption>
</figure>

-   For those working in big data analysis, many of the ethical
    issues, particularly stemming from disparate impact, have been
    [outlined in a recent paper](https://papers.ssrn.com/sol3/papers.cfm?abstract_id%3D2477899&sa=D&ust=1550739471103000) in
    Big Data & Society as well as the report “[Unfairness By Algorithm: Distilling the Harms of Automated Decision-Making](https://fpf.org/2017/12/11/unfairness-by-algorithm-distilling-the-harms-of-automated-decision-making/&sa=D&ust=1550739471103000)”
    by the future of privacy forum.
-   Researchers in robotics, i.e. embodied AI systems, can start
    with a discussion on the issues of [robot ethics in a mechanized world](https://www.sciencedirect.com/science/article/pii/S0004370211000178&sa=D&ust=1550739471104000). 

{% figure caption:'The biased results found by ACLU - people of color who are members of congress were found to be disproportionately incorrectly classified as being criminals 39% of the time with facial recognition technology from Amazon. [(source)](https://www.aclu.org/blog/privacy-technology/surveillance-technologies/amazons-face-recognition-falsely-matched-28&sa=D&ust=1550739471104000) ' %}
[<img class="postimage" src="{{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image3.png"/>]({{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image3.png)
{% endfigure %}

-   Facial recognition systems, a
    prominent application of computer vision, have a host of related
    [ethical concerns](https://www.emeraldinsight.com/doi/pdf/10.1108/14779960480000246&sa=D&ust=1550739471105000).
    The issue of [emotional privacy](http://blog.practicalethics.ox.ac.uk/2014/03/computer-vision-and-emotional-privacy/&sa=D&ust=1550739471105000) while decoding facial pain expressions is also relevant. This [broader survey](https://jyx.jyu.fi/bitstream/handle/123456789/55806/1/URN%253ANBN%253Afi%253Ajyu-201711084167.pdf&sa=D&ust=1550739471105000) highlights a range of additional ethical challenges for computer vision.

{% figure %}
[<img class="postimage" src="{{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image10.png"/>]({{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image10.png)
{% endfigure %}

-   For those who think about Natural Language Processing, it
    would be instructive to become familiar with some of the the
    significant [social impacts](https://aclanthology.coli.uni-saarland.de/papers/P16-2096/p16-2096&sa=D&ust=1550739471106000) of
    NLP, such as demographic misrepresentation, reinforcing linguistic
    biases, and topic overexposure and underexposure. Additionally, there has already been a paper titled
    [“Ethical by Design: Ethics Best Practices for Natural Language Processing”](http://aclweb.org/anthology/W17-1604.pdf) that is 
    of course relevant to the topic of the present article.
-   More broadly, many AI algorithms and developments are
    undeniably ‘[dual-use](https://en.wikipedia.org/wiki/Dual-use_technology&sa=D&ust=1550739471106000)’
    technologies (technologies which are designed for civilian purposes
    but which may have military applications, or more broadly designed
    for certain beneficial uses but can be abused for negative impacts).
    This concept is far from new and discussion over how to deal with it
    has long been established in fields such as [software security](https://dl.acm.org/citation.cfm?id%3D2622674.2622860&sa=D&ust=1550739471107000), and given its relevance to AI it is important that we
    are aware of it as well. 

These are just a few useful starting points to demonstrate that AI
researchers can (and we think, should) actively educate themselves on
these topics. For a more comprehensive overview, [Eirini Malliaraki](https://medium.com/@eirinimalliaraki/toward-ethical-transparent-and-fair-ai-ml-a-critical-reading-list-d950e70a70ea&sa=D&ust=1550739471107000) has
compiled [a list of books and articles](https://medium.com/@eirinimalliaraki/toward-ethical-transparent-and-fair-ai-ml-a-critical-reading-list-d950e70a70ea&sa=D&ust=1550739471108000) to
get you up to speed on many of the relevant subjects (such as
algorithmic transparency, data bias, and the social impact of
AI). The list is long and may appear
intimidating, but we recommend starting with a few topics that are
directly or even indirectly related to your own research and start
adding them to your reading list. A number of
recent papers have also done fantastic reviews of a large amount of
information, and are longer reads worth looking at:

-   [AI4People—An Ethical Framework for a Good AI Society: Opportunities, Risks, Principles, and Recommendations](https://link.springer.com/article/10.1007%252Fs11023-018-9482-5&sa=D&ust=1550739471108000)
-   [50 Years of Test (Un)fairness: Lessons for Machine Learning](https://arxiv.org/abs/1811.10104&sa=D&ust=1550739471109000)
-   [The Malicious Use of Artificial Intelligence: Forecasting, Prevention, and Mitigation](https://arxiv.org/pdf/1802.07228.pdf&sa=D&ust=1550739471109000)
-   [Perspectives on Issues in AI Governance](https://ai.google/static/documents/perspectives-on-issues-in-ai-governance.pdf&sa=D&ust=1550739471110000)
-   [Thinking About Risks From AI: Accidents, Misuse and Structure](https://www.lawfareblog.com/thinking-about-risks-ai-accidents-misuse-and-structure&sa=D&ust=1550739471110000)

In addition to individual research works, there is a growing
amount of institutional interest. [Rachel Thomas](https://www.fast.ai/2018/09/24/ai-ethics-resources/&sa=D&ust=1550739471111000),
herself a prominent AI and AI ethics researcher, has compiled a list of
researchers and research institutes working on fairness and ethics in
AI. If you are at a university, or live close to one, it
might be worthwhile to take or audit a course on ethics in AI.
[Casey Fiesler](https://medium.com/@cfiesler/tech-ethics-curricula-a-collection-of-syllabi-3eedfb76be18&sa=D&ust=1550739471111000) has
crowdsourced over 200 courses on technology
and ethics in universities all around the world, along with
their class syllabi for reference. Besides classes, there
are also easy-to-digest online compilations of information, such as
[Ethics in Technology Practice](https://www.scu.edu/ethics-in-technology-practice/&sa=D&ust=1550739471111000) from
the Markkula Center for Applied Ethics at Santa Clara University.

<figure>
<iframe width="560" height="315"
src="https://www.youtube.com/embed/59bMh59JQDo" frameborder="0"
allow="accelerometer; autoplay; encrypted-media; gyroscope;
picture-in-picture" allowfullscreen></iframe>
</figure>

**Practice: Codes of ethics and pledges**

{% figure caption:'[Source](https://www.pinterest.com/pin/184225440982670733/?autologin%3Dtrue&sa=D&ust=1550739471112000)' %}
[<img class="postimage_50" src="{{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image8.png"/>]({{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image8.png)
{% endfigure %}

The vast amount of information online can sometimes feel unwieldy
to the uninitiated. Luckily, AI is far from the first
context in which academics interested in knowledge have had to deal with
ethical questions. Nor is the need to think about ethical questions new
to AI itself. Therefore, a number of distilled codes of ethics exist
that concisely summarize the key points one should keep in mind:


-   CS professionals

Though not specific to AI, the codes of ethics of
both [IEEE](https://www.ieee.org/about/corporate/governance/p7-8.html&sa=D&ust=1550739471113000) and
[ACM](https://www.acm.org/code-of-ethics&sa=D&ust=1550739471113000) are quick to read and entirely relevant. Many of the
principles in these codes, such as being honest or not taking bribes,
represent common sense. But, whenever you may be in doubt as to the
ethical nature of possible research actions, it is a good idea to review
them.

{% figure caption:'Attendees at the [Beneficial AI Conference 2017](https://futureoflife.org/bai-2017/), which led to the Asilomar AI Principles. ' %}
[<img class="postimage" src="{{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image4.png"/>]({{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image4.png)
{% endfigure %}

-   AI Researchers

Of course, academic research in AI has a set of issues and
concerns unique to it that these general codes of ethics for technology
and computing professional may not address. Fortunately, substantial
effort has been put into addressing this area as well. In particular,
["ETHICALLY ALIGNED DESIGN - a Vision for Prioritizing Human Well-being with Autonomous and Intelligent Systems"](https://ethicsinaction.ieee.org/&sa=D&ust=1550739471114000) by
the IEEE Global Initiative on Ethics of Autonomous and Intelligent
Systems and the [ASILOMAR AI Principles](https://futureoflife.org/ai-principles/?cn-reloaded%3D1&sa=D&ust=1550739471114000) by
the Future of Life Institute both provide concrete recommendations for
AI researchers. Last but not least, researchers should be familiar with
the codes of conducts at their universities and professional events they
attend. Not all AI conferences have explicit codes of conduct, so a good
baseline to be aware of is the [NeurIPS 2018 code of conduct](https://neurips.cc/public/CodeOfConduct&sa=D&ust=1550739471114000).
Likewise, the NeurIPS affiliated [ML Ally pledge](https://sites.google.com/view/ml-ally-pledge/pledge&sa=D&ust=1550739471115000) is worth reading and thinking about.


-   AI Influencers

    Many researchers may also have the potential to have an impact
    beyond academia, such as in policy or industry. In addition to the
    prior recommendations, the [Montréal Declaration for Responsible AI](https://www.montrealdeclaration-responsibleai.com/&sa=D&ust=1550739471115000) and
    [AI4People list of Principles and Recommendations](https://link.springer.com/article/10.1007%252Fs11023-018-9482-5&sa=D&ust=1550739471115000) offer
    an excellent overview of things to consider when developing AI
    in general. And, there are also the more specific [Lethal Autonomous Weapons Pledge](https://futureoflife.org/lethal-autonomous-weapons-pledge/&sa=D&ust=1550739471115000),the [Campaign to Stop Killer Robots](https://www.stopkillerrobots.org/&sa=D&ust=1550739471116000), and
    the [Safe Face Pledge](https://www.safefacepledge.org/&sa=D&ust=1550739471116000),
    which are likewise very relevant to anyone involved in the research
    and development of AI technology. Specific companies and research
    labs (such as [Google](https://www.blog.google/technology/ai/ai-principles/&sa=D&ust=1550739471116000),
    [OpenAI](https://blog.openai.com/openai-charter/&sa=D&ust=1550739471116000),
    [DeepMind](https://deepmind.com/applied/deepmind-ethics-society/principles/&sa=D&ust=1550739471116000),
    [Salesforce](https://www.salesforce.org/ai-good-principles-believe/&sa=D&ust=1550739471117000),
    and [IBM](https://www.ibm.com/blogs/policy/trust-principles/&sa=D&ust=1550739471117000))
    have also begun to specify their principles, and when considering
    joining these institutions it is worth reviewing these documents.
    
{% figure caption:'[Source](https://motherboard.vice.com/en_us/article/vvv559/the-campaign-to-stop-killer-robots-is-not-going-well&sa=D&ust=1550739471117000)' %}
[<img class="postimage" src="{{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image14.png"/>]({{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image14.png)
{% endfigure %}

### Communication & Distribution

As far as we have seen, the potential misuses and ethical
considerations of new AI algorithms and products are rarely identified
and pointed out in documentation or academic
papers. A prominent practical example is that
Amazon’s documentation for its Rekognition product [did not have warnings on changing the default parameters of the product](https://www.skynettoday.com/briefs/aclu-amazon-rekognition&sa=D&ust=1550739471118000) for law enforcement use cases until after the ACLU pointed
out that the product could be misused to classify US senators as
criminals. 

Perhaps even more importantly, researchers do not just communicate
ideas with their papers -- they also distribute code, data, and models
to the wider AI society. As the capabilities of AI systems continue to
become stronger, considerations of [dual-use](https://en.wikipedia.org/wiki/Dual-use_technology&sa=D&ust=1550739471118000) will
have to prompt us to develop a set of new best practices with regards to
distribution, some of which we discuss here.

**Practice: Ethical Consideration Sections**

A novel and impactful practice researchers can undertake
now is to include a section on Ethical Considerations in our
papers, something that machine learning researchers in the Fairness, Accountability and Transparency
sub-community have already started to do. A prominent example is [Margaret Mitchell](http://www.m-mitchell.com/&sa=D&ust=1550739471119000),
a Senior Research Scientist at Google AI, and Tech Lead of Google’s ML
fairness effort, who has included an ethical consideration
section in several of her recent papers. For instance,
[her 2017 paper](http://www.m-mitchell.com/publications/multitask-clinical.pdf&sa=D&ust=1550739471119000) on
predicting imminent suicide risk in a clinical care scenario using
patient writings, flagged the potential for abuse of the research by
singling out people, which the authors addressed by anonymizing the
data. Clearly, this is particularly relevant for research with potential
for [dual-use](https://en.wikipedia.org/wiki/Dual-use_technology&sa=D&ust=1550739471120000).
Her [blog post](http://www.m-mitchell.com/publications/multitask-blurb.html&sa=D&ust=1550739471120000) provides
even more details. 

{% figure caption:'An example of an ethical considerations section (from [Margaret Mitchell's blog post](http://www.m-mitchell.com/publications/multitask-blurb.html&sa=D&ust=1550739471120000))' %}
[<img class="postimage_50" src="{{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image12.png"/>]({{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image12.png)
{% endfigure %}

**Practice: Cards, Certificates, and Declarations**

{% figure caption:'A ‘Dataset nutrition label’ from [the project website](http://datanutrition.media.mit.edu/&sa=D&ust=1550739471121000)' %}
[<img class="postimage" src="{{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image17.png"/>]({{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image17.png)
{% endfigure %}

Recently, groups from Google and IBM
Research have proposed standardized means of
communicating aspects of new datasets and AI services
in the papers [Model Cards for Model Reporting](https://arxiv.org/abs/1810.03993&sa=D&ust=1550739471122000),
[Datasheets for Datasets,](https://arxiv.org/abs/1803.09010&sa=D&ust=1550739471122000) [Data Statements for
NLP](https://openreview.net/forum?id%3DBy4oPeX9f&sa=D&ust=1550739471122000),
[The Dataset Nutrition Label](http://datanutrition.media.mit.edu/&sa=D&ust=1550739471122000), [Policy Certificates: Towards Accountable Reinforcement Learning](https://arxiv.org/pdf/1811.03056.pdf&sa=D&ust=1550739471122000),
and [Increasing Trust in AI Services through Supplier's Declarations of Conformity](https://arxiv.org/abs/1808.07261&sa=D&ust=1550739471123000).
These methods allow researchers to communicate important
information about their work such as a model’s use cases, a dataset’s
potential biases, or an algorithm’s security considerations.
As the impact of AI research on society continues to grow,
we should consider adopting these new standards
of communication. 

**Practice: Approval and Terms of Access for Datasets, Code, and Models**

ImageNet has been among the most important datasets
in Computer Vision, but what many may not be aware is that being given
easy download access requires [going through an approval stage and agreeing with precise terms of access](http://image-net.org/download-faq&sa=D&ust=1550739471124000).
It’s far from the only dataset that mandates a
request before being shared, with a newer example being
[The Pilot Parliaments Benchmark](https://www.ajlunited.org/gender-shades&sa=D&ust=1550739471124000) (see
above). The benefit of such a process is clear for any dataset with a
potential for [dual-use](https://en.wikipedia.org/wiki/Dual-use_technology&sa=D&ust=1550739471124000), although it is admittedly not without some overhead for the
lab or organization distributing the dataset. The same process could
also be applied for code and pretrained models, which of course also
have potential for dual-use, though this precedent there is not as
established; in general, we believe the AI research community will need
to discuss and develop new best practices for distribution of data,
code, and models that are essential for reproducibility but may be put
to harmful use. 

{% figure caption:'A set of considerations related to distributing research results
Google highlighted in [Perspectives on Issues in AI Governance](https://ai.google/static/documents/perspectives-on-issues-in-ai-governance.pdf&sa=D&ust=1550739471125000)' %}
[<img class="postimage" src="{{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image11.png"/>]({{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image11.png)
{% endfigure %}

**Practice: Revise Peer Review**

A reasonable retort to the above suggestions might be that they
are not typically done today, and the effort needed to follow them may
not help or even hurt your paper’s chances of getting accepted. That is
why we endorse the position put forth in “[It’s Time to Do Something: Mitigating the Negative Impacts of Computing Through a Change to the Peer Review Process](https://acm-fca.org/2018/03/29/negativeimpacts/&sa=D&ust=1550739471126000)”.
 As summarized well in a [New York Times](https://www.nytimes.com/2018/10/22/business/efforts-to-acknowledge-the-risks-of-new-ai-technology.html&sa=D&ust=1550739471126000) opinion
piece, “46 academics and other
researchers, are urging the research community to rethink the way it
shares new technology. When publishing new research, they say,
scientists should explain how it could affect society in negative ways
as well as positive.”

**Practice: Use, share, and create emerging tools and datasets**

<figure>
<iframe width="560" height="315"
src="https://www.youtube.com/embed/hpYl8WLYeKo" frameborder="0"
allow="accelerometer; autoplay; encrypted-media; gyroscope;
picture-in-picture" allowfullscreen></iframe>
</figure>

Lastly, there are several new and emerging tools and
datasets that can help you determine if your models or your dataset have
unintended biases and so check for that prior to wider distribution.
Some that we would like to highlight are:

-   [AI Fairness 360](https://www.ibm.com/blogs/research/2018/09/ai-fairness-360/&sa=D&ust=1550739471128000), by IBM: an open-source toolkit of different metrics and
    algorithms developed from the broader Fairness AI community that
    checks for unfairness and biases in models and datasets
-   [Facets](https://pair-code.github.io/facets/&sa=D&ust=1550739471128000),
    by [Google’s AI + People Research
    group](https://ai.google/research/teams/brain/pair&sa=D&ust=1550739471128000): two robust visualization tools to aid in understanding
    and analyzing machine learning datasets, which might be helpful in
    quickly identifying biases in your datasets 
-   [What If](https://pair-code.github.io/what-if-tool/&sa=D&ust=1550739471129000),
    by [Google’s AI + People Research
    group](https://ai.google/research/teams/brain/pair&sa=D&ust=1550739471129000)::
    a neat tool to play “what if” with theories of fairness, see the
    trade-offs, and make the difficult decisions that only humans
    can make.
-   [gn_glove](https://github.com/uclanlp/gn_glove&sa=D&ust=1550739471129000),
    by the authors of [Learning Gender-Neutral Word Embeddings (EMNLP 2018)](https://arxiv.org/abs/1809.01496&sa=D&ust=1550739471129000): a set of gender-neutral word vectors meant to remove
    stereotypes that have been shown to exist in prior word vectors.
   
<figure>
    <img class="postimagehalf" src="{{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image9.png"/> 
    <img class="postimagehalf" src="{{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image15.png"/> 
    <figcaption> 
      Google’s What If tool can be used to visualize inference results (left) and test algorithmic fairness (right)
    </figcaption>
</figure>

-   [WinoBias](https://github.com/uclanlp/corefBias&sa=D&ust=1550739471130000),
    by the authors of [Gender Bias in Coreference Resolution: Evaluation and Debiasing Methods](https://arxiv.org/abs/1804.06876&sa=D&ust=1550739471130000): a benchmark for coreference resolution focused on
    gender bias.
-   [Wikipedia Toxicity Dataset](https://meta.wikimedia.org/wiki/Research:Detox/Data_Release&sa=D&ust=1550739471131000), by Wikimedia: an annotated dataset of 1m crowd-sourced
    annotations targeting personal attacks, aggression,
    and toxicity.
-   [The Pilot Parliaments Benchmark (PPB)](https://www.ajlunited.org/gender-shades&sa=D&ust=1550739471131000), by the Algorithmic Justice League: a dataset of human
    faces meant to achieve better intersectional representation on the
    basis of gender and skin type. It consists of 1,270 individuals,
    selected for gender parity, in the national parliaments of three
    African countries and three European countries..
-   [The diverse facial recognition dataset](https://www.research.ibm.com/artificial-intelligence/trusted-ai/diversity-in-faces/&sa=D&ust=1550739471131000), by IBM: a dataset of 36,000 images that is equally
    distributed across skin tones, genders, and ages.

These are examples we have been able to find, but in
general keeping an eye out for such datasets and tools and considering
them for your own research is a sensible idea.

### Advocacy 

Fairness and ethical AI is a growing field in and of
itself. If you would like to go beyond the basic expected ethical
practices of educating yourself and communicating potential misuses of
your creations, here are suggestions on how you can help make the field
of AI and the tools our peers create to become more ethical, inclusive,
and fair. 

**Practice: Bring up Concerns in Teaching and Talks**

{% figure caption:'Image from one of Stanford AI Lab’s ‘[AI Salon](http://ai.stanford.edu/events/ai-salon/&sa=D&ust=1550739471132000)’ events on Best Practices in doing Ethical AI Research' %}
[<img class="postimage" src="{{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image7.png"/>]({{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image7.png)
{% endfigure %}


One straightforward principle is to explicitly communicate the
ethical implications of our research whenever we get the opportunity. We
can easily start in our classrooms, by dedicating parts of the syllabus
to address ethical concerns in the field and bring up historical or
current examples of misuse. For instance, we can bring up
the [possibility of unintended bias](https://www.theverge.com/2018/10/10/17958784/ai-recruiting-tool-bias-amazon-report&sa=D&ust=1550739471133000) and
how to guard against it when teaching Machine Learning.  When assigning
large projects, we can provide guidelines for how students can identify
and express concerns about the implications of their work. Further yet,
we can advocate for courses that delve deeper
into these topics, such as Stanford’s [CS122: Artificial Intelligence - Philosophy, Ethics, and Impact](http://web.stanford.edu/class/cs122/&sa=D&ust=1550739471133000), [CS181: Computers, Ethics, and Public Policy](http://web.stanford.edu/class/cs181/&sa=D&ust=1550739471133000),
 and [CS 521: Seminar on AI Safety](https://dorsa.fyi/cs521/&sa=D&ust=1550739471133000). A similar approach could be taken with talks and
interviews: simply allocate a portion of them to explicitly addressing
any ethical concerns in the research. 

**Practice: Take a Stand**

As you develop your own code of ethics,
you might start noticing when other researchers and institutions make
unethical decisions that make you feel uncomfortable. Sometimes, those
institutions might be the company or university you work for, or your
own government. If your institution goes against your code of ethics, we
want to remind you that as an AI researcher, you are not powerless.

{% figure caption:'A snapshot of the open letter by Google employees about Maven' %}
[<img class="postimage_75" src="{{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image1.png"/>]({{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image1.png)
{% endfigure %}

People can influence the sorts of research their institution
pursues through collective action, protests, and even activism.
Recently, more than 4000 employees at Google
signed an open letter against Project Maven, the company’s
contract to develop AI technology for drone video analysis for the
Pentagon, over fears that such technology will be used in drone strikes.
Google [announced](https://www.bloomberg.com/news/articles/2018-10-08/google-drops-out-of-pentagon-s-10-billion-cloud-competition&sa=D&ust=1550739471135000) soon
after that they would not continue with the project, and that they would
not participate in JEDI, the $10 billion cloud contract with the
Department of Defense, citing their AI principles. Similarly, employees
have also protested against [Microsoft](https://medium.com/s/story/an-open-letter-to-microsoft-dont-bid-on-the-us-military-s-project-jedi-7279338b7132&sa=D&ust=1550739471135000)’s
bid in JEDI, and [employees](https://medium.com/s/story/im-an-amazon-employee-my-company-shouldn-t-sell-facial-recognition-tech-to-police-36b5fde934ac&sa=D&ust=1550739471135000)of
[Amazon](https://www.washingtonpost.com/news/the-switch/wp/2018/06/22/amazon-employees-demand-company-cut-ties-with-ice/&sa=D&ust=1550739471135000) have
protested against the company’s work with the US Immigration and Customs
Enforcement (ICE).


**Practice: Obtain and promote more diverse research perspectives**

{% figure caption:'Joy Buolamwini’s TED Talk discussing [her research](http://gendershades.org/index.html&sa=D&ust=1550739471136000) on biases in facial recognition algorithms' %}
[<img class="postimage_75" src="{{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image18.png"/>]({{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/image18.png)
{% endfigure %}

In 2017, Joy Buolamwini discovered that
state-of-the-art facial recognition algorithms performed very poorly for
people of color after testing the algorithms
on herself. The fact that some of our best
algorithms cannot perform well for those that are underrepresented in
the field of AI should not be surprising: as researchers, our approaches
to and perspectives for research can often be limited by our
experiences, histories, and identities. This is why increasing diversity
and making AI a more inclusive place for underrepresented talents can
help our field become more ethical and less biased. 

To be clear, the general benefit of diversity in
research and academia is widely studied and accepted, and is well
explained in this [column](https://www.psychologicalscience.org/observer/diversity-makes-better-science&sa=D&ust=1550739471137000) by the Association for Psychological Science. Our specific
point is that having a more intellectually and experientially diverse
research community will make us much better at identifying the
implications of our research on the world at large, which in turn will
make it easier for the policymakers to come up with more equitable use
cases of AI algorithms. Therefore, we should strive to encourage and
nurture diversity and inclusiveness in our institutions and teams, not
just for better hiring and enrollment numbers, but for richer and more
profound perspectives on our work. 

If you are not in the position to hire and diversify your team, we
have two suggestions. First, you can expand your own intellectual
circle. It is easy to surround yourself with your colleagues and other
AI researchers, and only have research discussions with people in your
field. Try reaching out to thinkers and researchers in other fields,
especially ones in fields that think deeply about ethics and societal
implications of technology, such as philosophy, law, or sociology.
Second, consider mentoring underrepresented researchers. Through
mentorship, you can encourage more diverse talents to join the field and
to give more resources to the people historically underrepresented in
AI. We highly recommend getting involved with programs such as
AI4ALL, Women in AI, and Black in AI. 

**Practice: Large Scale Initiatives**

While all the above steps by us as individuals are collectively
powerful, more directed efforts aimed at dealing with ethical issues in
AI research are also useful. Therefore, we conclude by highlighting some
of the emerging institutions committed to this cause: 

-   [AI](https://ainowinstitute.org/&sa=D&ust=1550739471138000) [Now](https://ainowinstitute.org/&sa=D&ust=1550739471138000) Institute: An NYU research institute that focuses on four domains: rights and liberties, labor and automation, bias and
    inclusion, and safety and critical domains. 
-   [Human-Centered AI Institute](https://hai.stanford.edu/&sa=D&ust=1550739471138000) (HAI):
    A Stanford institute that works on advancing AI research that
    benefits humanity. HAI funds human-centered AI research, and has
    several [fellowships](https://hai.stanford.edu/career/&sa=D&ust=1550739471138000) available. 
-   [Ada Lovelace Institute](https://www.adalovelaceinstitute.org/&sa=D&ust=1550739471139000): An independent research group set up by the Nuffield
    Foundation that examines ethical and social issues arising from the
    use of data, algorithms, and AI. 
-   [The Ethics and Governance of Artificial Intelligence Initiative](https://aiethicsinitiative.org/&sa=D&ust=1550739471139000): A joint project of the MIT Media Lab and the Harvard
    Berkman-Klein Center for Internet and Society, that does both
    philanthropic work as well as research in justice, information
    quality, and autonomy & interaction.  
-   [The IEEE Global Initiative on Ethics of Autonomous and Intelligent Systems](https://standards.ieee.org/industry-connections/ec/autonomous-systems.html&sa=D&ust=1550739471140000):
    An IEEE committee working on the overarching principles of the
    ethical design and use of autonomous and intelligent systems
-   [Partnership on AI](https://www.partnershiponai.org/&sa=D&ust=1550739471140000): A multilateral organization that brings together
    companies, academics, researchers, and civil society organizations
    to formulate best practices for AI technologies. 
-   [The Future of Life Institute](https://futureoflife.org/):  An organization with the mission to catalyze and support research and initiatives for safeguarding life and developing optimistic visions of the future, which organizes [many events and discussions about the future of AI](https://futureoflife.org/ai-activities/). 
    
{% figure caption:'Group photo from [the Future of Life Beneficial AGI 2019 event](https://futureoflife.org/beneficial-agi-2019/)' %}
[<img class="postimage" src="{{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/asimolar.jpg"/>]({{ site.baseurl }}/assets/img/posts/2019-02-21-ethical_best_practices/asimolar.jpg)
{% endfigure %}

There are of course many other institutions and labs that do work
focusing on the ethics and policy of AI, and this list is
non-comprehensive.

### Conclusion
The question of the relationship between creations and
creators is not a new one. Regardless of one’s opinions on this
question, the undeniable fact is that AI has more potential to change
the landscape of our civilization than perhaps any other human
invention. The nature of this change will depend greatly on our
collective understanding of the benefits and limitations of various AI
technologies, and this understanding depends greatly on the engagement
of researchers with policymakers, legislators and the broader public as
well. The time for taking shelter behind the ivory tower of academic
immunity and washing our hands of the implications of our work is over.
Let us all come together, and in whatever way we find suitable,
contribute to a more profound and honest understanding of our field, for
the betterment of human civilization.


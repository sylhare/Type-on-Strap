---
layout: post
title: "Deep Learning, Structure and Innate Priors"
excerpt: "A Discussion between Yann LeCun and Christopher Manning about the future of Deep Learning and NLP"
summary: "A Discussion between Yann LeCun and Christopher Manning about the future of Deep Learning and NLP"
feature-img: "assets/img/posts/2018-12-03-deep-learning-structure-and-innate-priors/feature.png"
thumbnail: "assets/img/posts/2018-12-03-deep-learning-structure-and-innate-priors/thumb.png"
author: <a href='https://cs.stanford.edu/people/abisee/'>Abigail See</a>
tags: [AI Salon, video]
---

*This blog post was originally posted on [Abigail See's blog](http://www.abigailsee.com/2018/02/21/deep-learning-structure-and-innate-priors.html) on February 21 of 2018. The event it covers happened on February 2nd of 2018, as part of [SAIL's regular AI Salon](http://ai.stanford.edu/events/ai-salon/) discussion series. Look forward to more blog posts on new AI Salon events soon!*

---

<!--excerpt.start-->
{% include yt.html %}
<br><br>
Earlier this month, I had the exciting opportunity to moderate a discussion between Professors [Yann LeCun](http://yann.lecun.com/) and [Christopher Manning](https://nlp.stanford.edu/manning/), titled **_"What innate priors should we build into the architecture of deep learning systems?"_** The event was a special installment of [AI Salon](http://ai.stanford.edu/events/ai-salon/), a discussion series held within the Stanford AI Lab that often features expert guests.
<!--excerpt.end-->

This discussion topic – about the structural design decisions we build into our neural architectures, and how those correspond to certain assumptions and inductive biases – is an important one in AI right now. In fact, [last year I highlighted](http://www.abigailsee.com/2017/08/30/four-deep-learning-trends-from-acl-2017-part-1.html) "the return of linguistic structure" as one of the top four NLP Deep Learning research trends of 2017.

On one side, Manning is a prominent advocate for incorporating _more_ linguistic structure into deep learning systems. On the other, LeCun is a leading proponent for the ability of simple but powerful neural architectures to perform sophisticated tasks _without_ extensive task-specific feature engineering. For this reason, anticipation for disagreement between the two was high, with [one Twitter commentator](https://twitter.com/saiabishek1/status/959025926737670144) describing the event as "the AI equivalent of Batman vs Superman".

However, LeCun and Manning agreed on more than you may expect.
LeCun's most famous contribution (the Convolutional Neural Network) is _all about_ an innate prior -- the assumption that an image processing system should be [translationally invariant](https://www.quora.com/Why-and-how-are-convolutional-neural-networks-translation-invariant) -- which is enforced through an architectural design choice (weight sharing).
For his part, Manning has spoken publicly to say that the Deep Learning renaissance is [A Good Thing for NLP](https://www.mitpressjournals.org/doi/pdf/10.1162/COLI_a_00239).

While the two professors agreed on many other things during the discussion, certain key differences emerged -- you can watch the full video above. **The rest of this post is a summary of the main themes that emerged throughout the discussion**, plus some links to relevant further materials.

### Structure: a necessary good or a necessary evil?

In their opening statements, Manning and LeCun quickly established their main difference of opinion.

Manning described structure as a "necessary good" ([9:14](javascript:goTo(9,14))), arguing that we should have a positive attitude towards structure as a good design decision. In particular, structure allows us to design systems that can learn more from less data, and at a higher level of abstraction, compared to those without structure.

Conversely, LeCun described structure as a "necessary evil" ([2:44](javascript:goTo(2,44))), and warned that imposing structure requires us to make certain assumptions, which are invariably wrong for at least some portion of the data, and may become obsolete within the near future. As an example, he hypothesized that ConvNets may be obsolete in 10 years ([29:57](javascript:goTo(29,57))).

Despite this disagreement, we should note that LeCun and Manning did at least agree that structure is "necessary" – they just have different attitudes towards that necessity.

Manning views it as the right and principled thing to do – for example, language is fundamentally recursive, so NLP architectures should be too ([23:40](javascript:goTo(23,40)))!
He did acknowledge, however, that in practice it's difficult to make the correct structural assumptions, and those assumptions don't always translate to comprehensive performance gains (see for example, the mixed success of the [Recursive Neural Network](https://en.wikipedia.org/wiki/Recursive_neural_network), aka Tree-RNN, which imposes recursive compositionality as an innate prior).

LeCun has a much less idealized view of structure. Several times during the discussion, he referred to various types of structure (e.g. residual connections, convolutions), as merely "a meta-level substrate" ([53:33](javascript:goTo(53,33))) that is required for optimization to work. A similar network without the structural constraints, he claimed, would work just as well, except it would take longer to train.


### The limitations of today's AI

LeCun and Manning noted the historical trajectory that has brought us to this present moment in AI research. Over the last few decades, innate priors have gone out of fashion, and today Deep Learning research prizes closely-supervised end-to-end learning (supported by big-data and big-compute) as the dominant paradigm.

Both LeCun and Manning repeatedly highlighted the limitations of this paradigm – for example the progress that remains to be made on memory, planning, transfer learning, world knowledge, and multi-step reasoning – and expressed positivity ([22:17](javascript:goTo(22,17)), [37:20](javascript:goTo(37,20)), [57:28](javascript:goTo(57,28))) towards current research that aims to tackle these problems via structural design decisions.
<!-- For example, the past few years has seen a rapidly expanding body of work (including [Memory Networks](https://arxiv.org/abs/1503.08895), [Neural Turing Machines](https://arxiv.org/abs/1410.5401), [Differentiable Neural Computers](https://deepmind.com/blog/differentiable-neural-computers/), and others) equipping neural architectures with memory in order to effectively perform multi-step reasoning tasks. -->

However, Manning went further, asserting that the big-data big-compute paradigm of modern Deep Learning has in fact "perverted the field" (of computational linguistics) and "sent it off-track" ([10:48](javascript:goTo(10,48))). If you have access to huge amounts of data and computation, he argued, you can succeed by building simple but inefficient systems that perform "glorified nearest neighbor learning" at a superficial level ([43:20](javascript:goTo(43,20))). This disincentivizes researchers from building good learning systems – ones which learn representations at a higher level of abstraction, and do not require huge amounts of data. This, he said, is bad for the field as a whole. The answer? Impose the _right_ kind of innate structure, that enables systems to learn concepts efficiently at the right level of abstraction.

Despite my attempt to prod the two into conflict ([33:15](javascript:goTo(33,15))), I'm unsure what exactly LeCun thought of Manning's claim that Deep Learning has in some sense "perverted the field". However, LeCun did agree ([34:30](javascript:goTo(34,30))) that Deep Learning _is_ missing basic principles (to read more on that topic, see his CVPR'15 keynote, _[What's Wrong With Deep Learning?](http://www.pamitc.org/cvpr15/files/lecun-20150610-cvpr-keynote.pdf))_.


### The importance of unsupervised learning

While the discussion touched upon many core limitations of today's AI techniques, one particular challenge – which may be loosely described as Unsupervised Learning, or at least Less-Supervised Learning – emerged as a matter of particular urgency.

Both professors gave examples ([9:48](javascript:goTo(9,48)), [30:30](javascript:goTo(30,30))) of humans' ability to do few-shot learning; to learn about the world via observation, without a task or an external reward; and to learn abstract concepts with discrete structure (for example, categorization of objects) without explicit supervision.

These unsupervised learning abilities, they agreed, are essential to progress in AI.
But when it came to the role _structure_ should play in the [Unsupervised Revolution](https://twitter.com/rgblong/status/916062474545319938), however, LeCun and Manning disagreed.

Manning argued that imposing structure is the key to unlock unsupervised learning ([35:05](javascript:goTo(35,05))). If we provide machines with the right structural tools to learn at an appropriate level of abstraction, he said, then they can learn with less supervision.

By contrast, LeCun argued that if you can perform unsupervised learning, you don't _need_ to impose structure. As an example ([28:57](javascript:goTo(28,57))), he described how the human brain does not have any innate convolutional structure – but it doesn't need to, because as an effective unsupervised learner, the brain can learn the same low-level image features (e.g. oriented edge detectors) as a ConvNet, even without the convolutional weight-sharing constraint. He concluded that imposing _more_ structure on our current neural architectures may be futile, because once we have developed better methods for unsupervised learning, those structural design decisions may be obsolete.

The difference between the two positions was subtle; and perhaps mostly a chicken-and-egg distinction. Manning regards structure as an important key to achieve unsupervised learning, whereas LeCun regards unsupervised learning as the only long-term way to learn structure.


### Structure as a hard-wired prior, or learned from the environment?

During the discussion, it became clear that there are at least two types of "structure": structure baked into the model as an innate prior (for example, the convolutional assumption in ConvNets, or the recursive assumption in Recursive Neural Networks), and structure learned and computed dynamically by the machine (for example, the structure computed by dynamic routing in [Capsule Networks](https://arxiv.org/abs/1710.09829), or the alignments computed by the [attention mechanism](https://arxiv.org/pdf/1409.0473.pdf)).
There is no easy distinction between the two, and at one point Manning and LeCun differed on whether ConvNets' hierarchical structure should be regarded as one or the other ([25:55](javascript:goTo(25,55))).

LeCun repeatedly spoke against what he called hard-wired priors, arguing that all structure should instead be learned from the environment ([30:42](javascript:goTo(30,42)), [34:14](javascript:goTo(34,14))). Though Manning agreed that much structure should be learned from the environment, he also argued that we (the designers of AI systems) should play _some_ part in providing that structure. While we shouldn't return to the days of intricately human-designed systems (such as Chomskyan grammars), he said, we should provide machines with the right "primitives and scaffolding" to learn more effectively ([11:37](javascript:goTo(11,37))).


### Reward as an innate prior
LeCun and Manning agreed that ideally, reward should be _innate_ – that is, understanding the world correctly should be its own reward ([46:03](javascript:goTo(46,03))). For example, humans are constantly building their own internal model of the world, and revising it in response to external observations.

By contrast, most Machine Learning systems today learn from externally-provided rewards that are closely related to a particular task. Manning described these objective functions as too superficial – noting that we will never build AI systems that learn abstract concepts if the objective function is defined at such a low level ([37:55](javascript:goTo(37,55))). LeCun agreed that reward needs to be intrinsic, and rich – rather than learning from occasional task-specific rewards, AI systems should learn by constantly predicting "everything from everything", without requiring training labels or a task definition ([49:14](javascript:goTo(49,14))).


### On language
In the final minutes of the discussion, LeCun, perhaps being a little provocative, claimed language is "not that complicated", nor that crucial to achieving general intelligence ([59:54](javascript:goTo(59,54))). To support this, he appealed to the fact that orangutans are almost as intelligent as humans, yet they have no language. In response, Manning leaped to the defense of language – which, he claimed, is crucial to general intelligence, because language is the conduit by which individual intelligence is shared and transformed into societal intelligence!

### Miscellaneous notes and further reading

For convenience, here is a (non-comprehensive) list of some papers, ideas and resources mentioned or otherwise relevant to the discussion. There were some references mentioned in the discussion that I was unable to find, so please contribute any further links in the comments!

* At [19:17](javascript:goTo(19,17)), Manning discusses the paper _[Grammar as a Foreign Language](https://arxiv.org/abs/1412.7449)_, which tackled a highly _recursive_ linguistic task (parsing) with a surprisingly _unstructured_ method (sequence-to-sequence).
* The question at [39:15](javascript:goTo(39,15)) refers to the idea that Stochastic Gradient Descent acts as a kind of implicit regularization. To read more about this idea, see for example the work of Tomaso Poggio and his collaborators. [Here](https://stats385.github.io/assets/lectures/StanfordStats385-20171025-Lecture05-Poggio.pdf) is a set of slides he presented at Stanford's [Theory of Deep Learning class](https://stats385.github.io/) last year – slide 44 shows the connection between SGD and implicit regularization. More generally, Poggio and his collaborators are one of the many theorists LeCun mentions as investigating "the theoretical mystery" of why neural nets work ([41:29](javascript:goTo(41,29))).
* At [40:52](javascript:goTo(40,52)), I question whether bigger networks are necessarily better, and mention a paper that shows this is not always true. I was referring to the [ResNet paper](https://arxiv.org/pdf/1512.03385.pdf), which demonstrates that deeper networks can be harder to train than shallower networks, and thus sometimes achieve worse results. However, the same paper then shows that residual connections (the paper's main contribution) provide a way to train deep nets much more effectively. So perhaps "bigger isn't always better" isn't a fair conclusion to draw from the paper – "bigger is only better if you can train it effectively" would be more precise! For more thoughts on whether bigger is better, see _[Do Deep Nets Really Need to be Deep?](https://papers.nips.cc/paper/5484-do-deep-nets-really-need-to-be-deep.pdf)_
* The "where do rewards come from?" question at [45:52](javascript:goTo(45,52)) mentions [Andrew Barto](http://www-all.cs.umass.edu/~barto/), who has written a paper [with that exact title](http://www-all.cs.umass.edu/pubs/2009/singh_l_b_09.pdf).
* At [57:46](javascript:goTo(57,46)) LeCun mentions a paper by [Leon Bottou](http://leon.bottou.org/) on the idea of mapping representations back to the same space, thus enabling chains of reasoning. The paper is called _[From Machine Learning to Machine Reasoning](https://arxiv.org/abs/1102.1808)_.
* In October 2017, Yann LeCun took part in a debate with Gary Marcus at NYU, with a similar discussion topic to ours – _"Does AI Need More Innate Machinery?"_. It is a highly interesting discussion, and I recommend you watch it [here](https://www.youtube.com/watch?v=vdWPQ6iAkT4&feature=youtu.be). The two have since had [further](https://twitter.com/ylecun/status/952587501037916161) [disagreement](https://twitter.com/ylecun/status/953033413807755264) on Twitter on the subject of Deep Learning.

---

*Thanks to both Yann LeCun and Christopher Manning for sharing their perspectives with us in this discussion. Special thanks to [Siva Reddy](http://sivareddy.in/) for organizing much of the event.*

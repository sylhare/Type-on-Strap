---
layout: post
title: Powerful Abstractions for Programmatically Building and Managing Training Sets
short-summary: "Machine learning practitioners are spending less time on model architectures and hardware optimizations and, instead, focusing on training data. We describe three powerful abstractions that practitioners can use to programmatically build and manage their training data."
summary: "Machine learning practitioners are spending less time on model architectures and hardware optimizations and, instead, focusing on training data. We describe three powerful abstractions that practitioners can use to programmatically build and manage their training data."
thumbnail: "assets/img/posts/2019-06-21-training-data-abstractions/fig_abstractions_thumbnail.png"
author: <a href='https://vincentsc.com'>Vincent S. Chen</a>, <a href='https://stanford.edu/~senwu/'>Sen Wu</a>, <a href='https://www.bradenhancock.com'>Braden Hancock</a>, <a href='https://ajratner.github.io'>Alex Ratner</a>, <a href='https://cs.stanford.edu/people/chrismre/'>Chris Ré</a>, and&nbsp;<a href='https://cs.stanford.edu/people/chrismre/#students'>other&nbsp;members&nbsp;of&nbsp;Hazy&nbsp;Lab</a>
tags: [ml,systems,research]
---

## **Overview**

Machine learning practitioners are spending less time on model architectures and hardware optimizations and, instead, focusing on training data. As a result, programmers are relying on different abstractions—high-level design patterns—to build machine learning pipelines for their applications. In this post, we describe three powerful abstractions that practitioners can use to programmatically build and manage their training data.

We ran an experiment to test the effectiveness of basic training data operations—applying a handful of these using our framework, [Snorkel](http://snorkel.stanford.edu), and a standard NLP model (i.e. BERT) yields a state-of-the-art result on [SuperGLUE](https://super.gluebenchmark.com/)[^superglue]—a newly curated benchmark with six tasks for evaluating “general-purpose language understanding technologies”. Compared to the recent advances in natural language pretraining (i.e. BERT), we achieve a _<mark>new state-of-the-art score overall and the highest reported score anywhere on a majority of component tasks</mark>_.

Beyond SuperGLUE, we also highlight updates on Snorkel’s use in the real world with even more applications—from industrial scale at [Google’s Snorkel Drybell](https://ai.googleblog.com/2019/03/harnessing-organizational-knowledge-for.html) to scientific work in [MRI classification](https://nature-research-under-consideration.nature.com/users/37265-nature-communications/posts/38921-weakly-supervised-classification-of-rare-aortic-valve-malformations-using-unlabeled-cardiac-mri-sequences) and [automated Genome-wide association study (GWAS) curation](https://ai.stanford.edu/~kuleshov/papers/gwaskb-manuscript.pdf) (both accepted in [Nature Comms](https://www.nature.com/ncomms/))!

We will be releasing code in the [Snorkel repo](https://github.com/HazyResearch/snorkel) for reproducing and building on our results in conjunction with a 2-day **Snorkel workshop** during the last week  of June with collaborators from science, industry, and government. This workshop is unfortunately already completely full, but if you would like to be notified of future Snorkel workshops, please provide your name and contact information [here](https://docs.google.com/forms/d/e/1FAIpQLScOpiImyBA3uk_CnJ03R1b7Ese9VA3XjfLnemCO76WyTwrO5Q/viewform?usp=sf_link).


## **Three key abstractions**

In our SuperGLUE result, as well as more generally, we find that spending our time programmatically building and manipulating the training data—rather than the models— provides a powerful and effective strategy to achieve high performance in ML pipelines. In a past [post](https://dawn.cs.stanford.edu/2019/03/22/glue/), we talked about the value of incorporating more supervision signal from more sources, e.g. multi-task learning and transfer learning, as we achieved state-of-the-art results on the GLUE Benchmark (a precursor to SuperGLUE). In this post, we focus on three key abstractions for building and modifying training datasets:


1. **Labeling data** with labeling functions (LFs) [^dp]
2. **Transforming data** with transformation functions (TFs) [^tanda] [^autoaugment]
3. **Slicing data** with slicing functions (SFs) [_technical report + blog post coming soon!_]

{% figure %}
[<img class="postimage_75" src="{{ site.baseurl }}/assets/img/posts/2019-06-21-training-data-abstractions/fig_abstractions.png"/>]({{ site.baseurl }}/assets/img/posts/2019-06-21-training-data-abstractions/fig_abstractions.png)
{% endfigure %}


### **Running Example**

For the remainder of this post, we use a running example from the Words in Context (WiC) task from SuperGLUE: _is the target word being used in the same way in both sentences?_

{% figure %}
[<img src="{{ site.baseurl }}/assets/img/posts/2019-06-21-training-data-abstractions/example.png"/>]({{ site.baseurl }}/assets/img/posts/2019-06-21-training-data-abstractions/example.png)
{% endfigure %}


## **1. Weak labeling with labeling functions**

In many applications, unlabeled data is abundant—it may come from fleets of autonomous vehicles, or large corpora of unstructured data. Modern architectures are largely unable to take advantage of such potentially rich datasets because labeling them is intractable due to time or cost. With [Snorkel](https://hazyresearch.github.io/snorkel/), we’ve studied for years the use of **labeling functions (LFs)** for heuristically labeling training examples. LFs provide domain experts or machine learning practitioners with an intuitive interface for denoising and combining supervision sources from existing datasets, models, or crowd labelers.

{% figure %}
[<img src="{{ site.baseurl }}/assets/img/posts/2019-06-21-training-data-abstractions/lf_ex.png"/>]({{ site.baseurl }}/assets/img/posts/2019-06-21-training-data-abstractions/lf_ex.png)
<figcaption>
  For the WiC task (identifying whether a target word is used with the same "sense" in two sentences) we might consider weakly labeling examples based on whether or not they share a trigram including the target word.
</figcaption>
{% endfigure %}


## **2. Augmenting data with transformation functions**

Often, people think about data augmentation in terms of simple transformations—randomly rotating or stretching images—but they can refer to much more diverse range of operations. We see **transformation functions (TFs)** as a powerful abstraction that heuristically generates new, modified examples from existing ones. For instance, for a medical imaging task, we might write TFs to perform transformations that are specific to our imaging modality—e.g. resampling segmenting tumor masses or resampling background tissue. We have explored this abstraction in our own work, TANDA [^tanda-2], which seeks to learn compositions of transformations across domain-specific tasks. AutoAugment [^autoaugment-2] from Google builds on this work to automatically learn policies for augmentation strategies.

{% figure %}
[<img src="{{ site.baseurl }}/assets/img/posts/2019-06-21-training-data-abstractions/tf_ex.png"/>]({{ site.baseurl }}/assets/img/posts/2019-06-21-training-data-abstractions/tf_ex.png)
<figcaption>
  Given that “Sunday” does not change the word sense of “invite”, we can transform an example that includes the word “Sunday” into many copies of that example with different days of the week so that our model is likely to overfit.
</figcaption>
{% endfigure %}


## **3. Slicing data with slicing functions (<mark>new idea</mark>!)**

In many datasets, especially in real-world applications, there are subsets of the data that our model underperforms on, or that we care more about performing well on than others. For example, a model may underperform on lower-frequency healthcare demographics (e.g. younger patients with certain cancers) or we may care extra about model performance on safety-critical but rare scenarios in an autonomous driving setting, such as detecting cyclists. We call these data subsets _slices_. The technical challenge often faced by practitioners is to improve performance on these slices while maintaining overall performance.

**Slicing functions (SFs)** provide an interface for users to coarsely identify data subsets for which the model should commit additional representational capacity. To address slice-specific representations, practitioners might train many models that each specialize on particular subsets, and then combine these with a mixture-of-experts (MoE) approach [^moe]. However, with the growing size of ML models, MoE is often impractical. Another strategy would be to train a single model in the style of multi-task learning (MTL) with hard parameter sharing [^mtl]. While more computationally efficient, this approach expects representation bias across many slice-specific tasks to improve performance—an often unreliable approach. As a quick overview (_technical report + blog post coming soon!_)— we model slices in the style of multi-task learning, in which slice-based “expert-heads” are used to learn slice-specific representations. Then, an attention mechanism is learned over expert heads to determine when and how to combine the representations learned by these slice heads on a per-example basis.

We consider the following properties of our approach:



*   Our approach is **model-agnostic** — expert heads are learned on top of any backbone architecture (e.g. BERT, ResNET). As a result, practitioners improving performance with slicing functions can focus on the data rather than the model architecture.
*   By learning in a multi-task fashion, we **efficiently learn representations** without the need to make many copies of the model (i.e. MoE requires too much memory)!
*   By incorporating the attention mechanism, we **avoid manual tuning** of expert-heads—an otherwise significant developer cost.


{% figure %}
[<img src="{{ site.baseurl }}/assets/img/posts/2019-06-21-training-data-abstractions/sf_ex.png"/>]({{ site.baseurl }}/assets/img/posts/2019-06-21-training-data-abstractions/sf_ex.png)
<figcaption>
  From WiC error analysis, we might find that our model appears to perform worse on examples where the target word is a noun instead of a verb. Using an SF, we tell the model to pay attention to the differences between these slices and use a slightly different representation when making predictions for target words that it believes are nouns.
</figcaption>
{% endfigure %}


## **Key properties of LFs, TFs, and SFs**

*   **Intuitive interfaces**: These abstractions provide intuitive interfaces to existing practitioner workflows. They allow insights from debugging/error analysis to be directly encoded to improve models.
*   **Programming abstractions as weak supervision**: In practice, many of these techniques can be viewed as a form of weak supervision, as users specify them in noisy, heuristic, and imprecise ways. Dealing with this is one of the core technical challenges we tackle with Snorkel.
*   **Supervision as code**: These types of inputs are ways of supervising a model (i.e. they specify training sets). Concretely, they are also code, and thus carry many of the advantages of code—reusability, modifiability, etc.


## **SuperGLUE Results**

Using these programming abstractions, we achieved a new state-of-the-art score on the SuperGLUE Benchmark and 4 of its components tasks. SuperGLUE is similar to [GLUE](https://gluebenchmark.com/), but contains “more difficult tasks...chosen to maximize difficulty and diversity, and...selected to show a substantial headroom gap between a strong BERT-based baseline and human performance.” After reproducing the BERT++ baselines, we minimally tuned these models (baseline models, default learning rate, etc.) and found that with a handful of applications of the above programming abstractions, we saw improvements of +4.0 points on the SuperGLUE benchmark (21% reduction of the gap to human performance).


## **Snorkel in the Real World**

These Snorkel programming abstractions have also been used to fuel progress in high-impact real-world applications.

In March of this year, we published a [paper](https://arxiv.org/pdf/1812.00417.pdf) and [blog post](https://ai.googleblog.com/2019/03/harnessing-organizational-knowledge-for.html) with Google on the lessons learned from deploying Snorkel at industrial scale. Relying on diverse sources of knowledge across the organization—heuristics, taggers, knowledge graphs, legacy systems, etc.—they saw significant improvements in quality, by as much as 17.5 F1 points.

{% figure %}
[<img src="{{ site.baseurl }}/assets/img/posts/2019-06-21-training-data-abstractions/bav.jpg"/>]({{ site.baseurl }}/assets/img/posts/2019-06-21-training-data-abstractions/bav.jpg)
<figcaption>
  The Snorkel pipeline, deployed on the BAV classification task for large collections of up to 4,000 unlabeled MRI sequences. Figure credit to Fries et. al 2018.
</figcaption>
{% endfigure %}

In [recent work](https://www.biorxiv.org/content/10.1101/339630v4.full) that was accepted to Nature Communications, Snorkel was deployed in an ongoing collaboration with [Stanford University Pediatric Cardiology](https://priestlab.stanford.edu/), where labeled training data is a significant practical roadblock to developing automated methods. We focused on bicuspid aortic valve (BAV), the most common congenital heart malformation (with an incidence rate of 0.5-2% in the general population), with risk of adverse downstream health effects. Instead of relying on costly MRI labels from cardiologists, we worked directly with domain experts to develop LFs to generate large-scale training sets for downstream deep learning models. In patients identified by our end-to-end approach, an independent evaluation determined a 1.8-fold increase in risk for major adverse cardiac events.

In another forthcoming Nature Communications [paper](https://ai.stanford.edu/~kuleshov/papers/gwaskb-manuscript.pdf), we showed how Snorkel can be used to automate Gene-Wide Association Study (GWAS) curation. On a collection of hundreds of previously published studies reporting significant genotype-phenotype pairs, we auto-labeled a large training set using only labeling functions. The resulting classifier applied to a collection of 598 studies recovered over 3,000 previously documented open-access relations (with an estimated recall of 60-80%) as well as over 2,000 associations not present in existing human curated repositories (with an estimated precision of 82-89%). The resulting database is available for exploration with a user interface at [http://gwaskb.stanford.edu/](http://gwaskb.stanford.edu/).


## **Stay Tuned**

The Snorkel project is active and ongoing! We have a number of exciting, ongoing collaborations—from follow-on work at Stanford’s School of Medicine, to deployments at the [International Consortium of Investigative Journalists (ICIJ)](https://www.icij.org/blog/2019/03/how-artificial-intelligence-can-help-us-crack-more-panama-papers-stories/) to help journalists organize, index, and understand millions of unstructured documents.

A code release later this month will include significant infrastructural improvements and tutorials for how to apply LFs, TFs, and SFs to SuperGLUE and other tasks. If you've used Snorkel for your own applications, we'd love to hear about it! For updates on Snorkel developments and applications, you can always visit the Snorkel [landing page](http://snorkel.stanford.edu/) or [open-source repository](https://github.com/HazyResearch/snorkel).

## Acknowledgements
The authors would like to thank Feng Niu and Charles Srisuwananukorn for many helpful discussions, tests, and collaborations throughout the development of slicing!

<!-- ##### Footnotes -->
<!-- * footnotes will be placed here. This line is necessary -->
<!-- {:footnotes} -->
[^superglue]: Wang, Alex, et al. ["SuperGLUE: A Stickier Benchmark for General-Purpose Language Understanding Systems."](https://arxiv.org/abs/1905.00537). 2019.  _SuperGLUE_ consists of 6 datasets: the Commitment Bank (CB, [De Marneffe et al., 2019](https://github.com/mcdm/CommitmentBank/), Choice Of Plausible Alternatives (COPA, [Roemmele et al., 2011](https://www.aaai.org/ocs/index.php/SSS/SSS11/paper/viewPaper/2418)), the Multi-Sentence Reading Comprehension dataset (MultiRC, [Khashabi et al., 2018](https://www.aclweb.org/anthology/papers/N/N18/N18-1023/)), Recognizing Textual Entailment (merged from RTE1, [Dagan et al. 2006](https://link.springer.com/chapter/10.1007/11736790_9), RTE2, [Bar Haim et al., 2006](http://u.cs.biu.ac.il/~nlp/downloads/publications/RTE2-organizers.pdf), RTE3, [Giampiccolo et al., 2007](https://dl.acm.org/citation.cfm?id=1654538), and RTE5, [Bentivogli et al., 2009](http://www.cs.utexas.edu/users/pclark/papers/RTE6_overview.proceedings.pdf)), Word in Context (WiC, [Pilehvar and Camacho-Collados, 2019](https://www.aclweb.org/anthology/papers/N/N19/N19-1128)), and the Winograd Schema Challenge (WSC, [Levesque et al., 2012](https://www.aaai.org/ocs/index.php/KR/KR12/paper/viewPaper/4492)).

[^dp]: Ratner, Alexander J., et al. ["Data programming: Creating large training sets, quickly."](http://papers.nips.cc/paper/6523-data-programming-creating-large-training-sets-quickly) 2016.

[^tanda]: Ratner, Alexander J., et al. ["Learning to compose domain-specific transformations for data augmentation."](http://papers.nips.cc/paper/6916-learning-to-compose-domain-specific-transformations-for-data-augmentation) 2017.

[^tanda-2]: Ratner, Alexander J., et al. ["Learning to compose domain-specific transformations for data augmentation."](http://papers.nips.cc/paper/6916-learning-to-compose-domain-specific-transformations-for-data-augmentation) 2017.

[^autoaugment]: Cubuk, Ekin D., et al. ["Autoaugment: Learning augmentation policies from data."](https://arxiv.org/abs/1805.09501) 2018.

[^autoaugment-2]: Cubuk, Ekin D., et al. ["Autoaugment: Learning augmentation policies from data."](https://arxiv.org/abs/1805.09501) 2018.


[^moe]: Robert A Jacobs, Michael I Jordan, Steven J Nowlan, and Geoffrey E Hinton. ["Adaptive mixtures of local experts."]([http://www.csri.utoronto.ca/~hinton/absps/jjnh91.ps](http://www.csri.utoronto.ca/~hinton/absps/jjnh91.ps)) 1991.

[^mtl]: Rich Caruana. ["Multitask learning."](https://link.springer.com/article/10.1023/A:1007379606734) 1997.

[^bav]: Fries, Jason A., et al. ["Weakly supervised classification of rare aortic valve malformations using unlabeled cardiac MRI sequences."](https://www.biorxiv.org/content/10.1101/339630v4.full)2018.

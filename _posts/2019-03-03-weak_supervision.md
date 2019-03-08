---
layout: post
title: "Weak Supervision: A New Programming Paradigm for Machine Learning"
short-summary: "Using weak supervision to efficiently label training data."
summary: "Using weak supervision, or high-level noisy sources of labels, to efficiently label training data."
thumbnail: "assets/img/posts/2019-02-26-beyond_local_pattern_matching/thumb.png"
author: <a href='https://www.paroma.xyz'>Paroma Varma</a>, <a href='https://ajratner.github.io'>Alex Ratner</a>, <a href='https://www.bradenhancock.com'>Braden Hancock</a>, <a href='https://cs.stanford.edu/people/chrismre/'>Chris Ré</a>, and <a href='https://cs.stanford.edu/people/chrismre/#students'>other members of Hazy Lab</a>
tags: [nlp,ws]
---

In recent years, the real-world impact of machine learning has grown in leaps and bounds. In large part, this is due to the advent of deep learning models, which allow practitioners to get state-of-the-art scores on benchmark datasets without any hand-engineered features. Given the availability of multiple open-source machine learning frameworks like TensorFlow and PyTorch, and an abundance of available state-of-the-art models, it can be argued that high-quality machine learning models are almost a commoditized resource now. There is a hidden catch, however: **the reliance of these models on massive sets of hand-labeled training data.**

These hand-labeled training sets are expensive and time-consuming to create—often requiring person-months or years to assemble, clean, and debug --- a prohibitively expensive and slow task, especially when domain expertise is required. On top of this, tasks often change and evolve in the real world. For example, labeling guidelines, granularities, or downstream use cases often change, necessitating re-labeling (e.g., instead of classifying reviews only as positive or negative, introducing a neutral category). For all these reasons, practitioners have increasingly been turning to **weaker forms of supervision**, such as heuristically generating training data with external knowledge bases, patterns or rules, or other classifiers. Essentially, these are all ways of programmatically generating training data—or, more succinctly, **programming training data.** 

We begin by reviewing areas of machine learning that are motivated by the problem of labeling training data, and then describe our research that addresses challenges of modeling and integrating a diverse set of supervision sources. We also discuss our vision for building data management systems for the massively multi-task regime with tens or hundreds of weakly supervised *dynamic* tasks interacting in complex and varied ways. Check out the [our research blog](snorkel.stanford.edu) for detailed discussions of these topics and more!

## How to Get More Labeled Training Data? A Review
Many traditional lines of research in machine learning are similarly motivated by the insatiable appetite of modern machine learning models for labeled training data. We start by drawing the core distinction between these other approaches and weak supervision at a high-level: **weak supervision is about leveraging higher-level and/or noisier input from subject matter experts (SMEs).**

{% figure %}
[<img class="postimage" src="{{ site.baseurl }}/assets/img/posts/2019-03-03-weak_supervision/WS_mapping.png"/>]({{ site.baseurl }}//assets/img/posts/2019-03-03-weak_supervision/WS_mapping.png)
{% endfigure %}

The problem is that this is expensive: for example, unlike grad students, radiologists don’t generally accept payment in burritos and free T-shirts! Thus, many well-studied lines of work in machine learning are motivated by the bottleneck of getting labeled training data:

* In **active learning**, the goal is to make use of subject matter experts more efficiently by having them label data points which are estimated to be most valuable to the model (for a good survey, see (Settles 2012)). In the standard supervised learning setting, this means selecting new data points to be labeled. For example, we might select mammograms that lie close to the current model decision boundary, and ask radiologists to label only these. However, we could also just ask for weaker supervision pertinent to these data points, in which case active learning is perfectly complementary with weak supervision; as one example of this, see [(Druck, Settles, and McCallum 2009)](https://aclanthology.info/pdf/D/D09/D09-1009.pdf).  
&nbsp;
* In the **semi-supervised learning** setting, the goal is to use both a small labeled training set and a much larger unlabeled data set. At a high level, we then use assumptions about smoothness, low dimensional structure, or distance metrics to leverage the unlabeled data (either as part of a generative model, as a regularizer for a discriminative model, or to learn a compact data representation); for a good survey see [(Chapelle, Scholkopf, and Zien 2009)](http://www.acad.bg/ebook/ml/MITPress-%20SemiSupervised%20Learning.pdf). Broadly, rather than soliciting more input from subject matter experts, the idea in semi-supervised learning is to leverage domain and task-agnostic assumptions to exploit the unlabeled data that is often cheaply available in large quantities. More recent methods use generative adversarial networks [(Salimans et al. 2016)](http://papers.nips.cc/paper/6125-improved-techniques-for-training-gans.pdf), heuristic transformation models [(Laine and Aila 2016)](https://arxiv.org/pdf/1610.02242), and other generative approaches to effectively help regularize decision boundaries.   
&nbsp;
* In a typical **transfer learning** setting, the goal is to take one or more models already trained on a different dataset and apply them to our dataset and task; for a good overview see [(Pan and Yang 2010)](https://www.cse.ust.hk/~qyang/Docs/2009/tkde_transfer_learning.pdf). For example, we might have a large training set for tumors in another part of the body and classifiers trained on this set and wish to apply these to our mammography task. A common transfer learning approach in the deep learning community today is to “pre-train” a model on one large dataset, and then “fine-tune” it on the task of interest. Another related line of work is multi-task learning, where several tasks are learned jointly ([Caruna 1993](https://www.semanticscholar.org/paper/Multitask-Learning%3A-A-Knowledge-Based-Source-of-Caruana/9464d15f4f8d578f93332db4aa1c9c182fd51735); [Augenstein, Vlachos, and Maynard 2015](http://www.aclweb.org/anthology/D15-1086)). 
<!-- Some transfer learning approaches take one or more pre-trained models (potentially with some heuristic conditioning of when they are each applied) and use these to train a new model for the task of interest; in this case, we can actually consider transfer learning as a type of weak supervision.   -->

The above paradigms potentially allow us to avoid asking our SME collaborators for additional training labels. But what if we could ask them for various types of higher-level, or otherwise less precise, forms of supervision, which would be faster and easier to provide? For example, what if our radiologists could spend an afternoon specifying a set of heuristics or other resources, that, if handled properly, could effectively replace thousands of training labels?

### Injecting Domain Knowledge into AI
From a historical perspective, trying to “program” AI (i.e., inject domain knowledge) is nothing new—the change is that AI has never before been so powerful, nor such a “black box” in terms of interpretability and control.
{% figure %}
[<img class="postimage" src="{{ site.baseurl }}/assets/img/posts/2019-03-03-weak_supervision/ai_bg.png"/>]({{ site.baseurl }}//assets/img/posts/2019-03-03-weak_supervision/ai_bg.png)
{% endfigure %}

In the 1970’s and 1980's, the focus in AI was on **expert systems**, which combined manually-curated *knowledge bases* of facts and rules from domain experts with *inference engines* to apply them. 
<!-- The port of input was simple: just enter new facts or rules into the knowledge base. However, this very simplicity also belied the brittleness of these systems. Entering rules and facts by hand was neither sufficiently exhaustive nor scalable enough to handle the long-tail, high-dimensional data (e.g. text, images, speech, etc.) present in many real world applications. -->
In the 1990’s, machine learning began to take off as the vehicle for integrating knowledge into AI systems, promising to do so automatically from **labeled training data** in powerful and flexible ways. Classical (non-representation-learning) machine learning approaches generally had two ports of domain expert input. First, these models were generally of much lower complexity than modern ones, meaning that smaller amounts of hand-labeled data could be used. Second, these models relied on hand-engineered features, which provided a direct way to encode, modify, and interact with the model’s base representation of the data. However, feature engineering was and still is generally considered a task for ML experts, who often would spend entire PhDs crafting features for a particular task.

Enter **deep learning models**: due to their impressive ability to automatically learn representations across many domains and tasks, they have largely obviated the task of feature engineering. However, they are for the most part complete black boxes, with little control for the average developer other than labeling massive training sets. In many senses, they represent the opposite extreme of the brittle but easily-controllable rules of old expert systems—they are flexible but hard to control. This leads us back to our original question from a slightly different angle: How do we leverage our domain knowledge or task expertise to program modern deep learning models? Is there any way to combine the directness of the old rules-based expert systems with the flexibility and power of these modern machine learning methods?

## Code as Supervision: Training ML by Programming
[Snorkel](http://snorkel.stanford.edu) is a system we have built to support and explore this new type of interaction with ML. In Snorkel, we use no hand-labeled training data, but instead ask users to write labeling functions (LFs), black-box snippets of code which label subsets of unlabeled data. For example, suppose we were trying to train a machine learning model to extract mentions of adverse drug reactions from the scientific literature. To encode a heuristic about negation, for example, we could try writing the LF below:

{% figure %}
[<img class="postimage_75" src="{{ site.baseurl }}/assets/img/posts/2019-03-03-weak_supervision/snorkel_lf.png"/>]({{ site.baseurl }}//assets/img/posts/2019-03-03-weak_supervision/snorkel_lf.png)
{% endfigure %}

We could then use a set of such LFs to label training data for our machine learning model. Since labeling functions are just arbitrary snippets of code, they can encode arbitrary signals: patterns, heuristics, external data resources, noisy labels from crowd workers, weak classifiers, and more. And, as code, we can reap all the other associated benefits like modularity, reusability, debuggability. If our modeling goals change, for example, we can just tweak our labeling functions to quickly adapt!

{% figure %}
[<img class="postimage" src="{{ site.baseurl }}/assets/img/posts/2019-03-03-weak_supervision/dp.png"/>]({{ site.baseurl }}//assets/img/posts/2019-03-03-weak_supervision/dp.png)
{% endfigure %}

One problem, of course, is that the labeling functions will produce noisy outputs which may overlap and conflict, producing less-than-ideal training labels. In Snorkel, we de-noise these labels using our data programming approach, which comprises three steps:

* First, we apply the labeling functions to unlabeled data.
* Next, we use a generative model to learn the accuracies of the labeling functions without any labeled data, and weight their outputs accordingly. We can even learn the structure of their correlations automatically.
* Finally, the end output is a set of probabilistic training labels, which we can use to train a powerful, flexible discriminative model (such as a deep neural network) that will generalize beyond the signal expressed in our LFs.

This whole pipeline can be seen as providing a simple, robust, and model-agnostic approach to “programming” an ML model!

### Notes from Snorkel in the Wild!
{% figure %}
[<img class="postimage" src="{{ site.baseurl }}/assets/img/posts/2019-03-03-weak_supervision/snorkel_system.png"/>]({{ site.baseurl }}//assets/img/posts/2019-03-03-weak_supervision/snorkel_system.png)
{% endfigure %}

In our recent [paper on Snorkel](https://arxiv.org/abs/1711.10160), we find that in a variety of real-world applications, this new approach to interacting with modern machine learning models works very well! Some highlights include:

* In a user study, conducted as part of a [two-day workshop on Snorkel](http://mobilize.stanford.edu/events/snorkelworkshop2017/) hosted by the [Mobilize Center](http://mobilize.stanford.edu), we compared the productivity of teaching subject matter experts to use Snorkel, versus spending the equivalent time just hand-labeling data. We found that they were able to build models not only **2.8x faster** but also with **45.5% better predictive performance** on average.

* On two real-world text relation extraction tasks---in collaboration with researchers from Stanford, the U.S. Dept. of Veterans Affairs, and the U.S. Food and Drug Administration---and four other benchmark text and image tasks, we found that Snorkel leads to an average **132% improvement over baseline techniques**.

* We explored the novel tradeoff space of how to model the user-provided labeling functions, leading to a rule-based optimizer for accelerating iterative development cycles.

## Next Steps: Massively Multi-Task Weak Supervision
Various efforts in our lab are already underway to extend the weak supervision interaction model envisioned in Snorkel to other modalities such as [richly-formatted data](https://arxiv.org/abs/1703.05028) and [images](https://arxiv.org/abs/1709.02477), supervising tasks with [natural language](https://arxiv.org/abs/1805.03818) and [generating labeling functions automatically](http://www.vldb.org/pvldb/vol12/p223-varma.pdf)! On the technical front, we’re interested in both extending the core data programming model at the heart of Snorkel, making it easier to specify labeling functions with higher-level interfaces such as natural language, as well as combining with other types of weak supervision such as [data augmentation](https://arxiv.org/abs/1709.01643).

The increasing prevalence of multi-task learning (MTL) scenarios also invites the question: what happens when our noisy, possibly correlated label sources now label multiple related tasks? Can we benefit by modeling the supervision for these tasks jointly? We tackle these questions in a new multitask-aware version of Snorkel, [Snorkel MeTaL](https://arxiv.org/abs/1810.02840), which can support multi-task weak supervision sources that provide noisy labels for one or more related tasks. 

One example we consider is the setting of label sources with different granularities. For example, suppose we are aiming to train a fine-grained named entity recognition (NER) model to tag mentions of specific types of people and locations, and we have some noisy labels that are fine-grained—e.g. Labeling “Lawyer” vs. “Doctor” or “Bank” vs. “Hospital”—and some that are coarse-grained, e.g. labeling “Person” vs. “Location”. By representing these sources as labeling different hierarchically-related tasks, we can jointly model their accuracies, and reweight and combine their multi-task labels to create much cleaner, intelligently aggregated multi-task training data that improves the end MTL model performance.

{% figure %}
[<img class="postimage" src="{{ site.baseurl }}/assets/img/posts/2019-03-03-weak_supervision/mmtl.png"/>]({{ site.baseurl }}//assets/img/posts/2019-03-03-weak_supervision/mmtl.png)
{% endfigure %}

We believe that the most exciting aspects of building data management systems for MTL will revolve around handling what we refer to as the *massively* multi-task regime, where tens to hundreds of weakly-supervised (and thus highly dynamic) tasks interact in complex and varied ways. While most MTL work to date has considered tackling at most a handful of tasks, defined by static hand-labeled training sets, the world is quickly advancing to a state where organizations --- whether [large companies](https://arxiv.org/abs/1812.00417), academic labs, or online communities—maintain tens to hundreds of weakly-supervised, rapidly changing, and interdependent modeling tasks. Moreover, because these tasks are weakly supervised, developers can add, remove, or change tasks (i.e. training sets) in hours or days, rather than months or years, potentially necessitating retraining the entire model. 

In a [recent paper](http://cidrdb.org/cidr2019/papers/p58-ratner-cidr19.pdf), we outlined some initial thoughts in response to the above questions, envisioning a massively multi-task setting where MTL models effectively function as a central repository for training data that is weakly labeled by different developers, and then combined in a central “mother” multi-task model. Regardless of exact form factor, it is clear that there is lots of exciting progress for MTL techniques ahead—not just new model architectures, but also increasing [unification](https://hazyresearch.github.io/snorkel/blog/mtl_systems.html) with transfer learning approaches, new weakly-supervised approaches, and new software development and systems paradigms. 

We’ll be continuing to post our thoughts and code at [snorkel.stanford.edu](http://snorkel.stanford.edu) --- feedback is always welcome!

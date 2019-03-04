---
layout: post
title: "Weak Supervision"
short-summary: "Using weak supervision to efficiently label training data."
summary: "Using weak supervision, or high-level noisy sources of labels, to efficiently label training data."
thumbnail: "assets/img/posts/2019-02-26-beyond_local_pattern_matching/thumb.png"
author: Us
tags: [nlp,ws]
---

<!-- <figure>
    <p>
    <img class="postimagehalf" src="{{ site.baseurl }}/assets/img/posts/2019-02-26-beyond_local_pattern_matching/img1.png"/>
    <img class="postimagehalf" src="{{ site.baseurl }}/assets/img/posts/2019-02-26-beyond_local_pattern_matching/img2.png"/>
    <figcaption>
        Example search results from Google, as of the writing of this article.
    </figcaption>
    </p>
</figure> -->


<!-- {% figure %}
[<img class="postimage_75" src="{{ site.baseurl }}/assets/img/posts/2019-02-26-beyond_local_pattern_matching/img3.png"/>]({{ site.baseurl }}/assets/img/posts/2019-02-26-beyond_local_pattern_matching/img3.png)
{% endfigure %} -->

<!-- [^1]: Danqi Chen, Adam Fisch, Jason Weston, Antoine Bordes. Reading Wikipedia to Answer Open-Domain Questions. ACL 2017. -->

<!-- ### Outline
(1.i) WS (brief version of the survey)

(1.ii) Snorkel / etc

(2) WS ==> more tasks, moar quickly ==> MMTL

For (1), probably just splice (https://hazyresearch.github.io/snorkel/blog/snorkel_programming_training_data.html- the Snorkel stuff) + a bit from (https://hazyresearch.github.io/snorkel/blog/ws_blog_post.html- the high level chart and some refs)? -->

<!-- introduction and outline -->
In recent years, the real-world impact of machine learning has grown in leaps and bounds. In large part, this is due to the advent of deep learning models, which allow practitioners to get state-of-the-art scores on benchmark datasets without any hand-engineered features. Given the availability of multiple professional-quality open-source machine learning frameworks like TensorFlow and PyTorch, and an abundance of available state-of-the-art models, it can be argued that high-quality machine learning models are almost a commoditized resource now. There is a hidden catch, however: the reliance of these models on massive sets of hand-labeled training data.

These hand-labeled training sets are expensive and time-consuming to create–taking months or years for large benchmark sets, or when domain expertise is required–and cannot be practically repurposed for new objectives. In practice, these training sets have to be assembled, cleaned, and debugged—a prohibitively expensive and slow task, especially when domain expertise is required. Even more importantly, tasks often change and evolve in the real world. For example, labeling guidelines, granularities, or downstream use cases often change, necessitating re-labeling. For all these reasons, practitioners have increasingly been turning to **weaker forms of supervision**, such as heuristically generating training data with external knowledge bases, patterns or rules, or other classifiers. Essentially, these are all ways of programmatically generating training data—or, more succinctly, **programming training data.** 

We begin by reviewing areas of machine learning that are motivated by the problem of labeling training data and describe our research that addresses challenges of modeling and integrating a diverse set of supervision sources. We also discuss our vision for building data management systems for the massively multi-task regime with hundreds of weakly supervised, *dynamic*, tasks interact in complex and varied ways. Check out the [Snorkel blog](snorkel.stanford.edu) for detailed discussions of these topics and more!

## How to Get More Labeled Training Data? A Review
Many traditional lines of research in machine learning are similarly motivated by the insatiable appetite of modern machine learning models for labeled training data. We start by drawing the core distinction between these other approaches and weak supervision at a high-level: **weak supervision is about leveraging higher-level and/or noisier input from subject matter experts (SMEs).**

The problem is that this is expensive: for example, unlike grad students, radiologists don’t generally accept payment in burritos and free T-shirts! Thus, many well-studied lines of work in machine learning are motivated by the bottleneck of getting labeled training data:

* In **active learning**, the goal is to make use of subject matter experts more efficiently by having them label data points which are estimated to be most valuable to the model (for a good survey, see (Settles 2012)). Traditionally, applied to the standard supervised learning setting, this means selecting new data points to be labeled–for example, we might select mammograms that lie close to the current model decision boundary, and ask radiologists to label only these. However, we could also just ask for weaker supervision pertinent to these data points, in which case active learning is perfectly complementary with weak supervision; as one example of this, see (Druck, Settles, and McCallum 2009).


* In the **semi-supervised learning** setting, we have a small labeled training set and a much larger unlabeled data set. At a high level, we then use assumptions about smoothness, low dimensional structure, or distance metrics to leverage the unlabeled data (either as part of a generative model, as a regularizer for a discriminative model, or to learn a compact data representation); for a good survey see (Chapelle, Scholkopf, and Zien 2009). More recent methods use adversarial generative (Salimans et al. 2016), heuristic transformation models (Laine and Aila 2016), and other generative approaches to effectively help regularize decision boundaries. Broadly, rather than soliciting more input from subject matter experts, the idea in semi-supervised learning is to leverage domain- and task-agnostic assumptions to exploit the unlabeled data that is often cheaply available in large quantities.


* In the standard **transfer learning** setting, our goal is to take one or more models already trained on a different dataset and apply them to our dataset and task; for a good overview see (Pan and Yang 2010). For example, we might have a large training set for tumors in another part of the body, and classifiers trained on this set, and wish to apply these somehow to our mammography task. A common transfer learning approach in the deep learning community today is to “pre-train” a model on one large dataset, and then “fine-tune” it on the task of interest. Another related line of work is multi-task learning, where several tasks are learned jointly (Caruna 1993; Augenstein, Vlachos, and Maynard 2015). Some transfer learning approaches take one or more pre-trained models (potentially with some heuristic conditioning of when they are each applied) and use these to train a new model for the task of interest; in this case, we can actually consider transfer learning as a type of weak supervision.

The above paradigms potentially allow us to avoid asking our SME collaborators for additional training labels. But what if–either in addition, or instead–we could ask them for various types of higher-level, or otherwise less precise, forms of supervision, which would be faster and easier to provide? For example, what if our radiologists could spend an afternoon specifying a set of heuristics or other resources, that–if handled properly–could effectively replace thousands of training labels? This is the key practical motivation for weak supervision approaches, which we describe next.
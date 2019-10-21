---
layout: post
title: "Answering Complex Open-domain Questions at Scale"
thumbnail: "assets/img/posts/2019-10-21-answering-complex-questions/needle-haystack.png"
excerpt: "The NLP community has made great progress on open-domain QA, but our systems still struggle to answer complex open-domain questions in an large collection of text. We present an efficient and explainable method for enabling multi-step reasoning in these systems."
summary: "The NLP community has made great progress on open-domain QA, but our systems still struggle to answer complex open-domain questions in an large collection of text. We present an efficient and explainable method for enabling multi-step reasoning in these systems."
author: <a href='http://qipeng.me/'>Peng Qi</a> 
tags: [nlp, ai, qa]
---
*This post was originally on [Peng Qi's website](http://qipeng.me/blog/answering-complex-open-domain-questions-at-scale.html) and has been replicated here (with minor edits) with permission.*

> **TL;DR:** The NLP community has made great progress on open-domain question answering, but our systems still struggle to answer complex questions over a large collection of text. We present an efficient and explainable method for enabling multi-step reasoning in these systems.


From search engines to automatic question answering systems, natural language processing (NLP) systems have drastically improved our ability to access knowledge stored in text, saving us countless hours spent memorizing facts and looking things up.

<figure>
<img src="https://live.staticflickr.com/2129/2239767394_bbd6cab970_z.jpg" width="640" height="425" alt="Card Catalog" style="padding:0">
<figcaption>
    Who's old enough to remember these indexes and not just the search engine ones? <br>(Photo credit: <a href="https://www.flickr.com/photos/reedinglessons/2239767394/">Reeding Lessons @ Flickr (CC BY-SA-NC 2.0)</a>)
</figcaption>
</figure>

Today, whenever we have a question in mind, the answer is usually one Google/Bing search away. For instance, _"Which U.S. state is the largest by area?"_

<figure>
    <img src='{{ site.baseurl }}/assets/img/posts/2019-10-21-answering-complex-questions/google-alaska.png' width='90%'>
<figcaption>
Alaska! But also, great job, Google!
</figcaption>
</figure>

Other questions, however, are less straightforward. For example, _"Who was the first to demonstrate that GPS could be used to detect seismic waves?"_ Google isn't of much help if we were to directly type this question as a search query. On the other hand, the Internet’s encyclopedia, Wikipedia, does have an answer for us:

<figure>
    <img src='{{ site.baseurl }}/assets/img/posts/2019-10-21-answering-complex-questions/dr-larson.png' width='90%'>
<figcaption>
Thank you, Dr. Larson!
</figcaption>
</figure>

Wouldn’t it be nice if an NLP system could answer this question for us, without us having to find the article ourselves? This problem, called _open-domain question answering (open-domain QA)_, is an active area of NLP research.

### Background: Open-domain QA
Before diving into our new method for open-domain QA, let us first take a moment to understand the problem setup, challenges, and why existing solutions are not quite enough to answer complex questions.

#### Open-domain vs Closed-domain / Restricted-context
The first question answering systems built by NLP researchers, such as [BASEBALL](https://web.stanford.edu/class/linguist289/p219-green.pdf) and [LUNAR](https://www.semanticscholar.org/paper/Lunar-rocks-in-natural-english%3A-explorations-in-Woods/6390e2772c3359e4f3b5430423ac996473449ebb), were highly domain-specific. They were adept at answering questions about US baseball players over the period of one specific year and about lunar rocks brought back to Earth, but not terribly helpful beyond the domains they were built for. In other words, they are _closed-domain_.

Since then, researchers have moved towards tackling open-domain QA. In open-domain QA, the questions are not limited to predefined domains and domain knowledge; ideally, the system should be able to sift through a very large amount of text documents to find the answer for us.

Single-document open-domain QA (also known as _reading comprehension_) is one of the research areas seeing recent breakthroughs in natural language processing, where an NLP system is given a single document (or just a paragraph) that might contain the answer to a question, and is asked to answer the question based on this context. Take our Dr. Larson question for an example (_"Who was the first to demonstrate that GPS could be used to detect seismic waves?"_). A single-document QA system might be trained to answer this question given only the Wikipedia page _"Kristine M. Larson"_. This is the format of many popular question answering datasets used in the NLP community today, e.g., [SQuAD](https://rajpurkar.github.io/SQuAD-explorer/).

Question answering systems trained on SQuAD are able to generalize to answering questions about personal biographies.

<figure>
    <img src='{{ site.baseurl }}/assets/img/posts/2019-10-21-answering-complex-questions/bio-peng.png' width='90%'>
<figcaption>
Recent reading comprehension systems can answer our question, given appropriate context. Demo credit: <a href="https://demo.allennlp.org/reading-comprehension/OTk1OTky">AllenNLP</a>.
</figcaption>
</figure>

However, such systems cannot help us answer our question about Dr. Larson if we didn’t already know to look at her biography, which is quite limiting.

To solve this, researchers are developing question answering systems over large text collections. Instead of provided with the exact context necessary to answer the question, the system is required to sift through a collection of documents to arrive at the answer, much like how we search for answers on the web. This setting, called _open-context open-domain QA_, is much more challenging than reading comprehension. But, it is also a lot more useful when we have a question in mind but don’t really have a good idea where the answer might be from. The main challenge, besides those of restricted-context QA, is to narrow down the large collection of texts to a manageable amount with scalable approaches, such that we can run reading comprehension models to arrive at the answer.

#### Open-domain QA Systems
Inspired by the [series of question answering competitions at the Text REtrieval Conference](https://trec.nist.gov/data/qamain.html) (TREC), researchers in recent years have started to look into adapting powerful neural-network-based QA models to the open-domain task.

[Danqi Chen](https://www.cs.princeton.edu/~danqic/) and collaborators first combined traditional search engines with modern, neural question answering systems to attack this problem. Their approach to open-domain QA, named [DrQA](https://arxiv.org/pdf/1704.00051.pdf), is simple and powerful: given a question, the system uses it to search a collection of documents for context documents that may contain the answer. Then, this reduced context is the input to a reading comprehension system, which predicts the final answer.

<figure>
    <img src='{{ site.baseurl }}/assets/img/posts/2019-10-21-answering-complex-questions/drqa.png' width='90%'>
<figcaption>
    Illustration of Chen et al.'s "DrQA" model, which was presented at ACL 2017. Figure from the official <a href="https://github.com/facebookresearch/DrQA">Github repo</a>.
</figcaption>
</figure>

Most of the recent research in open-domain QA has largely followed this two-stage approach of retrieving and reading, with added features such as reranking (see, for example, [(Wang et al., 2018)](https://arxiv.org/abs/1709.00023)) and neural retrieval systems and better joint training (see, for example, [(Das et al., 2019)](https://openreview.net/pdf?id=HkfPSh05K7) and [(Lee et al., 2019)](https://arxiv.org/pdf/1906.00300.pdf)).

#### The Challenge of Complex Open-domain Questions
All systems that follow this retrieve-and-read paradigm are ill-equipped to handle complex questions. Let’s walk through an illustrative example of why that is together.

We all forget the names of celebrities from time to time. Suppose, one day, you are curious: _"What is the Aquaman actor's next movie?"_ To answer this question, you would probably first search for _"Aquaman"_ or _"the Aquaman actor"_ to find out who he/she is. Hopefully after scrolling through a few top search results, you will realize the answer is _"Jason Momoa"_, and then move on to finding out what his next movie is.

In this simple example, not all of the supporting evidence needed to answer the question can be readily retrieved from the question alone, i.e., there's a knowledge discovery problem to solve.[^2] This makes these questions difficult for retrieve-and-read open-domain QA systems, because there is usually some evidence that lack a strong semantic overlap with the question itself. Below is a sketch of the relations between the real-world entities that illustrate the multiple steps of reasoning required to answer this question.

<figure>
    <img src='{{ site.baseurl }}/assets/img/posts/2019-10-21-answering-complex-questions/jason-momoa.png' width='90%'>
<figcaption>
Reasoning required to answer the question "What is the Aquaman actor's next movie?". In this case, "Jason Momoa" is the missing link that connects the question to its answer, but cannot be easily retrieved based on the question.
</figcaption>
</figure>

One solution to this problem might be to train neural retrieval and reading comprehension models jointly to update queries to find more evidence (Das et al. (2019) set out to do just that). While this might also work in our setting, pretraining the neural retriever with distant supervision to promote documents that contain the answer string will likely fail  because of the missing semantic overlap between the question and all necessary documents. End-to-end training will also be prohibitively expensive, because the search space for queries beyond the first step of reasoning is enormous. Even if one manages to train a neural system to accomplish this task, the resulting system is probably very computationally inefficient and not very explainable.

So, can we build an open-domain QA system that is capable of handling complex, multi-step reasoning questions, and doing so in an efficient and explainable manner? We present such a system in our new EMNLP-IJCNLP paper -- [Answering Complex Open-domain Questions Through Iterative Query Generation][paper-link].

### Answering Complex Open-domain Questions

To introduce our system, we start with the overall strategy we use to address the problem of mutli-step reasoning in open-domain QA, before moving on to the dataset we evaluate our system on and experimental results.

#### Overall Strategy

As we have seen, retrieve-and-read systems can't efficiently handle complex open-domain questions that require multiple steps of reasoning, because (a) these questions require multiple supporting facts to answer, and (b) it is usually difficult to find all supporting facts necessary with only the question. Ideally, we want a system to be able to iterate between "reading" the information retrieved and finding further supporting evidence if necessary, just like a human.

That is exactly where the "iterative query generation" part of the paper title comes into play. We propose an open-domain QA system that iteratively generates natural language queries based on the currently retrieved context and retrieves more information if needed before finally answering the question. This allows us to (a) retrieve multiple supporting facts with different queries, and (b) make use of documents retrieved in previous iterations to generate queries that wouldn’t have been possible from the question alone. Moreover, because our system generates natural language queries, we can still leverage off-the-shelf information retrieval systems for efficient retrieval. Furthermore, the steps our model follows are more explainable to a human, and allow for human intervention at any time to correct its course.

Given the English Wikipedia as our source of textual knowledge, the full system operates as follows to answer the question _"Which novel by the author of 'Armada' will be adapted as a feature film by Steven Spielberg?"_:

<figure>
    <img src='{{ site.baseurl }}/assets/img/posts/2019-10-21-answering-complex-questions/golden-retriever.png' width='90%'>
<figcaption>
The proposed model answers the question "Which novel by the author of 'Armada' will be adapted as a feature film by Steven Spielberg?". The system first iterates between reading and retrieving to gather supporting facts, then concatenates all the top retrieval results and feeds them into a restricted-context QA model with the question to generate the final answer.
</figcaption>
</figure>

To answer this question, the model starts by generating a query to search Wikipedia to find information about the novel _Armada_. After "reading" the retrieved articles, it then attempts to search for _Ernest Cline_ (the name of the author) for more information. Finally, when we have retrieved all the context necessary to answer the question, we concatenate the top retrieved articles from these retrieval steps, and feed them into a restricted-context QA system to predict the final answer.

The main challenge in building this model lies in training the query generators collaboratively to generate useful queries for retrieving all the necessary information. Our main contribution is an efficient method for training these query generators with very limited supervision about which documents to retrieve, yielding a competitive system for answering complex and open-domain questions. Our method is based on the crucial observation that, if the question can be answered with knowledge from the corpus, then there exists a progressive chain (or graph) of reasoning we can follow. In other words, we note that at any given time in the process of finding all supporting facts, there is some strong semantic overlap between _what we already know_ (the question text, plus what we have found so far), and _what we are trying to find_ (the remaining supporting facts).


<figure>
	<img src='{{ site.baseurl }}/assets/img/posts/2019-10-21-answering-complex-questions/needle-haystack.png' width='90%'>
<figcaption>
Finding the multiple supporting facts necessary to answer complex questions is much like finding multiple needles in a haystack. Instead of looking for them independently, we make use of the thread connecting these needles, which is the strong semantic overlap between what we know and what we are trying to find.
</figcaption>
</figure>

In the beginning, the question the system is asked is all the information _we already know_. We are _trying to find_ any document part of reasoning chain needed to answer this question. Based on our observation, at least one of the gold documents[^4] would have strong semantic overlap with the question, and our goal is to find one such document to bootstrap our chain of reasoning. In our Armada example, this document would be the Wikipedia page of Armada the novel, where the overlap is the name _"Armada"_, and the fact that it’s a novel. To find this document with the help of a text-based information retrieval (IR) system, we just need to identify this overlap and use it as the search query.

After one step of retrieval, we have hopefully retrieved the _"Armada (novel)"_ page from Wikipedia, among others. If, at training time, we also know that the _"Ernest Cline"_ page is the next missing link in our chain of reasoning, we can apply the same technique. Now, the semantic overlap between _what we know_ (the question, the _"Armada (novel)"_ page, plus some other Wikipedia pages), and _what we are trying to find_ (_"Ernest Cline"_) to generate the desired query, _"Ernest Cline"_. To find this semantic overlap, we simply employ  longest common substring or longest common subsequence algorithms between _the knowns_ and _the wanted_.

With the desired queries at each step of reasoning, we can then train a model to predict them from the retrieval context (question + already retrieved documents) at each step. We then use these query generators to complete the task of open-domain multi-step reasoning. We cast the query generation problem as one of restricted-context QA, since the goal is to map the given question and (retrieved) context to some target derived from the context.

We name the full system GoldEn (Gold Entity) Retriever, because the model-retrieved Wikipedia pages are mostly entities, and it's a fun name for a retrieval-oriented model! Below are some example questions and the desired queries we train the query generators with:

<figure>
    <table>
        <tr>
            <th>Question</th>
            <th>Step 1 Query</th>
            <th>Step 2 Query</th>
        </tr>
        <tr>
            <td>What government position was held by the woman who portrayed Corliss Archer in the film Kiss and Tell?</td>
            <td>Corliss Archer in the film Kiss and Tell</td>
            <td>Shirley Temple</td>
        </tr>
        <tr>
            <td>Are Giuseppe Verdi and Ambroise Thomas both Opera composers?</td>
            <td>Giuseppe Verdi</td>
            <td>Ambroise Thomas</td>
        </tr>
    </table>
<figcaption>
Example queries generated from our overlap-finding process to train the query generators in GoldEn Retriever. As you can see in the first example, the query at Step 2 reveals information we can only find through iterative retrieval, and is not contained in the original question.
</figcaption>
</figure>

Two practical notes should be mentioned here. First, it is not difficult to see that our observation that supervision signal for query generation can be derived from this semantic overlap generalizes to any number of supporting documents. It also requires no additional knowledge about how the question can or should be decomposed into sub-questions to answer (which previous work has studied, e.g., [(Talmor and Berant, 2018)](https://www.aclweb.org/anthology/N18-1059) and [(Min et al., 2019)](https://www.aclweb.org/anthology/P19-1613)). As long as the gold supporting documents are known at training time, we can use this technique to construct the chain of reasoning in an open-domain setting very efficiently at scale. Second, we further make no assumption about knowledge of the order in which documents should be retrieved. At any given step of open-domain reasoning, one can enumerate all of the documents that have yet to be retrieved, find its semantic overlap with the retrieval context, and launch searches with these generated queries. Documents that are in the immediate next step of reasoning will naturally be more discoverable, and we can choose the desired queries accordingly. In our Armada example, for instance, the overlap between the question and the Ernest Cline article is _"Steven Spielberg"_, _"film"_, etc, which lead us nowhere close to the _"Ernest Cline"_ page, thus these are not chosen as the first-step query at training time.

#### Dataset: HotpotQA
To test the performance of GoldEn Retriever, we evaluate it on [HotpotQA](https://hotpotqa.github.io/), a recent multi-hop question answering dataset presented at EMNLP 2018 (by me &amp; collaborators). HotpotQA is a crowd-sourced QA dataset on English Wikipedia articles, in which crowd-workers are presented the introductory paragraphs from two related Wikipedia articles and asked to generate questions that require reasoning with both paragraphs to answer. Our example question about the Armada novel is one such question from this dataset. To encourage the development of explainable QA systems, we also asked crowd workers to highlight the sentences from these paragraphs that support their answer (we call these "supporting facts"), and ask QA systems to predict them at test time.

HotpotQA features two evaluation settings: a few-document distractor setting, and an open-domain fullwiki setting, which we focus on, where the system is only given the question and the entire Wikipedia to predict the answer from. HotpotQA also features a diverse range of reasoning strategies, including questions involving missing entities (our Armada example, where Ernest Cline is not in the question), intersection questions (_What satisfies property A and property B?_), and comparison questions, where two entities are compared by a common attribute, among others.[^3]

QA systems on this dataset are evaluated on two aspects, answer accuracy and explainability. Answer accuracy is evaluated with answer exact matches (EM) and unigram F1, and explainability is similarly evaluated with EM and F1 by calculating the supporting fact overlap between predictions and annotations. These two aspects are unified by joint EM and F1 metrics, which encourage QA systems to work well on both.

For the final restricted-context QA component, we use a modified BiDAF++ model in this work. For more technical details, please refer to [our paper][paper-link].


#### Results

We evaluate the effectiveness of our GoldEn Retriever model on two aspects: its performance on retrieving the gold supporting documents, and it’s end-to-end performance in question answering.

For retrieval performance, we compare GoldEn Retriever to a retrieve-and-read QA system that just retrieves once with the question. We evaluate these approaches on the recall of the two gold paragraphs when a total of 10 paragraphs are retrieved by each system, because this metric reflects the ceiling performance of the entire QA system if the restricted-context QA component were perfect.


<figure>
    <img src='{{ site.baseurl }}/assets/img/posts/2019-10-21-answering-complex-questions/ir-recall.png' width='90%'>
<figcaption>
Retrieval performance of a retrieve-and-read system vs GoldEn Retriever on the gold paragraphs.
</figcaption>
</figure>

As can be seen from the figure, while both systems achieve decent recall on the paragraph that is usually more connected to the question ("Paragraph 1" in the figure), GoldEn Retriever obtains significant improvement through iterative retrieval with query generation on the other paragraph (~24% improvement). This means for about 24% of the questions, GoldEn Retriever is able to find both gold supporting documents while the retrieve-and-read system can't. Further analysis shows that this is mainly from the improved recall for non-comparison questions (for which recall improved by about 25%), where the retrieval problem is less trivial.

For end-to-end QA performance, we compare GoldEn Retriever against various retrieve-and-read baselines on the development set, as well as systems submitted to the public leaderboard on the hidden test set. 

<figure>
    <img src='{{ site.baseurl }}/assets/img/posts/2019-10-21-answering-complex-questions/fullwiki-joint-f1.png' width='90%'>
<figcaption>
Comparing GoldEn Retriever against various other systems on HotpotQA's fullwiki setting.
</figcaption>
</figure>

We first contrast the performance of the QA component when using the IR system originally used in HotpotQA (as reflected by the released fullwiki dev set) and Elasticsearch in an retrieve-and-read setting. As can be seen by the leftmost two bars in the figure, a better search engine does improve end-to-end performance (from 22.75% F1 to 27.11%). However, this is still far from the best previously published system (34.92% F1 on the test set, which is empirically &plusmn;2% from the model’s dev set performance). With GoldEn Retriever, we improve this state of the art to 39.13% F1, which is significant especially if one considers that the previous state-of-the-art model uses BERT and we don’t. Although this doesn’t match the contemporaneous work which achieves 47.6% F1 with another BERT-based model, we see that if our query generators were able to faithfully reproduce the desired queries on the dev set, the performance of our system wouldn’t have been far off ("Oracle IR").

For explainability, aside from reporting supporting fact metrics that are part of HotpotQA's evaluation, we can also look at the search queries GoldEn Retriever generates on the dev set. As can be seen in the example below, the natural language queries generated by the model are very understandable. Furthermore, one can see where the model is making mistakes and correct it in the system if needed. 

<figure>
    <table>
        <tr>
            <th>Question</th>
            <th>Step 1 Predicted</th>
            <th>Step 2 Predicted</th>
        </tr>
        <tr>
            <td>What video game character did the voice actress in the animated film Alpha and Omega voice?</td>
            <td>voice actress in the animated film Alpha and Omega <span style="font-style: italic; color: #44aa33;">(animated film Alpha and Omega voice)</span></td>
            <td>Hayden Panettiere</td>
        </tr>
        <tr>
            <td>Yau Ma Tei North is a district of a city with how many citizens?</td>
            <td>Yau Ma Tei North</td>
            <td>Yau Tsim Mong District of Hong Kong <span style="font-style: italic; color: #44aa33;">(Hong Kong)</span></td>
        </tr>
    </table>
<figcaption>
Examples of queries generated by GoldEn Retriever on dev set examples. The model-generated queries are shown in black, and the heuristic-generated "desired queries" are shown in parenthesis in <span style="font-style: italic; color: #44aa33;">green italic font</span> when they differ from the model-generated ones. In the first example, we see that the model actually generates a constituent whereas the heuristics largely ignores constituency structure; in the second example, however, the model generated a Step 2 query that is overly specific.
</figcaption>
</figure>

### Resources
To help facilitate future research in open-domain multi-step reasoning, we make the following resources publicly available:

* The code to reproduce our results and our pretrained models
* Generated “desired” query files and modified HotpotQA training and development files generated from the heuristics to train GoldEn Retriever models
* Predicted search queries and dev/test set input for our restricted-context QA model

All of these can be found in our [code repository on GitHub](https://github.com/qipeng/golden-retriever).

**Language Note:** All datasets and most of the research mentioned in this post are collected/tested for the English language only, but our principle of semantic overlap is applicable to answering open-domain complex questions in other languages than English if suitably augmented with lemmatization for highly inflected languages.


[^2]: This is of course contingent on the fact that very few highly ranked articles on the Web mention Jason Momoa in his next movie in close proximity to stating that he’s the “Aquaman” star who played Aquaman in that movie. This is just an example to demonstrate that as simple as this question seems, it’s not too difficult to construct questions that require information from more than one document to answer.
[^3]: Comparison questions make up about 25% of the HotpotQA dataset. For more details please see [our HotpotQA paper](https://arxiv.org/pdf/1809.09600.pdf).
[^4]: By "gold documents" we mean the documents needed in the chain of reasoning to answer the question.

[paper-link]: https://nlp.stanford.edu/pubs/qi2019answering.pdf


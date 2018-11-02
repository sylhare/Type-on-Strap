---
layout: post
title: "Towards Vision-Based Smart Hospitals"
subtitle: "Building a computer vision system to detect when people wash their hands"
feature-img: "assets/img/posts/2018-11-02-measuring-hand-hygiene-in-hospitals/main.png"
thumbnail: "assets/img/posts/2018-11-02-measuring-hand-hygiene-in-hospitals/thumb.png"
author: Albert Haque and Michelle Guo
tags: [Test, Lorem]
---

Every year, ​more people​ die from hospital-acquired infections than ​car accidents​. This
means when you are admitted to a hospital, there is a ​1 in 30​ chance you will become
worse than had you not gone to the hospital at all.

A dire situation, but one we can improve through better hygiene. Hand hygiene is the first
line of defense in preventing the spread of infections not only in hospitals, but also in
public spaces like airports and restaurants. This is already well known, so the issue is not
one of ignorance but of vigilance; automated verification techniques are needed to keep
track of and check for hand washing. Among the many technological solutions to this
problem, perhaps the simplest is to adopt the most common human strategy to visually
confirm whether people are washing their hands or not, with computer vision.

_(Left) Camera above a hand hygiene dispenser in our computer science building. (Right)
Depth sensor in a hospital, above a dispenser, outside a patient room._

In this post, we will talk about the ​Stanford Partnership in AI-Assisted Care (PAC)​ project
to detect and measure proper hand hygiene. This multi-year project was done with many
collaborators around the world over several conference, workshop, and journal papers.
Much remains to be done, but we hope this technology can help hospitals decrease
infection rates and improve patients’ health.


## Motivation

Today, hospitals reinforce proper hand hygiene through educational tools such as
medical school classes, flyers posted on
bulletin boards, and weekly staff meetings. The
World Health Organization​ has even
proposed the ​Five Moments​ of hand hygiene.
It explicitly defines when a healthcare worker
should wash their hands.
They have also started to measure the
extent to which employees follow proper
hand hygiene, typically through RFID cards
or badges worn by employees. While this
works to some extent, there are workflow
disruptions such as swiping the RFID card by
the soap dispenser when rushing into a new
room. This stems from technical reasons: normal RFIDs have short range, and ‘active’
RFIDs with longer range are constrained by their directional antenna and need batteries.
Clearly, a new solution that does not have the drawbacks of RFID techniques is needed.

## Computer Vision and Hospitals

We have been developing such a technology over the past few years. As a collaboration
between the Stanford AI Lab, Stanford University School of Medicine, and Stanford
Healthcare, our team has made it possible to track hand hygiene in hospitals using
computer vision. Our solution does not have the drawbacks of RFID. Computer vision
works from tens of meters away, does not require clinicians to wear additional tags or
badges, and can preserve privacy by recording only depth data from which the identities

##### of people complying or failing to comply with hand hygiene rules cannot be known.

In healthcare, computer vision is already being used for medical imaging taskswith
researchers developing machine learning models in the comfort of their own lab using
pre-collected datasets. But, there hasn’t been much work on computer vision in the
physical hospital space. Luckily, computer vision has been applied in physical spaces for


another problem: self-driving cars. Self-driving cars use tons of sensors to understand
their environment: lidar, radar, sonar, and of course, regular cameras. Can we use some
of these sensors inside hospitals to better understand the ​ _healthcare ​_ environment?

### Depth Sensors

Depth sensors (e.g., Xbox Kinects) are like normal cameras but instead of recording color,
they record distance. In a normal color image, each pixel denotes a color. In a depth
image, each pixel denotes the “distance” to the pixel in real-world space. This is usually a
floating point number such as 1.337 meters.
_(Left) Color photo of the hospital, taken with a cell phone. (Right) Corresponding depth
image taken by our sensor on the ceiling. Darker colors indicate objects closer to the depth
sensor._
In the depth image above, notice how you can’t really see the person’s face. However, you
can still tell what they’re doing: There is a group of people watching a person reach for
the dispenser. From the depth image, we don’t know who specifically these people are.
This is important for privacy, especially in hospitals.
To develop and validate our computer vision technology, we installed depth sensors on
the ceiling at two hospitals. One was a children’s cardiovascular unit and the other was
an adult intensive care unit (ICU). For context, the ICU houses patients with
life-threatening conditions. These patients need constant monitoring, require specialized
equipment, and are given powerful medications just to stay alive. It is ​estimated​ that
0.66% of the United States GDP is spent in ICUs alone.


_(Above) Our depth sensors installed on the ceiling of a children’s hospital._
With depth sensors installed at two different hospitals, we then used 3D computer vision
tools to automatically measure hand hygiene. This involves three steps. First, we need to
detect the healthcare staff. Second, we need to track them as they walk around the unit.
Third, we need to classify their hand hygiene behavior. These three steps will allow us to
compute quantitative hand hygiene compliance metrics.

### Pedestrian Detection

Continuing with our self-driving car analogy: one of the first steps for understanding the
environment is to detect people. Although deep learning object detection methods would
clearly work here, most of methods are developed for color RGB images. So, we chose to
use a previously devised approach^1 that can run on any type of image by leveraging two
aspects of the problem: that people typically occupy small amounts of space in a given
image of a room, and that in depth images people typically look like ‘blobs’ that clearly
stand out against the background of the floor.

(^1) For more details, see ​“SCOOP: A Real-Time Sparsity Driven People Localization Algorithm”


One way to detect people is to determine an ​occupancy grid​ $$x$$ over the ground. An
occupancy grid is a binary matrix indicating whether a person is occupying a particular
position in the ground plane. Let $$y$$ denote our observation vector (i.e., the binary
foreground silhouettes of people). It is possible to compute $$y$$ using object
segmentation methods.
(Above) Entries of the dictionary $$ D $$. Each dictionary entry contains a synthetic
image, corresponding to how a person would look like, had they been standing at that
position.
To make the problem simpler, we can convert the ground (i.e., floor of a room) into a
discrete grid. For every point in this grid, we can “imagine” a person at that position by
rendering a blob, roughly the same height as a person. We can create a dictionary $$D$$
containing the blobs at every single point on the ground plane (remember: because we
synthetically created these blobs, we know their exact 2D and 3D position).
The problem reduces to solving the following optimization problem. We wish to solve for
the final occupancy grid $$ x $$:


$$ x = \argmin_{x\in\{0,1\}} ||y - D x||_2^2 \text{\; s.t. \;} ||x||_0 < \epsilon_p $$
Where $$ \epsilon_p $$ is a user-selected upper bound of the sparsity level. In summary,
the inputs to the above solver are: (i) an image $$y$$ containing segmentations of people
and (ii) a dictionary $$D$$ of synthetic blobs representing a person at every possible
discrete position. The output of the solver is a binary matrix (or vector) $$x$$ containing
1’s corresponding to positions of people in the input image $$y$$.

##### Tracking Across The Hospital Unit

In order to build a truly smart hospital, we need to use sensors spread out across the
entire hospital unit. Not everything happens in front of one sensor. However, this means
we need computer vision algorithms that can track people across sensors. Not only can
this provide details for hand hygiene compliance, but it can also be used for workflow and
spatial analytics. Formally, we want to find the set of trajectories $$ X $$, where each
trajectory $$ x \in X $$ is represented as an ordered set of detections, $$ L_x = (l_x^{(1)}
,...,l_x^{(n)} ) $$, representing the detected coordinates of pedestrians. The problem can
be written as a ​maximum a-posteriori (MAP) estimation​ problem.
Next, we assume a Markov-chain model connecting every intermediate detection $$
l_x^{(i)} $$ in the trajectory $$ X $$, to the subsequent detection $$ l_x^{(i+1)} $$ with a
probability given by $$ P(l_x^{(i+1)} | l_x^{i}) $$. We can now formulate the MAP task as a
linear integer program by finding the flow $$ f $$ that minimizes the cost $$ C $$:
$$ \min_{f} & \quad C = \sum_{x_i\in X}\alpha_i f_i + \sum_{x_i, x_j \in X}\beta_{ij} f_{ij}
\qquad \text{s.t} \qquad f_i, f_{ij} \in (0,1) \text{\quad and \quad} f_{i} = \sum_{j}{f_{ij}} $$
where $$ f_i $$ is the flow variable indicating whether the corresponding detection is a
true positive, and $$ f_{ij} $$ indicates if the corresponding detections are linked together.
The variable $$ \beta_{ij} $$ denotes the transition cost given by $$ \log P (l_i|l_j) $$ for
the detection $$ l_i,l_j\in L $$. The local cost $$ \alpha_i $$ is the log-likelihood of an
intermediate detection being a true positive. For simplicity, we assume that all detections
have the same likelihood. This is equivalent to the flow optimization problem, solvable in
real-time with ​k-shortest paths​.


### Hand Hygiene Activity Classification

So far, we have identified the tracks (i.e., position on the global hospital-unit ground
plane) of all pedestrians in the unit. The last step is to detect hand hygiene activity and
link it to a specific track. Hand hygiene activity is defined as ​ _positive​_ when a person uses
a wall-mounted alcohol-based gel dispenser. We then label each pedestrian track as
_clean​_ or ​ _not clean​_.
Deployment of sensors in real-world settings is often prone to installation constraints.
Whether intentional or not, construction and maintenance technicians install sensors that
vary in both angle and position. If our goal is to propose a non-intrusive vision-based
system for tracking hand hygiene, our model must be robust to such variances. We
cannot foresee every possible viewpoint and thus wish to train a single model to work
with any sensor viewpoint. However, traditional ​convolutional neural network​ (CNNs) are
generally ​not viewpoint invariant​ and so we use a ​spatial transformer network​ (STN) to
solve this.
_(Left) Data augmentation stage with a person segmentation. (Right) Hand hygiene activity
classifier: a ​spatial transformer network​ plus a ​densely connected​ convolutional network._
The input to the STN is any arbitrary image and the output is a warped image. To help our
model learn more quickly, we also provide a person segmentation (i.e., body mask) to the
STN. This body mask can be extracted using classical foreground-background
techniques or deep learning approaches. The STN warps the image into a learned,
“viewpoint-invariant” form. From this warped image, we use a standard CNN (i.e.,
DenseNet​) to perform binary classification of whether someone used the hand hygiene
dispenser or not.


### Spatio-Temporal Matching

At this point, we have a set of tracks and a separate set of hand hygiene detections. We
still need to combine them, which is tricky since this introduces two variables: space and
time.
For each hand hygiene classifier detection (i.e., dispenser is being used), we must match
it to a single track. A match occurs between the classifier and tracker when a track $$
\mathcal{T} $$ satisfies two conditions:

1. Track $$ \mathcal{T} $$ contains $$ (x,y) $$ points $$ \mathcal{P} $$ which
    occur at the same time as the hand hygiene detection event $$ \mathcal{E} $$,
    within some temporal tolerance level.
2. At least one point $$ p \in \mathcal{P} $$ is physically nearby the sensor
    responsible for the detection event $$ \mathcal{E} $$. This is defined by a
    proximity threshold around the patient's door.
If there are multiple tracks that satisfy these requirements, we break ties by selecting the
track with the closest $$ (x, y) $$ position to the door. The final output of our model is a
list $$ T $$ of tracks, where each track consists of an ordered list of $$ (t, x, y, a) $$
tuples where $$ t $$ denotes the timestamp, $$ x, y $$ denote the 2D ground plane
coordinate, and $$ a $$ denotes the latest action or event label. From $$ T $$, we can
compute the compliance rate or compare with the ground truth for evaluation metrics.

## Comparison to Human Auditors & RFID

Today, many hospitals measure hand hygiene compliance using ​ _secret shoppers​_. Secret
shoppers are trained individuals who walk around hospital units and watch if staff wash
their hands. They are called ​ _secret​_ because they do not reveal the true purpose of their
visit. A secret shopper could be a nurse, doctor, or even a visitor. We refer to this as a
_covert​_ observation, as opposed to an ​ _overt​_ observation performed by someone openly
disclosing their audit. The purpose of covert observations is to minimize the ​Hawthorne
effect​ (i.e., you change your behavior because someone is watching you).


(Above) Floor plan of the hospital ICU. We sent three undercover auditors to monitor
hand hygiene for two hours.
We conducted two covert observational studies: (i) a single auditor responsible for
monitoring the entire unit and (ii) a group of three auditors with a collective responsibility
of the entire unit. The reason for having two different covert groups is that it allows us to
see how many people you’d need to match the performance of our algorithm. Both covert
groups were disguised as hospital visitors. The group of three was spread out over the
unit, remaining stationary, while the individual auditor constantly walked around the unit
while monitoring hand hygiene compliance.
An important baseline we analyzed was the RFID method. In some sense, RFID can be
interpreted as a proximity algorithm. If a healthcare worker approaches a radio base
station (e.g., mounted on a hand hygiene dispenser) within some threshold, the RFID will
activate, indicating a hand hygiene event. The problem here is the healthcare worker may
not have actually washed their hands -- they were simply near the dispenser. Due to its
simplicity, we can simulate RFID: if a person approaches within one meter of a dispenser^2
, we consider the person clean. We have tracks and we know the position of dispensers --
therefore we can easily compute RFID detectors.

### Results

As you would suspect, RFID produces a lot of false positives. It got a compliance
accuracy of 18%. That is, it has a 18% chance of correctly predicting a track as clean or

(^2) ​A realistic distance for RFID activation, although less sensitive more near range badges
are also possible


dirty. One human auditor did much better at 63%, and three people did better yet at 72%.
However, our algorithm still does better and got 75% accuracy. This is not too surprising,
since the auditors are competing with the “global view” computer vision system.
Many people wonder: ground truth (the 100% correct tracking of hand washing) is usually
annotated by humans. How do the humans observers get less than 100%? And the
answer is: ground truth in our experiments was annotated remotely and not in real-time.
Remote annotators had access to all sensors and could play the video forward and
backward in time to ensure their annotations are correct. In-person auditors did not have
“access” to all sensors and they could not replay events in time.
(Above) Hand hygiene detections over time. Blue squares indicate someone using a hand
hygiene dispenser. Darker blue indicates more simultaneous events. The ground truth is
shown at the bottom. In general, more white space is bad.
Numbers aside, a more interesting result is a visual one. The image above shows how
infrequently the in-person auditors detect hand hygiene activity. Notice all the white
space? If you look at the ground truth row, there is usually no white space. That means
the observers are missing a lot of hand hygiene events. This is often due to observers
being distracted: they may doze off, may look at unrelated activity elsewhere in the unit,
or simply just not see hand hygiene events occuring.


(Above) Spatio-temporal heatmap of people walking in the intensive care unit. Yellow/red
colors indicate more people standing/walking in that area.
We conclude with one final visualization. The animation above shows a top view of the
hospital unit. Because we can track people across the entire unit, we know their specific
(x,y,z) position pretty much all the time. We plotted each point and created a heatmap
over time. This type of ​ _spatial analytics​_ can be useful for identifying traffic patterns and
potentially trace the spread of disease. Areas that are always yellow/red indicate crowded
spaces. These spaces are usually at hallway intersections or immediately outside patient
rooms. If you look carefully, you can spot our stationary auditors in red.

## Future Directions

We’ve shown how computer vision and deep learning can be used to automatically
monitor hand hygiene in hospitals. At the ​Stanford Partnership in AI-Assisted Care​, hand
hygiene is just one use case of computer vision in healthcare. We are also developing
computer vision systems to ​monitor patient mobility levels​, analyze the ​quality of surgical
procedures​, and check for ​anomalies in senior living​. We hope this work sheds new light
on the potential and impact of AI-assisted healthcare.

## References

Viewpoint Invariant Convolutional Networks for Identifying Risky Hand Hygiene
Scenarios.​ M. Guo, A. Haque, S. Yeung, J. Jopling, L. Downing, A. Alahi, B. Campbell, K.
Deru, W. Beninati, A. Milstein, L. Fei-Fei. ​ _Workshop on Machine Learning for Health (ML4H),
Neural Information Processing Systems (NIPS), Long Beach, CA, December 2017._


Towards Vision-Based Smart Hospitals: A System for Tracking and Monitoring Hand
Hygiene Compliance.​ A. Haque, M. Guo, A. Alahi, S. Yeung, Z. Luo, A. Rege, A. Singh, J.
Jopling, L. Downing, W. Beninati, T. Platchek, A. Milstein, L. Fei-Fei. ​ _Machine Learning in
Healthcare Conference (MLHC), Boston, MA, USA, August 2017._
Vision-Based Hand Hygiene Monitoring in Hospitals.​ **​** S. Yeung, A. Alahi, Z. Luo, B. Peng,
A. Haque, A. Singh, T. Platchek, A. Milstein, L. Fei-Fei. ​ _American Medical Informatics
Association (AMIA) Annual Symposium, Washington, DC, USA, November 2016._




---
layout: post
title: "Gitlab improvements: caches and Docker"
description: This post describe how I can manage to reduce the time consumed by Gitlab pipelines by caching the dependencies and using Docker images
author-id: "galera"
categories: [devops, ci/cd]
tags: [devops, ci/cd, gitlab, react, docker,frontend]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/docker-gitlab/featured.jpg"
thumbnail: "assets/img/posts/docker-gitlab/featured.jpg"
image: "assets/img/posts/docker-gitlab/featured.jpg"
---

Gitlab or any other CI/CD system works really great to have an automated build system. However, you can waste lots of time if you don't think about carefully. Every job needs to download the build tools and dependencies, so that's a lot of time that could be reduced.

In this article I describe two techniques to avoid wasting that much time. First one is using a cache for the dependencies and the second is using a pre-built Docker image with the required build tools.

<p><!--more--></p>

## Scenario

We will be building a simple react application to be deployed into an S3 Bucket. This simple react application contains a bunch of npm dependencies to simulate a real application. You can find the code here: <a href="https://gitlab.com/adrian.galera/gitlab-docker-react">https://gitlab.com/adrian.galera/gitlab-docker-react</a>.

To simulate a real work environment, let's define a pipeline consisting in three steps:

- test: runs the tests defined in the project
- build: produce the artifact to be deployed
- deploy: deploy the artifact to an AWS S3 Bucket

## Basic pipeline

In this basic pipeline, no cache nor Docker is configured. Each job needs to download everything:

```yaml
test frontend:
  image: node:12.13-alpine
  stage: test
  script:
    - npm install
    - npm test

build frontend:
  image: node:12.13-alpine
  stage: build
  script:
    - npm install
    - npm build
  only:
    - master

deploy frontend:
  image: node:12.13-alpine
  stage: deploy
  script:
    - apk add python curl build-base zip
    - curl "https://s3.amazonaws.com/aws-cli/awscli-bundle.zip" -o "awscli-bundle.zip"
    - unzip -o awscli-bundle.zip
    - ./awscli-bundle/install -i /usr/local/aws -b /usr/local/bin/aws
    - aws --version
    - aws s3 sync build s3://random-bucket
  only:
    - master
```

`test` and `build` jobs use a nodejs Docker image and install all the dependencies and run the scripts. `deploy` job uses the same image and downloads and install the awsclient in order to upload the built artifact to a S3 bucket.

Next step is to measure the time:

| Test | Build | Deploy | Total |
| 1m 24s | 1m 26s | 57 s | 3m 48 s |

That's a lot of time wasted downloading dependencies or installing tools. In the next section we will reduce it by using caches

## Caching node_modules

We can use Gitlab cache feature to keep the contents of the node_modules folder instead of downloading them every time. It's really simple to set it up:

```yaml
stages:
  - test
  - build
  - deploy

test frontend:
  image: node:12.13-alpine
  stage: test
  script:
    - npm ci --cache .npm --prefer-offline
    - npm test
  cache:
    key: "node-modules"
    paths:
      - .npm/

build frontend:
  image: node:12.13-alpine
  stage: build
  script:
    - npm ci --cache .npm --prefer-offline
    - npm build
  only:
    - master
  cache:
    key: "node-modules"
    paths:
      - .npm/

deploy frontend:
  image: node:12.13-alpine
  stage: deploy
  script:
    - apk add python curl build-base zip
    - curl "https://s3.amazonaws.com/aws-cli/awscli-bundle.zip" -o "awscli-bundle.zip"
    - unzip -o awscli-bundle.zip
    - ./awscli-bundle/install -i /usr/local/aws -b /usr/local/bin/aws
    - aws --version
    - aws s3 sync build s3://random-bucket
  only:
    - master
  cache:
    key: "node-modules"
    paths:
      - .npm/
```

We only need to tell `npm` to work with the cache with this line:

`npm ci --cache .npm --prefer-offline` instead of the traditional `npm install`

and configure the cache in each job:

```yaml
cache:
  key: "node-modules"
  paths:
    - .npm/
```

We are telling gitlab to store the `.npm` folder in the cache named `node-modules`

The time improvement can observed in the following table:

| Test | Build | Deploy | Total |
| 1m 3s | 57s | 56s | 2m 54s |

We can observe a huge decrease in the `test` and `build` jobs. The time spent on deploy is pretty similar to the one before.

## Using Docker image

The `deploy` job is downloading and installing the awsclient to upload the built artifact to S3 bucket. We could improve that time by using a Docker image which already contains
the awslcient. In order to do so, we can build our own image and store it in Gitlab Container Registry. The Container registry can be found in the side bar: `Packages & Resgistries -> Container Registry`

We will build an image with nodejs 12:13 and the aws client. In order to do so, we will use the following Dockerfile:

```
FROM node:12.13-alpine
RUN apk add python curl build-base zip
RUN curl "https://s3.amazonaws.com/aws-cli/awscli-bundle.zip" -o "awscli-bundle.zip"
RUN unzip -o awscli-bundle.zip
RUN ./awscli-bundle/install -i /usr/local/aws -b /usr/local/bin/aws
RUN aws --version
```

Basically, the Dockerfile is running the same commands as the pipeline, but only once. Once is the image is built, we push it to the registry and we start using it. No need to install anything!

To build and push the image, you only need to run the following commands:

```bash
docker login registry.gitlab.com
docker build -t registry.gitlab.com/<user>/<project>/node_python_aws .
docker push registry.gitlab.com/<user>/<project>/node_python_aws
```

Once is pushed we can use it in the pipeline in combination with the cache:

```yaml
stages:
  - test
  - build
  - deploy

test frontend:
  image: registry.gitlab.com/<user>/<project>/node_python_aws
  stage: test
  script:
    - npm ci --cache .npm --prefer-offline
    - npm test
  cache:
    key: "node-modules"
    paths:
      - .npm/

build frontend:
  image: registry.gitlab.com/<user>/<project>/node_python_aws
  stage: build
  script:
    - npm ci --cache .npm --prefer-offline
    - npm build
  only:
    - master
  cache:
    key: "node-modules"
    paths:
      - .npm/

deploy frontend:
  image: registry.gitlab.com/<user>/<project>/node_python_aws
  stage: deploy
  script:
    - aws s3 sync build s3://random-bucket
  only:
    - master
  cache:
    key: "node-modules"
    paths:
      - .npm/
```

Note that in the `deploy` job, we're only using the aws command directly.

We can notice and decrease in time on the deploy pipeline, as we were expecting:

| Test | Build | Deploy | Total |
| 1m 9s | 1m 2s | 35s | 2m 46s |

The decrease of time, might not look like it's very big, but that's because that's a very simple example.

For bigger projects, the improvements can be huge, since there are lots of dependencies and build tools that could be cached or pre-installed

## Wrap up

We have seen different techniques to speedup Gitlab pipeline: using caches for dependencies and Docker image to pre-install build tools.

You can find the comparisson tables of the different approaches hereunder:

Improvement | Test | Build | Deploy | Total |
None | 1m 24s| 1m 26s | 57s | 3m 48s |
Cache | 1m 3s | 57s | 56s | 2m 54s |
Cache + Docker | 1m 9s | 1m 2s | 35s | 2m 46s |

### Total

<div style="height: 400px">
<canvas id="total-results"></canvas>
</div>

### Jobs

<div class="results-job">
    <div style="height: 400px">
    <canvas id="test-results"></canvas>
    </div>
    <div style="height: 400px">
    <canvas id="build-results"></canvas>
    </div>
    <div style="height: 400px">
    <canvas id="deploy-results"></canvas>
    </div>
</div>

In the charts separated by job, we can see that the time improvement for test and build comes from using caches. Regarding deploy job, the big improvement comes we have pre-installed aws client provided by the Docker image.

<script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/2.9.3/Chart.bundle.min.js"></script>
<script>
new Chart( "total-results", {
    type: 'bar',
    data: {
        labels: ["Simple", "Cached", "Cached + Docker"],
        datasets : [
            {
                label: "Pipeline total elapsed time (s)",
                data: [228,174,166],
                "fill": false,
                "backgroundColor": [
                    "rgba(255, 99, 132, 0.2)",
                    "rgba(75, 192, 192, 0.2)",
                    "rgba(153, 102, 255, 0.2)",
                ],
                "borderColor": [
                    "rgb(255, 99, 132)",
                    "rgb(75, 192, 192)",
                    "rgb(153, 102, 255)",
                ],
                "borderWidth": 1
            }
        ]
    },
    options: {
        responsive: true,
        maintainAspectRatio: false,
        scales: {
            yAxes : [{
                "ticks" : {
                    "beginAtZero" : true
                }
            }]
        }
    }
});
new Chart( "test-results", {
    type: 'bar',
    data: {
        labels: ["Simple", "Cached", "Cached + Docker"],
        datasets : [
            {
                label: "Test job elapsed time (s)",
                data: [84,63,69],
                "fill": false,
                "backgroundColor": [
                    "rgba(255, 99, 132, 0.2)",
                    "rgba(75, 192, 192, 0.2)",
                    "rgba(153, 102, 255, 0.2)",
                ],
                "borderColor": [
                    "rgb(255, 99, 132)",
                    "rgb(75, 192, 192)",
                    "rgb(153, 102, 255)",
                ],
                "borderWidth": 1
            }
        ]
    },
    options: {
        responsive: true,
        maintainAspectRatio: false,
        scales: {
            yAxes : [{
                "ticks" : {
                    "beginAtZero" : true
                }
            }]
        }
    }
});
new Chart( "build-results", {
    type: 'bar',
    data: {
        labels: ["Simple", "Cached", "Cached + Docker"],
        datasets : [
            {
                label: "Build job elapsed time (s)",
                data: [86,57,62],
                "fill": false,
                "backgroundColor": [
                    "rgba(255, 99, 132, 0.2)",
                    "rgba(75, 192, 192, 0.2)",
                    "rgba(153, 102, 255, 0.2)",
                ],
                "borderColor": [
                    "rgb(255, 99, 132)",
                    "rgb(75, 192, 192)",
                    "rgb(153, 102, 255)",
                ],
                "borderWidth": 1                
            }
        ]
    },
    options: {
        responsive: true,
        maintainAspectRatio: false,
        scales: {
            yAxes : [{
                "ticks" : {
                    "beginAtZero" : true
                }
            }]
        }
    }
});
new Chart( "deploy-results", {
    type: 'bar',
    data: {
        labels: ["Simple", "Cached", "Cached + Docker"],
        datasets : [
            {
                label: "Deploy job elapsed time (s)",
                data: [57,56,35],
                "fill": false,
                "backgroundColor": [
                    "rgba(255, 99, 132, 0.2)",
                    "rgba(75, 192, 192, 0.2)",
                    "rgba(153, 102, 255, 0.2)",
                ],
                "borderColor": [
                    "rgb(255, 99, 132)",
                    "rgb(75, 192, 192)",
                    "rgb(153, 102, 255)",
                ],
                "borderWidth": 1                
            }
        ]
    },
    options: {
        responsive: true,
        maintainAspectRatio: false,
        scales: {
            yAxes : [{
                "ticks" : {
                    "beginAtZero" : true
                }
            }]
        }
    }
});

</script>

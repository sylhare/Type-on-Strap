---
layout: post
title: Mocking external API with wiremock
description: Mocking external API with wiremock. Prepare a docker-compose bootstrap project for wiremock. https://github.com/adriangalera/docker-compose-wiremock
author-id: "galera"
categories: [testing]
tags: [java,testing,e2e,docker,wiremock]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/wiremock/featured.jpg"
thumbnail: "assets/img/posts/wiremock/featured.jpg"
image: "assets/img/posts/wiremock/featured.jpg"
redirect_from:
  - /2018/11/27/mocking-external-apis-with-wiremock/
---
<p class="p1">What happens when you are developing a component that heavily rely on an external API you do not control? Or even worst, that still does not exist. How could you test your component without connecting the external dependency? When we don't have control over the API that we need to integrate, we need a tool like a "mock server".</p>
<p class="p1">This article will discover and provide a bootstrap project for wiremock. More info: <a href="http://wiremock.org/">wiremock.org</a></p>
<p><!--more--></p>
<p>Quoting from their website:</p>
<blockquote><p>WireMock is a simulator for HTTP-based APIs. Some might consider it a <strong>service virtualization</strong> tool or a <strong>mock server</strong>.</p></blockquote>
<p>At its core is a Java software that receives HTTP requests with some mapped requests to responses</p>
<p>TL;DR: <a href="https://github.com/adriangalera/docker-compose-wiremock/" target="_blank" rel="noopener">https://github.com/adriangalera/docker-compose-wiremock/</a></p>
<h3>Configuring wiremock</h3>
<p>Configuring wiremock only consists on defining the requests to be mocked and the response that should be answered on the presence of the mocked request.</p>
<h3>Docker</h3>
<p>One nice way of integrate wiremock with your current testing environment is using it inside docker. There's this project <a href="https://github.com/rodolpheche/wiremock-docker">https://github.com/rodolpheche/wiremock-docker</a> that provides the wiremock service to docker.</p>
<p>In order to configure it, you must create the following folder structure:</p>
```
.
├── Dockerfile
└── stubs
    ├── __files
    │   └── response.json
    └── mappings
        └── request.json
```

<p>The mappings folder contains all the mocked requests definitions and __files contains the response JSON for the mocked requests as shown before.</p>
<h3>Example</h3>
<p>Let's say we have an external API developed by another team in the company under the host externalapi.com and is not yet finished. The call that our service needs to perform is externalapi.com/v1/resource/resource1 and will respond:</p>

```
{
    "hello":"world"
}
```
<p>Let's configure wiremock, so we can start working on our service in parallel with the other team.</p>

- Configure the request mapping

```
{
    "request":{
        "method":"GET",
        "urlPathPattern":"/v1/resource/([a-zA-Z0-9-\\_]*)"
    },
    "response":{
        "status":200,
        "bodyFileName":"response.json",
        "headers":{
            "Content-Type":"application/json"
        }
    }
}
```
- Configure the response

```
{
    "hello":"world"
}
```

- Test it

```
docker-compose up --build -d
curl http://localhost:7070/v1/resource/resource1
{
"hello" : "world"
}
```

<p>Yay! It worked!</p>

<p>The only missing point is configure the actual component to point to the mocked server. For example with <a href="https://github.com/Netflix/ribbon">ribbon</a>:</p>

```
externalservice.ribbon.listOfServers=http://localhost:7070
```
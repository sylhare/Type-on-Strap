---
layout: post
title: Unix useful commands
description: Compilation of unix (Linux/Mac OS X) useful commands, for instance, getting env variables by variables
author-id: "galera"
categories: [unix,scripting]
tags: [unix,scripting,linux,macosx]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/unix/featured-image.jpg"
thumbnail: "assets/img/posts/unix/featured-image.jpg"
image: "assets/img/posts/unix/featured-image.jpg"
---
<p>This article is a compilation of Unix useful commands to solve multiple issues found in my career.</p>
<p><!--more--></p>
<h2>Environment variables by variable</h2>

Image we need to access the contents of an environment variable, but its name is stored in another variable.

For example:
```bash
DEV_AWS_ACCESS_KEY_ID="1234-dev"
RD_AWS_ACCESS_KEY_ID="abcd-prd"
```
This could happen for instance while configuring multiple AWS in a CI system. 

Let's continue with the example, the CI system provide a variable called "stage", which can be `dev` or `prd`; then we want to prepend the content of this variable to get the credentials to the proper account:

```bash
CUR_ENV=`echo ${stage} | tr a-z A-Z`
ENV_ACCESS_KEY="${CUR_ENV}_AWS_ACCESS_KEY_ID"
```

Now the magic comes, if this the shell is based in `bash`, we could the technique called as "variable indirection", like this:

```bash
ACTUAL_KEY=echo ${!ENV_ACCESS_KEY}
```

However, this will not work on all the shells, a more general solution could be:

```bash
eval ACTUAL_KEY=\$$ENV_ACCESS_KEY
```

However, there might be security implications by using `eval`

<h2>Check if there are git changes in script</h2>

In a CI pipeline, you might want to check if there are changes to create an automatic commit, etc.

You can do that by running the following snippet:

```
git diff-index --quiet HEAD
ANY_CHANGE=$?
[ $ANY_CHANGE -ne 0 ] && echo "Do something with the change"
```
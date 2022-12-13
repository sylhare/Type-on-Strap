---
layout: post
title: Password protection for Github Pages
description: How to protect a static page served from Github pages with password and avoid link sharing without using any backend
author-id: "galera"
categories: [javascript, browser, github-pages, security]
tags: [javascript, browser, github-pages, security]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/github-password/featured.jpg"
thumbnail: "assets/img/posts/github-password/featured.jpg"
image: "assets/img/posts/github-password/featured.jpg"
---
In this article I implement a workaround to protect with password a static page stored in Github pages.

<p><!--more--></p>

I am using Github Pages to store static pages without any backend. That's super nice, but now I need to serve a password protected page. What kind I do?

I have found a guy that asked the same question and have a nice proposal: <b>use hashes</b>. You can find his code here: <a href="https://github.com/chrissy-dev/protected-github-pages">https://github.com/chrissy-dev/protected-github-pages</a>. His solution is older than walking.

## Hash to the rescue

The workaround is really simple, you choose a password that can be hard to guess. Hint: use some page that computes the password strength like: <a href="https://www.idstrong.com/tools/password-strength-checker/">https://www.idstrong.com/tools/password-strength-checker/</a>.

With that word, you generate the sha1 hash:

```
echo -n "<your-word>" | openssl sha1
cb1dc474e185777dad218b7d60f2781723d8190b
```

Now generate a folder with that text in the root of the repo and place all the password protected content there.

Then in the root of the repository place an index.html that will have a password form. 

When the user enters the password, compute the sha1 hash of the field they just enter. Then, perform a redirection to the URL, if the answer is different than 200, the folder has not been found, so the password is invalid.

That's the code that does the magic:

```javascript
function login(secret) {
            var hash = sha1(secret)
            var url = hash + "/index.html"
            var alert = document.querySelectorAll('[data-id="alert"]')

            var request = new XMLHttpRequest()
            request.open('GET', url, true)

            request.onload = function () {
                if (request.status >= 200 && request.status < 400) {
                    window.location = url
                } else {
                    parent.location.hash = hash
                    alert[0].style.display = 'block'
                    password[0].setAttribute('placeholder', 'Incorrect password')
                    password[0].value = ''
                }
            }
            request.onerror = function () {
                parent.location.hash = hash
                alert[0].style.display = 'block'
                password[0].setAttribute('placeholder', 'Incorrect password')
                password[0].value = ''
            }
            request.send()
        }

button[0].addEventListener("click", function () {
    login(password[0].value)
})
```
That works really nice, however, once the authentication is passed, a user can share the link and the authentication will be bypassed. We need an extra layer of security

## Avoid sharing the link

You can make a very easy implementation to make short-lived links. The page that provides the autentication can add a `not-valid-after` parameter to the URL. For example:

```javascript
request.onload = function () {
        if (request.status >= 200 && request.status < 400) {
            let nva = new Date().getTime() + 1_000
            window.location = url + "?nva="+nva
...
```
With this piece of code, the links are only valid for 1 second.

In the index.html of the protected content, we only need to check for the `not-valid-after` parameter:

```javascript
let paramString = window.location.search.split('?')[1];
let queryString = new URLSearchParams(paramString);
let nva = parseInt(queryString.get("nva"))
let now = new Date().getTime()
if (Number.isNaN(nva) || now > nva) {
    console.log("not-valid-after invalid, going to redirect to /")
    window.location = "/"
}
```
This piece of the code reads the `nva` parameter from the URL and checks if it's not present or if it's too late. In any of those cases, it redirects to the root; where the password form appears.
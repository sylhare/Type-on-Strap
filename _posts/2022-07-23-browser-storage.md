---
layout: post
title: Browser storage
description: How to store data in browser in order to make standalone apps without backend. This articles explores using browser localstorage and caches.
author-id: "galera"
categories: [javascript, browser, persistence]
tags: [javascript, browser, persistence]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/browser-storage/featured-image.jpg"
thumbnail: "assets/img/posts/browser-storage/featured-image.jpg"
image: "assets/img/posts/browser-storage/featured-image.jpg"
---
In this article I discuss some techniques to store data in the browser. This way the web applications do not require a expensive backend, every customer stores its information in the browser.

<p><!--more--></p>

This is part of my series of articles about leaflet:

- <a href="/leaflet-fog-of-war">Leaflet fog of war</a>
- <a href="/leaflet-draw-polygon-markers">Draw a polygon from markers in leaflet</a>
- <a href="/leaflet-load-gpx">Load and display GPX in leaflet</a>
- <a href="/browser-storage">Browser storage</a>

For the implementation of fog of war map, I do not want to spend any time dealing with the backend. Besides that, I don't want to spent not even a cent on the storage of data.

In this scenario, the visited areas are stored into a huge GeoJSON document. The persistence should store that document so the user does not need to re-create it every time.

However, how the can the data be persisted without any backend? It turns out the browser offers some persistence capabilities, let's analyse them.

## Local storage

The first approach is to use browser's local storage. This storage is a key-value storage which has a very simple synchronous contract:

```javascript
window.localStorage.setItem("a","b")
window.localStorage
> Storage {a: 'b', length: 1}
window.localStorage.getItem("a")
> 'b'
```

This storage is really simple and easy to use, however it comes at the cost of have a very limited space in the order of few MBs.

The original approach was to use this type of storage, however I reached the size limit really fast when I started to import GPX files.

## Caches

Reading a little bit more on browser storage capabilities, I discovered the caching mechanism. This is designed to store the answers from HTTP calls, hence its name. However, its original purpose can be violated to store any kind of data, not only HTTP responses.

This new API is asynchronous, that make the transition from local storage to caches a little bit painful, but it's a price we have to pay for having a massive amount of storage capability. According to <a href="https://web.dev/cache-api-quick-guide/">this article</a> the storage availability is based on the amount of storage available on the disk.

```javascript
function GeoJsonStorage() {
    const CACHE_NAME = "geojson"
    const CACHE_KEY = "https://xxxx/geojson.json"
    return {
        set: function (geojson) {
            caches.open(CACHE_NAME)
                .then(function (cache) {
                    cache.put(CACHE_KEY, new Response(JSON.stringify(geojson)));
                })
                .catch(err => console.log(`Cannot open the cache, error: ${err}`))
        },

        get: async function () {
            return caches.open(CACHE_NAME)
                .then(cache => cache.match(CACHE_KEY))
                .then(response => {
                    if (response)
                        return response.json()
                    return undefined
                })
                .catch(err => console.log(`Cannot get the contents from the cache, error: ${err}`))
        },
        clear: function () {
            caches.delete(CACHE_NAME)
        }
    }
}
```
The key to store arbitrary data into the caches mechanism is to trick the system saying that the cache key is an HTTP request: `https://xxxx/geojson.json`. This way, you can put and retrieve a JSON inside the caching mechanism

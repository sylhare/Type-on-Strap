---
layout: post
title: Load GPX in leaflet
description: Following up from the previous article about implementing the fog of war in leaflet, I want to be able to load a GPX file and display it in the map
author-id: "galera"
categories: [leaflet, javascript, browser, gis, gpx]
tags: [leaflet, javascript, browser, gis, gpx]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/leaflet-load-gpx/featured-image.jpg"
thumbnail: "assets/img/posts/leaflet-load-gpx/featured-image.jpg"
image: "assets/img/posts/leaflet-load-gpx/featured-image.jpg"
---
<style type="text/css">
.image-table td{
    border: 0px;
}
.image-table .center{
    text-align: center;
}
</style>
Following up from the previous article about implementing the fog of war in leaflet, I want to be able to load a GPX file and display it in the map.

<p><!--more--></p>

This is part of my series of articles about leaflet:

- <a href="/leaflet-fog-of-war">Leaflet fog of war</a>
- <a href="/leaflet-draw-polygon-markers">Draw a polygon from markers in leaflet</a>
- <a href="/leaflet-load-gpx">Load and display GPX in leaflet</a>
- <a href="/browser-storage">Browser storage</a>

Now, that I can draw paths with some distance in the map, I want to be able to load <a href="https://en.wikipedia.org/wiki/GPS_Exchange_Format">GPX files</a>. Those files store geographical information in the way of latitude,longitude and optionally elevation, time and other information. Those documents are XML based, so it should be easy to parse.

## Load file

How to load any file from the computer from the browser? `window.showOpenFilePicker` allows the browser to opens a file picker. That file picker can be customised to load files of only one type, to allow multiple files, etc. Take into account that is an async operation, so you'll have to deal with promises:

```javascript
const gpxPickerOpts = {
    types: [
        {
            description: 'GPX Files',
            accept: {
                'application/gpx+xml': ['.gpx']
            }
        }
    ],
    multiple: true
};
const fileHandlers = await window.showOpenFilePicker(gpxPickerOpts);
for (let fh of fileHandlers) {
    const file = await fh.getFile();
    const content = await file.text();
}
```
In this case, we want to load multiple GPX files. The `window.showOpenFilePicker` methods returns an array of file handlers which later we need to open and consume. `file.text()` operation returns the full contents of the file in text.

## Parse the GPX

Now that we have the text contents of the file, we must extract the relevant information, i.e. latitude and longitude of the points stored. Since GPX is XML-based, we can use two approaches to parse the file:

- SAX: Simple API for XML. Event based, when a new node is detected, an event is generated and passed to the event handler. Extremely efficient but complicated to implement.
- DOM: Document Object Mapper. It parses the file in one go. If the documents are big can lead to performance decrease, but it's extremely easy to implement.

Since this is a toy project, let's use DOM because of its simplicity:

```javascript
const _xmlTrackPointToLatLng = (trkpoint) => {
    return [parseFloat(trkpoint.attributes.lat.nodeValue), parseFloat(trkpoint.attributes.lon.nodeValue)]
}
var gpxDom = (new DOMParser()).parseFromString(content, 'text/xml');
const trackPoints = Array.from(gpxDoc.getElementsByTagName("trkpt"));
const latlngs = trackPoints.map((trkpnt) => container._xmlTrackPointToLatLng(trkpnt))
```

This snippet generated the DOM from the contents of the GPX file. Later extract the elements in the DOM tree that belongs to points in the track: `<trkpt>`. Then, every tag is processed to extract latitude and longitude.

In this point, we have an array of points with latitude/longitude pairs.

## Display the GPX

GPX can have a massive amounts of points, I have some of them with 20k points. In order to display it smoothly, I'm grouping the points in groups of 200 to draw the polygon like in the previous article (<a href="/leaflet-draw-polygon-markers">Draw a polygon from markers in leaflet</a>):

```javascript
const  _group = (arr, n) => {
    const res = [];
    let limit = 0;
    while (limit + n <= arr.length) {
        res.push(arr.slice(limit, n + limit));
        limit += n
    }
    return res
}
const groups = container._group(latlngs, 200)
var polygonGeoJSON = undefined
for (let group of groups) {
    const polLatLng = _joinLinesInPolygon(group) //from previous article
    const pol = L.polygon(polLatLng).toGeoJSON()
    if (!polygonGeoJSON) {
        polygonGeoJSON = pol
    } else {
        polygonGeoJSON = turf.union(pol, polygonGeoJSON)
    }
}
```

This `_group` function creates batches of 200 points in order to perform the joinin of those points into a polygon. Once a polygon is created, I'm merging them with the <a href="https://turfjs.org/docs/#union">turf library</a> by performing a union.

This way of displaying GPX files produces very appealing representation such as:

<img src="/assets/img/posts/leaflet-load-gpx/gpx.png" alt="GPX representation in leaflet map"/>

Here you can see it in action: <a href="https://www.agalera.eu/leaflet-fogofwar/" target="_blank" rel="noopener">https://www.agalera.eu/leaflet-fogofwar/</a>
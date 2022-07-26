---
layout: post
title: Draw a polygon from markers in leaflet
description: How to join multiple markers in leaflet into a polygon
author-id: "galera"
categories: [leaflet, javascript, browser, gis]
tags: [leaflet, javascript, browser, gis]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/draw-polygon-markers/featured-image.jpg"
thumbnail: "assets/img/posts/draw-polygon-markers/featured-image.jpg"
image: "assets/img/posts/draw-polygon-markers/featured-image.jpg"
---
<style type="text/css">
.image-table td{
    border: 0px;
}
.image-table .center{
    text-align: center;
}
</style>
Following up from the previous article about implementing the fog of war in leaflet, I want to add some markers to the map and to create a polygon that joins them.

Let's see how I manage to do that.

<p><!--more--></p>

The user can click the map to generate markers that follow a route, e.g. a road, a trail, etc.. and later join those markers to create a complex polygon. This reflects the fact that you have visited the road, but you only have certain visibility of the environment (maybe 10 meters or so):

<table class="image-table">
<tr>
<td>
<img src="/assets/img/posts/draw-polygon-markers/1.png" alt="Markers"/>
</td>
<td>
<img src="/assets/img/posts/draw-polygon-markers/2.png" alt="Joined markers into a polygon"/>
</td>
</tr>
<tr>
<td class="center">
<small>Markers</small>
</td>
<td class="center">
<small>Joined markers into a polygon</small>
</td>
</tr>
</table>

## Join markers into a polygon

In order to compute all the geographic information, I'll use the <a href="https://github.com/bjornharrtell/jsts">jsts</a> library.

When a new marker is added, it added to the map and to a internal array:
```javascript
        onAdd: function (e) {
            const marker = new L.Marker(e.latlng)
            container.markers.push(marker.addTo(map))
        }
```
When the user click on a button, those markers are processed and joined into a polygon:

```javascript
const _joinLinesInPolygon = (points) => {
    const pointToGeomCoordinate = (p) => {
        if (p.lat && p.lng)
            return new jsts.geom.Coordinate(p.lat, p.lng)
        return new jsts.geom.Coordinate(p[0], p[1])
    }

    const toLeafletPoint = (p) => {
        return [p.x, p.y]
    }

    const meters = 40 //the user can selected the width of the generated polygon
    const distance = (meters * 0.0001) / 111.12; //Geometry aproximations
    const geometryFactory = new jsts.geom.GeometryFactory();
    const pathCoords = points.map((p) => pointToGeomCoordinate(p));
    const shell = geometryFactory.createLineString(pathCoords);
    const polygon = shell.buffer(distance);
    const polygonCoords = polygon.getCoordinates();
    return polygonCoords.map((coord) => toLeafletPoint(coord))
}
```
This method converts the leaflet points into a format the jsts libray can understand and perform the `buffer` operation in jsts which does all the magic. Later the coordinates are transformed into the leaflet format
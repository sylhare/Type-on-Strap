---
layout: post
title: Latitude,longitude bounds
description: This post describe the basic algorithm to calculate the bounds of a set of points coordinates with latitude and longitude.
author-id: "galera"
categories: [algorithms, typescript, latlng, maps]
tags: [algorithms, typescript, latlng, maps]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/bounds-lat-lng/featured-image.jpg"
thumbnail: "assets/img/posts/bounds-lat-lng/featured-image.jpg"
image: "assets/img/posts/bounds-lat-lng/featured-image.jpg"
---
I'm currently developing an application based on maps. In that application I want to represent a set of markers. In order to do so, the map library I'm using it has a `fitBounds` method. However, you need to compute the bounds of the map that allow all the markers to be visible. I describe in this article the implemented algorithm.

<p><!--more--></p>

## Abstraction of latitude longitude

First of all, we need to do a nasty approximation. We can represent the earth globe, as a 2-D cartesian axis. In order to do that, we can consider latitude as the y axis and longitude as the x axis.

Longitude will be in range [-180,180] and latitude in the range [-90,90]. We can define the cardinal points in the chart:

- North -> Point in (0,90)
- East -> Point in (180,0)
- South -> Point in (0,-90)
- West -> Point in (-180,0)

We can define the bounds of a set of points using two points only: North-East point and South

- North-East -> Point in (180,90)
- South-West -> Point in (-180,-90)

We can see a visual representation of those points in the following chart:

<div style="height: 400px">
<canvas id="myChart"></canvas>
</div>


Taking this into account we can set up some algorithm to calculate the bounds.

## Algorithm to compute the bounds

The first basic algorithm is to iterate over all points and compare the longitude (x) and latitude (y) and obtain the point with higher x and y and the point with lower x and y.

```typescript
const SW: LatLngTuple = [-90, -180]
const NE: LatLngTuple = [90, 180]
export const ALL_WORLD_BOUNDS: LatLngBoundsExpression = [NE, SW]
export const getBoundsFromPoints = (points: Point[]): LatLngBoundsExpression => {

    if (points.length === 0)
        return ALL_WORLD_BOUNDS

    let nex = 0, swx = 0, ney = 0, swy = 0
    points.forEach((point) => {
        if (nex === 0 && swx === 0 && ney === 0 && swy === 0) {
            nex = swx = point.longitude
            ney = swy = point.latitude
        } else {
            if (point.longitude > nex) nex = point.longitude;
            if (point.longitude < swx) swx = point.longitude;
            if (point.latitude > ney) ney = point.latitude;
            if (point.latitude < swy) swy = point.latitude;
        }
    })
    return [[ney, nex], [swy, swx]]
}
```

Bear in mind that the map library expects and array of `[lat,lng]`, that's why we are switching the natural order of x,y and we're using y,x that corresponds to `[lat,lng]`.

## Unit test

These are the unit tests that check this algorithm behaves correctly:

```typescript
import {ALL_WORLD_BOUNDS, getBoundsFromPoints} from "./BoundCalculator";
import {Point} from "../../types/Point";

test("should compute bounds of an empty list", () => {
    const bounds = getBoundsFromPoints([])
    expect(bounds).toEqual(ALL_WORLD_BOUNDS)
})

test("should compute bounds of one point", () => {
    const lat = 1, lng = 1
    const point = new Point(lat, lng)
    const bounds = getBoundsFromPoints([point])
    expect(bounds).toEqual([[lat, lng], [lat, lng]])
})

test("should compute bounds a list of points, each point per quadrant", () => {
    const lat1 = 30, lng1 = 90
    const lat2 = 30, lng2 = -90
    const lat3 = -30, lng3 = -90
    const lat4 = -30, lng4 = 90
    const point1 = new Point(lat1, lng1)
    const point2 = new Point(lat2, lng2)
    const point3 = new Point(lat3, lng3)
    const point4 = new Point(lat4, lng4)
    const bounds = getBoundsFromPoints([point1, point2, point3, point4])
    expect(bounds).toEqual([[lat1, lng1], [lat3, lng3]])
})
```

<script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/2.9.3/Chart.bundle.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/chartjs-plugin-datalabels"></script>
<script>
var ctx = 'myChart';

var north = {x: 0, y: 90}
var east = {x: 180, y: 0}
var south = {x: 0, y: -90}
var west = {x: -180, y: 0}
var northEast = {x: 180, y: 90}
var southWest = {x: -180, y: -90}
const labels = ["N","E","S","W","NE","SW"]
Chart.helpers.merge(Chart.defaults.global, {
			plugins: {
				legend: false,
				title: false
			}
		});

var myLineChart = new Chart(ctx, {
    type: "scatter",
    data: {
        labels: labels,
        datasets : [
            {
                label: "Cardinal points",
                data: [north, east, south,west, northEast, southWest]
            }
        ]
    },
    options: {
        responsive: true,
        layout: {
                padding: {
                    left: 0,
                    right: 0,
                    top: 50,
                    bottom: 50
                }
            },       
        maintainAspectRatio: false,
        plugins: {
            datalabels: {
                align: 'end',
                anchor: 'end',
                color: function(context) {
                    return context.dataset.backgroundColor;
                },
                font: function(context) {
                    var w = context.chart.width;
                    return {
                        size: w < 512 ? 12 : 14
                    };
                },
                formatter: function(value, context) {
                    return context.chart.data.labels[context.dataIndex];
                }
            }
        },
    }    
})

</script>

---
layout: post
title: Download excel file using AngularJS
tags: [javascript, AngularJS]
author-id: oppalove
excerpt_separator: <!--more-->
---
# Overview
How to download excel file using Angular-JS
<!--more-->


### Create button and add ng-click
```html
<button type="button" class="btn btn-success" ng-click="export()"> 
<span class="glyphicon glyphicon-export"/> Export </button>
```


### Define function
```javascript
$scope.export = function () {
    $http.post('/table/export.do', req,{responseType: 'arraybuffer'}).then(function (response) {
        var header = response.headers('Content-Disposition');
        var fileName = header.split("=")[1].replace(/"/gi,'');
        console.log(fileName);
        
        var blob = new Blob([response.data],
            {type : 'application/vnd.openxmlformats-officedocument.presentationml.presentation;charset=UTF-8'});
        var objectUrl = (window.URL || window.webkitURL).createObjectURL(blob);
        var link = angular.element('<a/>');
        link.attr({
            href : objectUrl,
            download : fileName            
        })[0].click();
    })
};
```
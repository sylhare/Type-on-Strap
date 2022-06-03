---
layout: post
title: Publishing versions in Gradle BOM
description: With the current behaviour of the gradle BOM generation the version of the dependencies are not available from the client project. In this article we make it available.
author-id: "galera"
categories: [java]
tags: [gradle,groovy]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/gradle-bom/featured.png"
thumbnail: "assets/img/posts/gradle-bom/featured.png"
image: "assets/img/posts/gradle-bom/featured.png"
---
<p>We are planning to use the concept of Bill Of Materials (BOM) to define the version of our dependencies. However this externalization do not allow us to have the version variable in the projects that import the bom. Here you will find out how we managed to get the version in our components.</p>
<p><!--more--></p>
<h2>Bill of Materials</h2>
Everyone knows dependencies are hell, specially dealing with versions. That's why we want to centralize all the version definition in a Bill Of Materials file. This approach is used for instance by <a href="https://spring.io/blog/2015/02/23/better-dependency-management-for-gradle">Spring framework</a>.

However in our projects we need to provide the current version of libraries for some external component configuration.

In this sense, the two requirements collide because the versions cannot be extracted from the bom file.

<h2>Publishing the versions</h2>

Here's our bom gradle script:

```groovy
plugins {
    id 'java'
    id 'maven-publish'
    id 'io.spring.dependency-management' version '1.0.6.RELEASE'
}

group 'com.example'

def versions = [
        library1: '5.14.13',
        library2: '1.0.6.RELEASE',  
]


dependencyManagement {
    dependencies {
        dependency group: 'com.example', name: 'library1', version: versions.library1
        dependency group: 'com.example', name: 'library2', version: versions.library2    
    }
}

publishing {
    repositories {
        maven {
            url 'https://repo.com/maven-releases'
            credentials credentials
        }
    }

    publications {
        mavenBom(MavenPublication) {
            artifacts = []
        }
    }
}
```

When we execute the publish it correctly generates the bom file. However the versions are not accesible from the projects that import the bom.

In order to extract the version, we can add it as properties to the bom file:

```groovy
 publications {
        mavenBom(MavenPublication) {
            artifacts = []
            pom.withXml {
                def propertiesNode = new Node(null, "properties")
                versions.each { entry -> propertiesNode.appendNode(entry.key, entry.value) }
                asNode().append(propertiesNode)
            }
        }
    }
```

It took like 3 hours to get these 3 lines working because Groovy sucks ... After that we get the versions in the pom.xml:

```xml
<?xml version="1.0" encoding="UTF-8"?>
<project xmlns="http://maven.apache.org/POM/4.0.0" xsi:schemaLocation="http://maven.apache.org/POM/4.0.0 http://maven.apache.org/xsd/maven-4.0.0.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <modelVersion>4.0.0</modelVersion>
  ...
  <properties>
    <library1>5.14.13</library1>
    <library2>1.0.6.RELEASE</library2>
  </properties>
</project>

```

From the project that want to use the bom we can access these versions like this:

```groovy
dependencyManagement.importedProperties['library1']
```

---
layout: post
title: "Relational model in DynamoDB"
description: In this post I implement the persistence of a relational model in DynamoDB. Fix some challenges using ECS to perform long running DynamoDB operations
author-id: "galera"
categories: [aws]
tags: [aws,cloud,dynamodb,relational,persistence,ecs,python]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/relationalmodel/featured.jpg"
thumbnail: "assets/img/posts/relationalmodel/featured.jpg"
image: "assets/img/posts/relationalmodel/featured.jpg"
redirect_from:
  - /2019/01/28/relational-model-in-dynamodb/
---
<p>For a new project that I am starting, I have the requirement to implement a highly relational model in AWS. The most important condition is to keep costs as low as possible, since this is a personal project and I do not want to get poor.</p>
<p>Therefore I will persist the model in DynamoDB configured to use the minimum resources as possible.</p>
<p><!--more--></p>
<p>The application consist on three entities: User,Map and Points.</p>
<p>Users can create multiple maps that contain several points. The following UML schema explain the relationships:</p>

[ ![Relational model UML](/assets/img/posts/relationalmodel/estuve-model.png) ](/assets/img/posts/relationalmodel/estuve-model.png)
*Relational model UML*

<p>DynamoDB is a key-value store with support for range key. Thanks to that I am able to implement the following queries:</p>

- CRUD User,Map,Point
- Add a map for one user
- Add a point in a map
- Get points from a map
- Remove map from user
- Remove point from user
- Get maps for a user

<h2>The DynamoDB model</h2>
<h5>Users table</h5>
<p>The user table is straightforward, the only key is a unique identifier for the user.</p>

```python
{
    "TableName": USER_TABLE,
    "KeySchema": [
        {
            "AttributeName": "user_id",
            "KeyType": "HASH"
        }
    ],
    "AttributeDefinitions": [
        {
            "AttributeName": "user_id",
            "AttributeType": "S"
        }
    ]
}
```


<p>There are additional attributes that keep track of the number of points and maps stored for that user:</p>

```python
record = {
    'user_id': {
        'S': str(obj.user_id)
    },
    'num_points': {
        'N': str(obj.num_points)
    },
    'num_maps': {
        'N': str(obj.num_maps)
    }
}
```

<h5>Maps table</h5>
<p>The map table is a little bit more complex, because it has to keep relations between users and maps. Therefore, I use the range key to save the unique identifier of the map:</p>

```python
{
    "TableName": MAPS_TABLE,
    "KeySchema": [
        {
            "AttributeName": "user_id",
            "KeyType": "HASH"
        },
{
            "AttributeName": "map_id",
            "KeyType": "RANGE"
        },
    ],
    "AttributeDefinitions": [
        {
            "AttributeName": "map_id",
            "AttributeType": "S"
        },
{
            "AttributeName": "user_id",
            "AttributeType": "S"
        }
    ]
}
```

<p>There are additional attributes associated to the map (self-explanatory):</p>

```python
{
    'user_id': {
        'S': str(obj.user_id)
    },
    'map_id': {
        'S': str(obj.map_id)
    },
    'name': {
        'S': str(obj.name)
    },
    'description': {
        'S': str(obj.description)
    },
    'num_points': {
        'N': str(obj.num_points)
    }
}
```

<h5>Points table</h5>
<p>This is most complex table. The keys are similar to the maps, the range key is used to store the unique identifier of the map:</p>

```python
{
    "TableName": POINTS_TABLE,
    "KeySchema": [
        {
            "AttributeName": "map_id",
            "KeyType": "HASH"
        },
        {
            "AttributeName": "point_id",
            "KeyType": "RANGE"
        },
    ],
    "AttributeDefinitions": [
        {
            "AttributeName": "map_id",
            "AttributeType": "S"
        },
        {
            "AttributeName": "point_id",
            "AttributeType": "S"
        },
    ]
}
```
<p>And the additional parameters:</p>
```python
{
    'point_id': {
        'S': obj.point_id
    },
    'map_id': {
        'S': str(obj.map_id)
    },
    'lat': {
        'S': str(obj.lat)
    },
    'lon': {
        'S': str(obj.lon)
    },
    'date': {
        'N': str(obj.epoch)
    },
    'name': {
        'S': str(obj.name)
    },
}
```
<p>The challenge with this model is to be able to delete a map with a large number of points. It is counter-intuitive, because one might think that removing only the points with the primary key of the map will make the work but ...</p>
<blockquote><p><strong>THIS WILL NOT WORK!</strong></p></blockquote>
<p>Due to the way DynamoDB is implemented, this is not possible (<a href="https://stackoverflow.com/questions/34259358/dynamodb-delete-all-items-having-same-hash-key">https://stackoverflow.com/questions/34259358/dynamodb-delete-all-items-having-same-hash-key</a>). In that kind of tables, you need to provide the primary key and the range key in order to delete an item.</p>
<p>Since the number of items can be large, it could take a lot of capacity to delete a the points. I do not want to consume that capacity, so I will let DynamoDB throttle the deletes to adapt to the capacity.</p>
<p>The project is serverless (Lambda) based and trying to delete a large number of points will result in timeouts when DynamoDB throttle the deletes. There are two possible solutions here: increase the write capacity of the table (increase cost) or increase the Lambda timeout (increase cost).</p>
<p>After thinking a little bit, the valid solution I choose is to launch an ECS Task with the logic to delete the large number of maps:</p>

```python
client.run_task(
    cluster="arn:aws:ecs:eu-west-1:***:cluster/remove-points-cluster",
    taskDefinition="arn:aws:ecs:eu-west-1:***:task-definition/remove-points",
    overrides={
        "containerOverrides": [
            {
                "name": "remove-points",
                "environment": [
                    {
                        "name": "MAP_ID",
                        "value": map_id
                    }
                ]
            }
        ]
    },
    launchType="FARGATE",
    networkConfiguration={
        "awsvpcConfiguration": {
            "subnets": ["1", "2"],
            "assignPublicIp": "ENABLED"
        }
    }
)
```
<p>The best part of this ECS Task is that only took 5 minutes to use the same code base and Dockerize the logic of removing the points!</p>
<p>Now the long-running task of delete a large number of points is done in ECS, where the pricing model is pay per use. Since this is a feature that is not going to happen a lot, it's perfectly fine.</p>

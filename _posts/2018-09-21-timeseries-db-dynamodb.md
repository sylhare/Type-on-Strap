---
layout: post
title: Timeseries database with DynamoDB
description: Design, challenges and implementation of a DynamoDB timeseries database using AWS Lambda, Kinesis and DynamoDB Streams
author-id: "galera"
categories: [aws]
tags: [aws,cloud,dynamodb,kinesis,python]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/dynamodbts/featured.jpg"
thumbnail: "assets/img/posts/dynamodbts/featured.jpg"
image: "assets/img/posts/dynamodbts/featured.jpg"
redirect_from:
  - /2018/09/21/creating-a-timeseries-database-with-dynamodb/
---
<p>Working with sensor data, we need to store timestamped data in a way that was easily checked and show it in graphs. This post explains the process we follow to end up implementing a timeseries database on DynamoDB.</p>
<p><!--more--></p>
<h2>Requirements</h2>
<p>The problem was pretty straightforward at the beginning: store temporal data that came from multiple sensors where one sensor could have more than one metric. For instance: one battery sensor might monitor the battery level as well as the battery voltage. The sensor data data flows continuously at an undetermined interval, it can be 1 minute, 5 minutes, ... Generally speaking the data can be represented as:</p>
<p><strong>METRIC NAME - TIMESTAMP - VALUE</strong></p>
<p>In order to support the defined usage scenarios we need to provide a tool able to:</p>

* Perform temporal queries: give me the metric between 2017/05/21 00:00:00 and 2017/08/21 00:30
* Support multiple granularity: second, minute, hour, ...
* Support TTL for cold data: expire cold data
* Perform in a cost-efficient manner

<h2>Available technologies</h2>
<p>When one thinks about database at first the relational ones appear as the first option. The problem with relational databases such as MySQL, PostgreSQL, .. is that we the database becomes very unusable when the size of the tables grows. And in this case the data  will grow a lot with the usage.</p>
<p>Furthermore, when the data is so big the indexes start to generate headaches making the queries take a lot of time.</p>
<p>Finding these drawbacks in the traditional relational databases, we shift towards NoSQL databases.</p>
<p>The first one that came into our mind (because we had some previous experience with it) was whisper <a href="https://github.com/graphite-project/whisper">https://github.com/graphite-project/whisper</a>. This database is a small component of the graphite project, basically is a wrapper to write the temporal data to a file performing multiple roll-up aggregations on the fly. This looked promising, however, when we heavy loaded it performed very poorly.</p>
<p>Since the platform we were building it was AWS based, we decided to analyse what Amazon can provide us and finally found <strong>DynamoDB</strong>!</p>
<h2>DynamoDB at the rescue!</h2>
<p>DynamoDB is a key-value database that supports the configuration of item TTL and its costs are predictable because are in function of the required capacity.</p>
<p>One might ask: <em>how can I store the presented model in a key-value database?</em> The magic comes with the composite key feature: <a href="https://aws.amazon.com/es/blogs/database/choosing-the-right-dynamodb-partition-key/">https://aws.amazon.com/es/blogs/database/choosing-the-right-dynamodb-partition-key/</a></p>



<h2>Knowing DynamoDB</h2>
<p>Quoting AWS documentation:</p>
<blockquote><p>This type of key is composed of two attributes. The first attribute is the <em>partition key</em>, and the second attribute is the<em> sort key</em>.</p></blockquote>
<p>Therefore, we can use the metric name as the partition key and the timestamp and the sort key and the sensor value as the DynamoDB object:</p>

[ ![DynamoDB Key Schema](/assets/img/posts/dynamodbts/dynamo-db-key.png) ](/assets/img/posts/dynamodbts/dynamo-db-key.png)
*DynamoDB Key Schema*



<h2>Rollup aggregations?</h2>
<p>Now we can store the timestamped data in DynamoDB, however, what happens with the multiple granularity requirement? DynamoDB has a nice feature called DynamoDB Streams.</p>
<p>DynamoDB Streams sends the events generated in the database (new records, edit, delete, ...) to a an AWS Lambda. Hence, we can perform the aggregations for the multiple granularities as soon as a new data value arrives. In order to perform the storage of the multiple aggregations, we can define one table for each aggregation.</p>
<h2>The implementation: DynamoDB Timeseries Database</h2>
<p>Finally, in order to complete the setup we have used a serverless approach in order to allocate the cost of the project to the required capacity.</p>
<p>The final structure of the implemented solution looked like this:</p>

[![Components of DynamoDB TimeseriesDB](/assets/img/posts/dynamodbts/dynamo-db-database-2.png)](/assets/img/posts/dynamodbts/dynamo-db-database-2.png)
*Components of DynamoDB TimeseriesDB*

<p>1) The service that uses the DynamoDB Timeseries database is a serverless application with an API Gateway calling an Lambda function. One of the steps of the business logic is communicate with the DynamoDB Insert Lambda.</p>
<p>2) The insertion mechanism of the database is by invoking an AWS Lambda function. This function inserts the timestamped data into the lower granularity table.</p>
<p>3) When the Insert function inserts the data into the lower granularity table, DynamoDB Streams invokes the AWS Lambda involved in performing the roll-up aggregations.</p>
<p>4) The Aggregate function has the logic implemented on how to perform multiple aggregations (average, sum, count, ...).</p>
<p>5) Each time serie can be configured to have different aggregations, TTL, timezone to perform temporal aggregation, etc... There are an additional lambda in order to configure the timeserie parameters.</p>
<p>6) Once the aggregation is performed, the data is stored into the appropriate DynamoDB table according to the granularity: minute, hour, day, month, year, ...</p>
<p>With this solution we achieve a solution where we can fine tune the capacity of the database.</p>
<h2>Production time!</h2>
<p>When we put this system into production, some issues arises with the capacity of the DynamoDB. When the incoming data flows faster than the reserved capacity in DynamoDB, the lambda functions became very slow and resulting in a high number of timeout errors.</p>
<p>The reasons for this is that when DynamoDB is running at capacity, it throttle requests, making the client adjust its speed to the reserved capacity. This works good for a traditional client, however it breaks with the serverless approach, because the Lambda functions were taking too much time.</p>
<p>In order to fix situation, we add another component to the system: an AWS Kinesis Stream <a href="https://aws.amazon.com/kinesis">https://aws.amazon.com/kinesis</a>. Instead of writing directly to the DynamoDB table, the service Lambda function now writes the data into a Kinesis stream.</p>
<p>In the other side of the stream, we place a Kinesis consumer that is able to consume items from the stream in batches of items. Additionally, we are able to control the insertion speed of items in DynamoDB by sleeping some time between batch consumption. Since this consumer needs to be 24/7, it runs on a traditional EC2 instance.</p>
<p>Now the scheme looks like this:</p>

[![Adding Kinesis into the TimeseriesDB](/assets/img/posts/dynamodbts/dynamo-db-kinesis-database.png)](/assets/img/posts/dynamodbts/dynamo-db-kinesis-database.png)
*Adding Kinesis into the TimeseriesDB*

<p>The additional point (7) shows the Kinesis stream where the items are being inserted by the Service Lambda function and consumed by the DynamoDB Timeseries DB Kinesis Consumer. The configured batch size (B) and sleep time (T) allows the consumer to buffer the insertion of data up to the reserved DynamoDB capacity.</p>
<h2>Show me the code</h2>
<p>An open sourced version of the production code can be found here:</p>
<p><a href="https://github.com/adriangalera/dynamodb-timeseries">https://github.com/adriangalera/dynamodb-timeseries</a></p>

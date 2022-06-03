---
layout: post
title: Discovering terraform
description: Recently I have been playing with Terraform tool and I wrote some basic use cases for the sake of sharing the knowledge
author-id: "galera"
categories: [devops,cloud]
tags: [devops,aws,terraform,cloud]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/terraform/terraformed.jpg"
thumbnail: "assets/img/posts/terraform/terraformed.jpg"
image: "assets/img/posts/terraform/terraformed.jpg"
---

Recently I have been playing with Terraform tool and I wrote some basic use cases for the sake of sharing the knowledge
<p><!--more--></p>

## What is terraform?

Imagine you start a new cloud project in AWS. At the very begining, you will write some scripts to generate your AWS resources:

- DynamoDB tables
- S3 buckets
- RDS databases
- ...

However, lots of configuration are required to make it work properly, i.e. IAM roles, ARN of resources, etc... So, managing this burden manually or with some scripts rapidly could become an issue.

Even worst, imagine you have resources in both Google Cloud an AWS. The amount of scripts will be doubled as well as the issues.

Terraform was created to help with those issues. It is a tool written in Go language that enables managing the state of cloud infrastructure from configuration files with its own Domain Specific Language (DSL). 

Furthermore, it is a handy tool to define the infrastructure as code. In this approach the changes to the infrastructure can be pushed to a git repository for instance; this way the changes can be reviewed, etc.

Another fancy feature is that it keeps the state of the resources running in Cloud. This way terraform allows to make incremental changes (if possible).

## Basic usage

There are three main commands to use:

- `terraform plan`: process the configuration templates and present in stdout the actions that will be applied
- `terraform apply`: process the configuration templates and apply the required actions
- `terraform destroy`: process the configuration templates and delete the resources

## Scenarios

In this section, I prepared 4 scenarions to show some cool features that terraform can provide.

### Basic

This scenario create a DynamoDB table and an S3 bucket. There's nothing fancy here, only basic usage of terraform:

```
variable "app-prefix" {
  type = string
  default = "comms-ks-01"
}
resource "aws_dynamodb_table" "configuration" {
  name           = "${var.app-prefix}_timeserie_configuration"
  billing_mode   = "PAY_PER_REQUEST"
  hash_key       = "timeserie1"

  attribute {
    name = "timeserie1"
    type = "S"
  }

}
resource "aws_s3_bucket" "s3-bucket-rnd-name" {
    bucket = "${var.app-prefix}-timeserie-configuration"
} 
output "bucket-arn" {
  value = "${aws_s3_bucket.s3-bucket-rnd-name.arn}"
}
```
The template defines a variable ("app-prefix), a DynamoDB table ("configuration") and a S3 bucket ("s3-bucket-rnd-name"). The last line defines to output to print the ARN of the created bucket

### For each

In this scenario, the requirement is to create three tables with similar names. In order to reuse the code, we can use terraform's `for-each` statement:

```
resource "aws_dynamodb_table" "configuration" {

  for_each = {
    test1: "${var.app-prefix}_configuration_test_1",
    test2: "${var.app-prefix}_configuration_test_2",
    test3: "${var.app-prefix}_configuration_test_3",
  }

  name           = each.value
  billing_mode   = "PAY_PER_REQUEST"
  hash_key       = "timeserie"

  attribute {
    name = "timeserie"
    type = "S"
  }
}
```
This way with only one resource definition, we can create multiple resources.

### Linking resources

Up to here, we haven't done anything fancy, just creating resources. However, we can create more interesting environments, such as a serverless pipeline to process a file:

[![S3-lambda pipeline](/assets/img/posts/terraform/s3-lambda.png)](/assets/img/posts/terraform/s3-lambda.png)

```
# Lambda definition
resource "aws_lambda_function" "download-s3-lambda" {
  filename      = "download-s3-file-lambda.zip"
  function_name = "${var.app-prefix}-download-files-lambda"
  role          = "${aws_iam_role.iam_for_lambda.arn}"
  handler       = "receive-file-s3.handler"
  runtime       = "python3.7"
  depends_on    = ["aws_iam_role_policy_attachment.lambda_logs", "aws_cloudwatch_log_group.example"]
}

# Notify lambda when a file is created in the S3 bucket
resource "aws_s3_bucket_notification" "bucket_notification" {
  bucket = "${aws_s3_bucket.s3-files-bucket.id}"

  lambda_function {
    lambda_function_arn = "${aws_lambda_function.download-s3-lambda.arn}"
    events              = ["s3:ObjectCreated:*"]
  }
}
# S3 bucket to place files
resource "aws_s3_bucket" "s3-files-bucket" {
    bucket = "${var.app-prefix}-files-bucket"
    force_destroy = "true"
}
```
The template in this environment is more complex, because there's a lot of IAM permissions in place. For simplicity, not all the resources are included in this article. For more info: <a href="https://github.com/adriangalera/terraform-knowledge-sharing">https://github.com/adriangalera/terraform-knowledge-sharing</a>.

In the upper code snippet can be observed how the AWS lambda and the AWS S3 bucket are created. Additionally, the notification from S3 to Lambda is created. When an object is created  in the S3 bucket, the AWS lambda will be invoked.

### Templates
Terraform also supports templates to enable code reuse. In this environment, we will create two docker container to be ran on ECS. Those two docker containers will execute the same image and echo the value of an environment variable.

If we did not have templates, we would need to define the resources twice. Thanks to the templates, we don't need to define the container definition twice, we only need to pass the variables.

```
data "template_file" "container_backend" {
  template = "${file("container_definition.tpl")}"
  vars = {
    container_name = "${var.app-prefix}_backend_container"
    log_group = "${var.app-prefix}_backend_container"
    service_type = "backend"
  }
}
  family = "${var.app-prefix}_backend_task_definition"
  requires_compatibilities = [ "FARGATE" ]
  network_mode =  "awsvpc"
  execution_role_arn = "${aws_iam_role.ecs_container_iam_role.arn}"
  cpu = 256
  memory = 512
  container_definitions = "${data.template_file.container_backend.rendered}"
}
```
In order to get the rendered values of the template, we need to get the `rendered` field.

In this case, the template contains lots of configuration and IAM definitions. For more info, check the whole repository: <a href="https://github.com/adriangalera/terraform-knowledge-sharing">https://github.com/adriangalera/terraform-knowledge-sharing</a>
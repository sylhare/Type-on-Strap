---
layout: post
title: AWS scheduled lambda functions using CloudWatch.
tags: [AWS, Lambda, CloudWatch, SNS, Java]
published: false
---

This tutorial will help you to gain some hands on experience with AWS Lambda, CloudWatch Alarms, SNS (Simple Notification Service) and Java. In particular, we look at scheduling AWS Lambda functions to run on a fixed schedule. This is useful if you have any task that want to run on a schedule, for instance you want to order yourself your favourite bottle of whisky on Thursday to make sure it arrives for Friday, and you want to do this each and every week. Better examples might be to run a job which sends batches of emails at a particular time, or triggers a backup of your website. The task itself (The Lambda function) is written in Java. AWS Lambda functions are pretty much just code which sits there and does nothing until we invoke it. The best part about that is, if we don't invoke it, it does nothing, and that includes costing nothing. With AWS Lambda you pay for only the time your lambda function is running. This allows you to have a task which runs on a schedule without paying for a server sat waiting around for a particular time.

Lambda functions can be written in a whole range of languages, but for this one I chose Java. Why? Because it was one of the first lambda functions I ever wrote, and so I picked a familiar language so I can focus on AWS and infrastructure. You can do this in Java, Node.js, Python, C#; take your pick.

The scenario we will go through is this: I'm a fan of a football team (Soccer for any Americans), Sheffield Wednesday, but due to my schedule I struggle to keep track of when the games on. The outcome of this project is that each and every day my Lambda function will check if Sheffield Wednesday are playing, and if they are, it will trigger SNS to send a text message to my phone at 10am on the day. This can be easily adapted to any team, any sport, or pretty much anything you want notifying of. Better still, the final solution is all in the AWS free tier, so it shouldn't cost you a penny*.
*Shouldn't cost you a penny - you should check on this yourself though, as the pricing structure of AWS does change from time to time. Even if it was outside of the free tier, the Lambda function runs for only a couple of seconds at most, and so we would be paying for 2 seconds of compute time per day. One second of compute time costs $0.00001667 per GB memory used. This will use less than 1GB of memory and will run for a couple of seconds a day, I think we are safe that the costs won't get out of hand. Text messages cost $0.00205 and will only be sent on the days Sheffield Wednesday play, so in a month that might add up to 1 cent, but under the free tier the first 100 a month are free anyway.
 
Architecture
In this tutorial, you will create an AWS CloudWatch Event to be scheduled every day. This event triggers an AWS Lambda function running your Java code. 
It will look like this:

![image]({{ site.baseurl }}/assets/img/posts/AWS-Cloudwatch-Lambda.png)
 
Finding the data
There is an API which presents data regarding football team fixtures. Some friendly developers at football-data.org created it and it's easy to access. Taking a quick look at their documentation, they offer an endpoint which shows all fixtures for a certain team.
The endpoint we will use is: http://api.football-data.org/v1/teams/{teamID}/fixtures
You'll need to replace {teamID} with the ID of our team. At the current time, their API doesn't give you a list of teams, but it does give you teams in a competition. You can get a list of competitions here:
http://api.football-data.org/v1/competitions/
At the time of writing, Sheffield Wednesday are in the Championship. A quick look down the response shows the Championship of 2016/17 has an ID of 427. You can then visit the endpoint to get all teams in the Championship, this will give you the ID for your team.
http://api.football-data.org/v1/competitions/427/teams
The response gives you a link to the teams fixtures in the _links block: http://api.football-data.org/v1/teams/345/fixtures
Calling this endpoint returns an array of Fixture resources. A 'Fixture' looks like this:
 
    {
      "_links": {
        "self": {
          "href": "http://api.football-data.org/v1/fixtures/151400"
        },
        "competition": {
          "href": "http://api.football-data.org/v1/competitions/427"
        },
        "homeTeam": {
          "href": "http://api.football-data.org/v1/teams/345"
        },
        "awayTeam": {
          "href": "http://api.football-data.org/v1/teams/63"
        }
      },
      "date": "2017-05-07T11:00:00Z",
      "status": "FINISHED",
      "matchday": 46,
      "homeTeamName": "Sheffield Wednesday",
      "awayTeamName": "Fulham FC",
      "result": {
        "goalsHomeTeam": 1,
        "goalsAwayTeam": 2
      },
      "odds": {
        "homeWin": 2.45,
        "draw": 3.5,
        "awayWin": 3.1
      }
    }
 
You can simply search this array for a fixture with today's date. Alternatively, you could parse the response and search for events in the future (perhaps you want your text message a day early?) or events against certain teams. However, for this project you probably just want a simple text message to say if they are playing today.
 
Creating the Java project
I'll give you an option here - if you are the kind of person who just wants the solution, it's available on github here (It requires you to modify your phone number, team ID, etc in MatchAlert.java):
https://github.com/NutterzUK/swfc-notifier/
If not, here is a brief explanation. The project is a maven project, which means you will want maven installed. You'll still want to check the codebase out from github, and you'll need Maven installed. You can find out how to install maven if you haven't already here: https://maven.apache.org/install.html
The project depends on two dependencies, both from amazon. They are the Amazon Lambda SDK and the Amazon SNS SDK. These are defined in the pom.xml. These dependencies also pull in the Apache HTTP library which is used for sending the request.
If you go into the source code, you'll see handleRequest() is the method which is executed by Amazon Lambda.
It uses the Apache HTTP Client to send a HTTP request to http://api.football-data.org/v1/teams/" + TEAM_ID + "/fixtures. For this to work, you will need to put in the right values at the top of the class. After it gets the response, it takes the body as a String. It then uses the Strings .contains() method to check if today's date is anywhere in the response. If it is, that means your team is playing, so it calls the method to send a text message, and if not, it does nothing.
Creating the lambda function
The first step is to go into the root of the java project, open a console window, and type:
    mvn clean install
This will package up the project as a .jar file which can be uploaded to Amazon Web Services.
The next step is to go to the AWS console and log in, then go to "Lambda" under Compute > Services.
In the top right, select a region closest to you. You should remember this region, as you'll need to know it to find your lambda function in future.
Select "Create a Lambda Function", then click Blank Function. For now we have no trigger, so press "Next".
Give your function a name, description, and select Java 8 as the runtime.
Upload the jar file. This can be found in the target directory of your project after you ran the "mvn clean install" command. It will be named: lambda-0.0.1.jar
For the handler, you need to specify the fully qualified class name and the method name, in this format:
io.nutbrown.matchalerts.lambda.MatchAlert::RequestHandler
For Role, select "Create a new role". This will open a new page. Name it something sensible, e.g "FootballNotifierLambdaRole".
Press "Allow".
Note that you have created a role, but that role currently does not allow you to send SMS messages via SNS.
Press Next, and your lambda function will be created, but as you do not yet have permission to send texts, it will not work.
In the top corner, select "Services" and go to IAM > Roles. Here you will find the role you have created. Select it, then click "attach policy". Search for "SNS" and click the check box next to AmazonSNSFullAccess. Press attach policy. Now, your lambda function which has this role will have full access over SNS. You can create more granular policies to say exactly what the role can do, but for now this is enough.
Go back to services > lambda. Select your lambda application and we are ready to test!
Click "Test". Delete the default JSON structure there, and enter a simple string with any input (We don't use it at the moment). It should look something like this:
 
Press test, and if your team is playing today, you should receive a notification. If not, the resulting log will say that a message was not sent.
 
CloudWatch events
The final step is to trigger your lambda function daily. For this, we will use CloudWatch events.
Click Services, and go to "CloudWatch" (Please don't confuse it with CloudFront or CloudFormation).
In the left navigation, click events.
Click "Create Rule".
Select "Schedule"
Enter the cron expression:
0 10 * * ? *
This cron expression runs at 10am every day. You can change this to easily decide when you want your lambda function to run.
If you want more information on the format of cron expressions, see: https://docs.aws.amazon.com/AmazonCloudWatch/latest/events/ScheduledEvents.html
Under targets, select your lambda function, and select Constant (JSON text). In here, enter your test input as you did before (e.g "test").
Everything should be ready to go!
 
Results!
Here is a screenshot of the final result, sent to me during testing on my Android phone:
 
Perhaps you can alter to code to give more information (e.g when they are playing) or you can give notifications on other events. Let me know in the comments what you manage to create or if you have any questions.

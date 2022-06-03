---
layout: post
title: "Push notifications with SNS and Firebase"
description: This post describe how I setup SNS to communicate with Firebase to send push notifications to Android and iOS devices
author-id: "galera"
categories: [aws, mobile, backend, java]
tags: [aws, sns, mobile, java, backend, firebase]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/sns-push/featured.jpg"
thumbnail: "assets/img/posts/sns-push/featured.jpg"
image: "assets/img/posts/sns-push/featured.jpg"
---

Currently we're investigating usage of SNS to send Push notifications. We have some friends working already with Firebase and recommend us to check it out. We discovered that SNS supports sending messages to Firebase. However, not everything is as easy as it looks like.

<p><!--more--></p>

## How push notifications work

There are quite a lot of moving parts in the scenario:

- AWS SNS: Messaging system from AWS
- Firebase Cloud Messaging (FCM): Message system to connect with the devices
- Device: device that receives push notification

The communication works this way:

1. The mobile device is registered within FCM.
2. FCM answers with a token that identifies the user. For the sake of testing
this will be manually extracted, however when doing this in PRD,
we need to automate this process.
3. In AWS SNS, a platform app has been created connecting AWS SNS with
   FCM with the provided API credentials (see section below)
4. WIth the token obtained in step 2). A new app endpoint is created.
   This endpoint identifies the app that registered the token
5. When the backend wants to send the push notification, it uses the
   registered app endpoint for that token.
6. SNS forwards the message to FCM
7. FCM sends the message to the mobile device

![Push notifications communication](/assets/img/posts/sns-push/pn.png "Push push notifications communication")

## Testing FCM device configuration

The preliminary test are done using an android application because I'm more used to Android development. Once FCM is setup, we could obtain the token from the device and try to send a push notification from FCM console to the device. If we do that, we receive a message with an structure similar to:

```json
{  
    "notification" : {
        "title": "Test notification",
        "body": "Test notification body"
    },
    "data": {

    }
}
```

So far, so good. We can do the same test from FCM to IOS and we will receive the same payload and the notification will pop up.

## Testing SNS FCM connection

After following this <a href="https://aws.amazon.com/es/premiumsupport/knowledge-center/create-android-push-messaging-sns/">guide</a>, now we have connected SNS with FCM.

We can try sending a message to the mobile device by specifying the app endpoint. We can access through AWS Console to SNS UI and search for app plattform and the app endpoint that belongs to our device.

If we go the default way and send a message, our Android device will receive a message similar to:

```json
{  
    "data": {
        "default": "Test message"
    }
}
```

Comparing the received data with the previous one, we see a fundamental difference: the notification field is empty and the message is inserted in the data field inside a `default` field.

If we repeat the same operation for an iOS device, we will not receive the push notification.

To send in this default way from Java, we could use the following code:

```java
public void publishDefaultMessage(String endpoint) {
    PublishRequest publishRequest = PublishRequest.builder()
            .message(getTextMessage())
            .targetArn(endpoint)
            .build();

    PublishResponse result = snsClient.publish(publishRequest);
    System.out.println(result.messageId() + " Message sent. Status was " + result.sdkHttpResponse().statusCode());
}
```

## Sending to Android and iOS

If we want to be able to send both to Android and iOS through FCM, we need to send a custom payload to SNS. If we send with default configuration, Android can receive the notification but not iOS.

To do so, we need to provide a JSON message as the payload:

```json
{"GCM":"{\"notification\":{\"title\":\"Title\",\"body\":\"Notification body sent with custom payload\"},
\"data\":{\"orderId\":\"1234\",\"customerId\":\"1234\"}}"}
```

We can send that test message from AWS Console, specifying to send different payload per protocol. Or we can do it through code:

```java
public void publishCustomMessage(String endpoint) {
    PublishRequest publishRequest = PublishRequest.builder()
            .message(customFirebaseMessage())
            .targetArn(endpoint)
            .messageStructure("json") // Send custom payload per transport type
            .build();

    PublishResponse result = snsClient.publish(publishRequest);
    System.out.println(result.messageId() + " Message sent. Status was " + result.sdkHttpResponse().statusCode());
}
private String customFirebaseMessage() {
    Map<String, String> customMessage = new HashMap<>();
    final String FIREBASE_PROTOCOL = "GCM";
    customMessage.put(FIREBASE_PROTOCOL, getFirebaseMessage());
    return new Gson().toJson(customMessage);
}
private String getFirebaseMessage() {
    FirebaseMessage message = new FirebaseMessage()
            .withTitle("Title")
            .withBody("Notification body sent with custom payload")
            .withDataEntry("customerId", "1234")
            .withDataEntry("orderId", "1234");
    return message.toJson();
}
```

Where `FirebaseMessage` is an object we have created:

```java
public class FirebaseMessage {
    private final Map<String, Object> notification = new HashMap<>();
    private final Map<String, Object> data = new HashMap<>();

    public FirebaseMessage withTitle(String title) {
        this.notification.put("title", title);
        return this;
    }

    public FirebaseMessage withBody(String body) {
        this.notification.put("body", body);
        return this;
    }

    public FirebaseMessage withDataEntry(String key, String value) {
        this.data.put(key, value);
        return this;
    }

    public String toJson() {
        return new Gson().toJson(this);
    }
}
```

If we send the messages with this format, they will be received both in Android and iOS
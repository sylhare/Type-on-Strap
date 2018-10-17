---
layout: post
title: "SOAP : Simple Object Access Protocol"
tags: [SOAP]
subtitle: "SOAP의 개념 및 동작원리"
feature-img: "md/img/SOAP/soap-communication.png"              
thumbnail: "md/img/SOAP/soap-communication.png"
excerpt_separator: <!--more-->
sitemap:
changefreq: daily
priority: 1.0
---

<!--more-->

# SOAP의 개념 및 동작원리

---

### SOAP

SOAP은 Web Services Interfacesd의 일종으로  컴퓨터 네트워크(HTTP, HTTPS, SMTP 등)상에서 XML 기반의 메시지를 교환하는 통신 프로토콜이다. 
HTTP 특성상 Proxy와 방화벽에 제약을 받지 않고 쉽게 통신할 수 있다. 보통의 경우 SOAP은 RPC 형태의 메시지 패턴을 사용하여 통신한다.

>RPC(Remote Procedure Call) - 원격 프로시져 호출
>
>별도의 원격 제어를 위한 코딩 없이 다른 주소 공간에 있는 함수나 프로시저를 실행할 수 있게 하는 프로세스 간 통신 기술 뜻한다. 즉 어떤 프로그램에 대하여 그 프로그램이 로컬에 있건 원격에 있건 사용할 수 있다.

---

### 왜, 언제 SOAP을 사용할까?

다음의 예시를 보자.

>A사는 Window 환경에서의 JAVA로 개발한 프로그램을 소유하고 있다. A는 검색엔진이 필요로 하여 검색엔진이 유명한 B사의 검색엔진 프로그램을 사용하려 한다.
>하지만 B사의 검색엔진은 Linux 환경에서 C#으로 개발한 프로그램이어서 분산된 환경에서 기능을 호출하기 어려움을 겪는다.
>
>이때 SOAP을 통해 기능을 연동한다.

여기서 우리는 SOAP에 대해 짐작할 수 있다.

#### Widow에서 Linux로, JAVA에서 C#으로 통신이 가능할까?

첫 번째, 분산된 환경에서 SOAP이다. 'Widow에서 Linux로, JAVA에서 C#으로 통신이 가능할까?'라는 의문에 답은 'SOAP이면 가능하다'이다.

<img src="/md/img/SOAP/soap-communication.png">
<em>SOAP의 통신</em>

이처럼 분산된 환경에서 데이터 통신이 가능한 이유는 SOAP은 네트워크(인터넷)상에서 데이터를 교환하기 때문이다.

#### 프로그래밍 언어가 다를 때 데이터 교환이 가능할까?

두 번째, SOAP은 네트워크상에서 메시지를 교환하여 통신한다. 이때 메시지는 XML 기반으로 하고 있는데 XML은 W3C에서 개발된 마크업 언어로 표준을 따른다. 표준을 따르기 때문에 플랫폼에 종속적이지 않고 독립적이다.
즉 프로그래밍 언어가 다른 두 프로그램이 데이터 교환이 가능하다.


#### 다양한 전송 프로토콜 지원

세 번째, SOAP은 HTTP뿐만 아니라 다양한 전송 프로토콜들의 사용을 허용한다. 표준 스택에서는 전송 프로토콜로 HTTP를 사용하지만, 다른 프로토콜 역시 사용할 수 있다.

---

### 단점
   
SOAP은 분산 환경에서의 기능을 호출할 수 있다는 장점이 있지만, 다음과 같은 단점들이 있다.

 - 오버헤드가 발생할 확률이 높다.
 - 높은 개발 난이도를 요구한다. 
 - 무겁고 느리다.

SOAP은 복잡한 구조를 띠고 있어 오버헤드가 발생할 확률이 높다. 또한, XML 기반의 메시지가 교환되는 과정에서 인코딩/디코딩은 필수적이라 웹 서비스 개발이 어려운 편이다. 메시지는 XML 기반이기 때문에 XML 특성상 다른 미들웨어 기술(REST와 같은)보다 상대적으로 무겁고 속도도 느리다. 다만 전송할 메시지가 적을 때는 문제 되지 않을 수 있다.

>최근 XML의 성능을 향상시키기 위해, <br/>
>VTD-XML과 같은 emerging non-extractiv XML 처리 모델이 있다.
>바이너리 객체를 포함시킨 특별한 경우의 XML(바이너리 XML을 말하는듯)로 메시지 전송 최적화 메커니즘(Message Transmission Optimization Mechanism; MTOM)이 있다. 	

---

### 동작원리

<img src="/md/img/SOAP/soap-architecture.png">
<em>SOAP 아키텍처</em>

SOAP의 동작원리는 다음과 같다.

웹 서비스 제공자(`Serivce provider`)는 웹 서비스를 WSDL로 인코딩하여 WSDL을 관리하는 UDDI 레지스트리에 등록(`Publish`)한다. 웹 서비스 요청자(`Service requestor`)는 사용할 웹 서비스를 WSDL로 인코딩하여 UDDI 레지스트리를 통해 탐색(`Find`)한다. 
이때 해당 WSDL이 있다면 바인딩(`Bind`)하고 인코딩된 WSDL 데이터를 디코딩하여 해당 웹 서비스 로직(기능)을 수행한다.

---

### XML 메시지

<img src="/md/img/SOAP/soap-message-structural.png">
<em>SOAP 메시지 구조</em>

SOAP 메시지는 XML을 근간으로 `<Header>`와 `<Body>`를 포함한 `<Envelope>`로 구성된 디자인 패턴으로 설계되어 있다.


### `<Envelope>`
모든 SOAP 메시지의 루트 요소이며 두 개의 하위 요소인 `<Header>`요소 와 `<Body>`요소를 포함한다.

### `<Header>`
`<Envelope>`의 선택적인 하위 요소이다. 반복이나 보안 및 트랜잭션을 정보로 하는 메타 정보를 가지고 있고, 메시지 플로우에 따라 SOAP 노드로만 처리될 애플리케이션 관련 정보를 전달하는데 사용된다.

### `<Body>`
`<Envelope>`의 필수적인 하위 요소이다. 메시지의 최종 수신인을 대상으로 하는 주요한 정보(호출 및 응답에 필요한 정보)를 포함하고 있다.

|  Entry | 하위요소 | <center>Data</center>|<center>Description</center> |
|:---:|:---:|:---|:---|
|`<Header>`|선택| 메타정보|반복이나 보안 및 트랜잭션을 정보로 하는 메타 정보|
|`<Body>`|필수|  호출 및 응답 정보|메시지의 최종 수신인을 대상으로 하는 정보|


#### xml sample code
```xml
<soap : Envelope xmlns : soap = '.....' 
                   xmlns : rep = '....' 
                   xmlns : xmlmime = '....'>
                   
	<!-- Header -->
	<soap : Header>
    		<!-- Header Block Start -->
        	<rep : Representation resource = '....'>
          	<rep : Data xmlmime : contentType = 'img/png'>
          		...
          	</rep : Data>
        	</rep : Representation>
        	<!-- Header Block End -->
	</soap : Header>
	
	<!-- Body -->
	<soap : Body>
	        <x : MyData xmlns : x = '....'>
		          <x : name> 6 clock </ x : name>
		          <x : img src = '....'/>
	        </x : MyData>
    </soap : Body>
</soap : Envelope>
```

---

### 메시지 - Header

<img src="/md/img/SOAP/soap-message-structural-header.png">
<em>SOAP 메시지 - Header 구조</em>

#### Header - Header-Block

`<Header>`의 바로 아래 하위 요소를 `<Header-Block>`이라 한다.<br/>
SOAP 메시지 내부에 같은 메타정보를 포함한 `<Header>`가 여러 존재할 수 있는데  구분하기 위해 `<Header-Block>`를 통해 각각의 `<Header>`를 정의할 수 있다.

`<Header-Block>`은 요청자와 응답자 간 메시지 교환 과정에서 발생할 수 있는 모든 정보(엑세스 제한, 참조정보 등)를 정의한다.

`<Header-Block>`은  SOAP 중개 노드와 최종 SOAP 수신자 노드로 처리할 수 있다. 그러나 실제 애플리케이션에서는 모든 노드가 모든 `<Header-Block>`을 처리하는 것은 아니다. 각 노드는 일반적으로 특정 `<Header-Block>`을 처리하도록 디자인되고 각 `<Header-Block>`은 특정 노드에 의해 처리된다.

>***SOAP 노드** <br/>
>SOAP 노드는 웹 서비스 처리가 구성되고 적용되는 플로워에서 지점의 역할을 한다. - [IBM SOAP 노드](https://www.ibm.com/support/knowledgecenter/ko/SSMKHH_10.0.0/com.ibm.etools.mft.doc/ac55850_.htm)

```xml
	<!-- Header -->
	<soap : Header>
    		<!-- Header Block Start -->
    		<!-- 메시지의 메타 정보 -->
        	<rep : Representation resource = '....'>
          	<rep : Data xmlmime : contentType = 'img/png'>
          		...
          	</rep : Data>
        	</rep : Representation>
        	<!-- Header Block End -->
	</soap : Header>
```

#### Header - 제어정보

제어정보는 기능을 처리할 수 있는 속성이 필수인지 또는 선택적인지 정의한다. 이러한 제어 정보에는 지시문 전달 또는 메시지 처리에 관련된 컨텍스트 정보 등이 포함된다.

 _**Header 제어정보를 통해 SOAP 메시지를 통신 당사자 간 사전 합의 없이도 분산 방식으로 SOAP 메시지에 기능을 확장할 수 있다.**_

---

### 메시지 - Body

`<Body>` 요소와 그에 연관된 하위 요소는 요청자와 응답자 간 데이터를 교환하는 데 사용된다.

<img src="/md/img/SOAP/soap-message-structural-body.png">
<em>SOAP 메시지 - Body 구조</em>

#### Body - Fault

SOAP은 `<Body>`에 대해 한 개의 하위 요소인 `<Fault>`를 정의하며 오류 보고에 사용된다. `<Fault>`는 `<Body>`의 항목으로 표시되어야 하고 `<Body>`에 두번 이상 표시되어서는 안된다.

`<Fault>`의 하위 요소는 SOAP1.1에서와 SOAP1.2에서 다르므로 자세한 내용은 [IBM SOAP Fault](https://www.ibm.com/support/knowledgecenter/ko/SSMKHH_10.0.0/com.ibm.etools.mft.doc/ac55810_.htm)를 참고하자.

---

### 올바른 사용

SOAP는 인터넷 애플리케이션 계층에 있는 프로토콜을 전송계층의 프로토콜로 사용할 수 있게 만든다.<br/>
**SOAP의 반대 여론**들은 프로토콜의 의도된 목적과 역할이 맞지 않아 부정 이용이 된다고 비판하지만, **SOAP의 지지 여론**들은 터널링을 위한 다양한 계층(level)에 쓰이고 있는 다른 프로토콜들과 비슷하다고 말하고 있다. <br/>

<img src="/md/img/SOAP/network-layer.png">
<em>Network Layer</em>

SOAP의 올바른 이용은 SMTP와 HTTP에서 애플리케이션 계층 프로토콜로 트랜스포트 계층의 역할을 대신하는 것이다. 
이러한 이유는 HTTP는 오늘날 인터넷 인프라와 매우 잘 동작하여 더욱 폭넓은 지원을 가능하기 때문이다. 특히나, HTTP는 방화벽이 작동하는 네트워크 안에서도 문제 없이 작동한다.

SOAP는 암호화된 HTTPS(애플리케이션 계층에서는 HTTP와 동일하나 트랜스포트 계층 아래에서는 암호화됨)에서도 간략하게 또는 상호적으로 사용된다.

---

### REST(Representational State Transfer)

많은 개발자는 ‘SOAP(RPC)은 죽고 미래는 RESTFUL(REST의 호환 시스템)에 있다.’라고 말한다. 이러한 이유는 REST는 SOAP이 복잡해지면서 새롭게 개선한 웹 서비스 인터페이스이기 때문이다.

REST는 리소스 기반 (또는 데이터베이스)을 중점으로 작업하고 HTTP에서 기능(GET, PUT, POST, DELETE)을 상속받아 작업한다. 표준 HTTP 통신의 특성상 대부분 방화벽에 대한 제약이 없다. 프로토콜에 오버헤드가 없으므로 가볍고 XML 기반의 데이터를 사용하는 SOAP과 달리 REST는 일반적으로 JSON 형태의 데이터를 사용하기 때문에 빠르다. 이러한 Simplicity는 공개 API를 다룰 때 REST를 선호하는 경우가 많고 Amazon 및 Google과 같은 대기업이 API를 SOAP에서 REST로 옮기는 가장 큰 이유 중 하나이다.

---

### 참고

> [https://www.ibm.com/support/knowledgecenter/en/SSGMCP_5.3.0/com.ibm.cics.ts.webservices.doc/concepts/dfhws_definition.html](https://www.ibm.com/support/knowledgecenter/en/SSGMCP_5.3.0/com.ibm.cics.ts.webservices.doc/concepts/dfhws_definition.html) <br/>
> [https://ko.wikipedia.org/wiki/SOAP](https://ko.wikipedia.org/wiki/SOAP)<br/>
> [https://www.codecademy.com/articles/what-is-rest](https://www.codecademy.com/articles/what-is-rest) <br/>
> [https://docs.oracle.com/cd/E12890_01/ales/docs32/webservicesprogrammersguide/webservice_api.html](https://docs.oracle.com/cd/E12890_01/ales/docs32/webservicesprogrammersguide/webservice_api.html)<br/>
> [https://medium.freecodecamp.org/rest-is-the-new-soap-97ff6c09896d](https://medium.freecodecamp.org/rest-is-the-new-soap-97ff6c09896d)<br/>
> [https://stormpath.com/blog/rest-vs-soap](https://stormpath.com/blog/rest-vs-soap)<br/>
> [https://www.infoq.com/articles/rest-soap-when-to-use-each](https://www.infoq.com/articles/rest-soap-when-to-use-each)<br/>
> [https://www.tutorialspoint.com/soap/what_is_soap.htm](https://www.tutorialspoint.com/soap/what_is_soap.htm)<br/>
> [https://searchmicroservices.techtarget.com/definition/SOAP-Simple-Object-Access-Protocol](https://searchmicroservices.techtarget.com/definition/SOAP-Simple-Object-Access-Protocol)<br/>

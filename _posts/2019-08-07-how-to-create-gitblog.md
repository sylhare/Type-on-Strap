---
layout: post
title: "Github.io 깃블로그 만들기- 1)"
---

###  Github Blog 시작하기

개발 공부를 시작하고 많은 지식을 단기간에 쌓게 되며, 여러 사이트와 블로그를 참고했었다. 정리가 잘 된 블로그를 보며 공부를 하거나 에러를 해결하게 될 때, 속이 뻥 뚫리는 후련함과 동시에 고마움을 겪었던 적이 많다. 어느덧 개발자로 커리어를 시작하고 하나하나 커리어를 쌓아가고 있는 요즘, 나와 같은 문제를 겪게될 다른 사람들에게 나도 도움을 주고 싶다는 생각과 동시에 내가 실제로 부딪히며 배운 지식들을 흔적처럼 남겨놓고 싶다는 생각을 하게 되었다. 

<strong>이러한 흔적들이 모여 훗날 더 멋진 개발자가 되어있길 바라며! </strong>

그래서 내가 선택한 도구는 깃블로그이다. 이유는 단순하다. 만들기 간편하고, 다양한 Jekyll theme이 있기 때문에 원하는 디자인을 가져다 쓰면 되기 때문에 디자인에 대한 고민보다 컨텐츠를 채우는 것에 집중할 수 있기 때문이다. 

*물론 Vue 혹은 React로 멋진 디자인의 개인 페이지를 만드는 것도 하나의 목표이긴하나.. 그것은 나중을 위해 아껴두고..*

아무튼 이러한 이유로 블로그를 시작하기로 결심했고, 방금 만들었으니! 첫번째 글은 자연스럽게 github 블로그 만드는 방법에 대한 소개이다.

![](https://media.giphy.com/media/4nfa3HlyeEef6/giphy.gif)



### STEP 01. Jekyll Theme Repository Fork 

원하는 Jekyll 테마를 받기 위해 구글에 "Jekyll Theme"으로 검색 후 원하는 테마를 고른다. 

![image-20190807195538574](/assets/img/images/jekyll.png)

나는 다양한 테마 중 **Type on Strap** 을 선택했다. 각자 원하는 테마를 선택했다면, 클릭 후 "Homepage" 버튼을 클릭한다. 그럼 git repository로 연결될텐데, 

![image-20190807200235978](/assets/img/images/git-repository-fork.png)

여기서 fork 버튼을 눌러주면 소스코드를 본인의 repository로 가지고 올 수 있다. 

> 현재 모든 과정은 본인의 github 계정을 가지고 있다는 가정 하에 진행되고 있으니, 없다면 지금 만들고 다시 따라하면된다.

그럼, 본인의 repository 목록을 확인해보면, fork된 것을 확인할 수 있다. 





### STEP 02. Change Repository name

git blog를 만들기 위해서는 repository 이름을 적절히 바꿔주어야하는데, 이를 위해서 fork된 repository에서 "Settings"에 들어간다. 

![gihub-repname](/assets/img/images/gihub-repname.png)

깃블로그를 만들기위해 설정해야 하는 repository이름에는 정해진 규칙이 있다. ( **아이디.github.io** ) 따라서 이에 맞춰 Rename을 해주면 끝! 본인 아이디와 동일하게 맞춰주면 된다. 

 그 뒤, 같은 Settings 페이지에서 스크롤을 조금 내린 뒤 새로고침을 몇 번 하다보면, 다음과 같이 블로그가 publish 됐다는 초록색의 알림이 뜬다. 

![gitblog-ready](/assets/img/images/gitblog-ready.png)

이제 한 가지의 설정만 바꿔주면, 나만의 깃블로그가 만들어진다! 

> publish된 주소로 들어가도 html/css가 적용이 안 돼 모두 깨져보일 것이다. 당황하지말고 STEP 03를 따라해보자. 





### STEP 03. Edit <code>_config.yml</code>

웹의 디자인이 이상했던 것은, css파일 경로가 제대로 잡혀있지 않아서 발생한 문제이다. 이는 <code>_config.yml</code> 이라는 환경설정 파일에 들어가 url을 수정해주면 해결된다. 

![config-loc](/assets/img/images/config-loc.png)

이 곳으로 들어가, **baseurl** 부분을 찾는다. 각 테마에 따라 config 파일이 다르기 때문에, baseurl의 위치는 각자 다를 수 있다. 

![config-edit](/assets/img/images/config-edit.png)

파일을 수정하기 위해서는 연필 모양을 클릭하면 된다. 내 경우 baseurl이 맨 윗부분에 위치해있었고, 여기 써있던 경로를 지워준 뒤 저장("Commit changes")하면 된다. 



그리고나서, 다시 publish된 주소로 접속하면! 

![theme-screenshot](/assets/img/images/theme-screenshot.png)

#### Ta-dah-! 내가 선택했던 테마가 잘 적용된 나만의 블로그가 만들어졌다! 뿌-듯

![](https://media.giphy.com/media/OcZp0maz6ALok/giphy.gif)



틀은 모두 만들어졌으니, 이제 내가 채우고싶은 것들로 채우면 된다. 생각보다 글이 너무 길어질 것 같아, 이 내용은 바로 다음 글에 적도록 하겠다-! 




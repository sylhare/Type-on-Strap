---
layout: post
title: "GateXcanner FrontEnd Azure Pipeline 을 이용한 자동빌드"
author: Shinwusub
tags: [GateXcanner, FrontEnd]
---


GateXcanner FrontEnd를 자동 빌드 및 중앙 버전관리를 적용하기 위하여 Azure Pipeline을 사용해보았습니다.
참고로 GateXcanner 는 Docker 환경이 아닙니다.

# 목표

- main branch로 pull-request요청이 있을때 빌드 유효성 검사 및 자동 빌드와 버저닝이 되어야 한다.

---

## 세부 요구사항

### 버저닝

- 시맨틱 버저닝을 따라야한다. (ex : 1.0.0-beta+test)
- 사용자의 입력이 없다면 기본적으로 패치 버전만 증가한다 (ex : 1.0.0-beta+test, 1.0.1-beta+test, 1.0.2-beta+test, 1.0.3-beta+test...)
- 버전은 고유해야한다(중복불가).
- package.json 의 version 항목에 위 사항들이 적용되어야 한다.
- package.json 의 version 변경 시 git tag에 명시되어야 한다.

### 유효성 검사

- main branch로 pull-request요청이 있을때 해당 요청이 정상적으로 build되는 파일인지 유효성 검증이 가능해야한다.

---

# Pipelines
pipeline은 특정 레포지토리의 변경을 감지하여 개발자가 사전에 정해놓은 동작을 수행하게 해줍니다.

## azure-pipelines.yml
azure-pipelines.yml 는 main branch 에 pull-request요청이 발생 후 해당 pull-request가 Commit 되었을 때 자동으로 빌드 및 버전관리를 해줍니다.
```yml
# main branch 의 변경을 감지하여 변경이 있을때 pipeline이 작동합니다.
trigger:
  - main

# pipeline을 통하여 빌드, 테스트 등이 이루어질 일종의 Test PC환경의 Image를 MS에서 제공하는 Cloud Agent Image중 "ubuntu-latest" 라는 label값과 일치하는 값으로 셋팅합니다.
# MS에서 제공중인 Image들의 정보는 https://docs.microsoft.com/en-us/azure/devops/pipelines/agents/hosted?view=azure-devops&tabs=yaml 에서 확인할 수 있습니다.
pool:
  vmImage: ubuntu-latest

# 순차적으로 (위에서 아래) 실행될 프로세스들을 정의하는 영역입니다.
steps:
  # checkout 값을 self로 하면 현재 감지된 branch로 checkout합니다.
  - checkout: self
    # 스크립트(pipeline)이 git OAuth Token에 접근할 수 있도록 합니다.
    persistCredentials: true

  # 스크립트(pipeline)내부에서 git을 사용할 때 입력될 유저 정보를 설정합니다.
  # pipeline에서는프라이빗한 개인적인 작업이 아닌  테스트, 빌드, 버저닝 등 공개적인 작업이 이루어지기에 개발자 개인의 정보가 아닌 제품개발1팀 에 대한 정보를 기입하였습니다.
  - bash: |
      git config --global user.email "cb1-1@softcamp.co.kr"
      git config --global user.name "Product Development Team 1"
    displayName: "set git user.email & user.name"

  # ./templates/build.yml 경로의 build.yml template를 실행합니다.
  - template: ./templates/build.yml
    # build.yml template를 호출할때 useVersioning 파라미터에 값을 넘깁니다.
    parameters:
      useVersioning: true

  # ./templates/publish-build-artifacts.yml 경로의 publish-build-artifacts.yml template를 실행합니다.
  - template: ./templates/publish-build-artifacts.yml
    # publish-build-artifacts.yml template를 호출할때 contests 와 artifactName 이라는 파라미터에 값을 넘깁니다.
    parameters:
      contests: "dist/**"
      artifactName: "GateXcannerWeb"
```
</br>

## pr-build-validation.yml
pr-build-validation.yml 는 main branch 에 pull-request요청이 발생한 시점에 해당 pull-request가 실제로 빌드 가능한 유효한 코드들인지를 실제 build하여 확인합니다.
```yml

# build-validation 은 trigger를 통하여 감시 대상을 지정하지 않습니다.
# azure devops의 Project settings 탭에서 build validation UI를 통하여 별도로 지정해야합니다.
trigger: none

pool:
  vmImage: ubuntu-latest

# pr-build-validation.yml 의 목적 자체가 정상적으로 빌드가 되는 소스들인지를 Commit전에 확인하기 위함이기에 build 만 하는 build.yml template를 사용합니다.
# 버전관리, 빌드 결과물 퍼블리싱 등은 진행하지 않습니다.
steps:
  - checkout: self
    persistCredentials: true

  - template: ./templates/build.yml
    parameters:
      useVersioning: false
```
### build validation 설정

![build validation settings](/assets/img/posts/2022-05-19-GateXcanner-FrontEnd-Azure%20Pipeline-%EC%9D%84-%EC%9D%B4%EC%9A%A9%ED%95%9C-%EC%9E%90%EB%8F%99%EB%B9%8C%EB%93%9C/build-validation-%EC%84%A4%EC%A0%95%EB%B2%95.png)
![bypass policies when pushing](/assets/img/posts/2022-05-19-GateXcanner-FrontEnd-Azure%20Pipeline-%EC%9D%84-%EC%9D%B4%EC%9A%A9%ED%95%9C-%EC%9E%90%EB%8F%99%EB%B9%8C%EB%93%9C/build-validation-%EB%AC%B4%EC%8B%9C-%EA%B6%8C%ED%95%9C-%EB%B6%80%EC%97%AC.png)

#### 주의 사항
- main branch에 build validation 을 설정하면 main branch 로 pull request를 요청하는 경우의 모든 pull request는 build validation을 거치게 됩니다. (ex : feature/test -> main, dev -> main, hotfix -> main 모두 빌드 유효성 검사 대상)

- 특정 branch에 build validation 을 설정하면 push는 기본적으로 불가능하며, push가 필요한 경우 push를 진행하는 계정에 대하여 build validation 을 무시할 수 있는 권한을 부여해야합니다.

- Build expiration을 "Immediately when main is updated" 로 선택하였다면 pull-request의 from이 되는 branch가 push등으로 업데이트된다면 즉시 build validation 을 다시 진행합니다.

---

# templates
templates는 pipeline에서 목적을 가진 일정 부분을 별도의 yml 파일로 분리하여 컴포넌트 개념과 같이 여러 pipeline에서 가져다 사용할 수 있도록합니다.

## build.yml
build.yml 는 npm build와 versioning을 해주는 template입니다.
```yml
parameters:
  - name: useVersioning
    type: boolean
    default: false
  - name: packageVersion
    type: string
    default: 1.0.0
  
steps:
  # node 설치 template을 사용합니다.
  - template: ./install-node.yml
    # 설치 버전을 파라미터로 넘깁니다.
    parameters:
      versionSpec: "14.x"

  # npm 설치 template을 사용합니다.
  - template: ./install-npm.yml

  # npm version 명령어는 기본적으로 git tag도 생성해줍니다.
  # 하지만 원하는 타이밍에 tag를 생성하기 위해 tag 생성 기능을 제외시킵니다.
  - ${{ if eq(parameters.useVersioning, true) }}:
      - bash: |
          npm --no-git-tag-version version patch
        displayName: "version up"

  - bash: |
      npm run build
    displayName: "npm build"

  # package.json 의 version 정보를 가져옵니다.
  - ${{ if eq(parameters.useVersioning, true) }}:
      - bash: |
          v=`node -p "const p = require('./package.json'); p.version;"`
          echo "##vso[task.setvariable variable=packageVersion]$v"
          echo "$(packageVersion)"
        displayName: "get version"

  # 버전을 올리면서 변경된 package.json을 commit & push 합니다.
  - ${{ if eq(parameters.useVersioning, true) }}:
      - task: CmdLine@2
        inputs:
          script: |
            git add package.json
            git commit -a -m "***NO_CI***Azure Pipeline Auto Build : $(packageVersion)"
            git push -u origin HEAD:main
        displayName: "commit & push the version-changed package.json file"

  # 변경 된 버전 정보를 tag로 생성 및 push 합니다.
  - ${{ if eq(parameters.useVersioning, true) }}:
      - task: CmdLine@2
        inputs:
          script: |
            git tag "$(packageVersion)"
            git push -u origin HEAD:main "$(packageVersion)"
        displayName: "tag creation and version change push"
```

## install-node.yml
install-node.yml 는 node 설치를 해주는 template 입니다.
```yml
parameters:
  - name: versionSpec
    type: string
    default:
      - default

steps:
  - task: NodeTool@0
    inputs:
      # parameters.versionSpec 에 명시된 버전으로 설치합니다.
      versionSpec: ${{ parameters.versionSpec }}
    displayName: "install node.js"
```

## install-npm.yml
install-npm.yml 는 npm 설치를 해주는 template 입니다.
```yml
steps:
  - bash: |
      npm install
    displayName: "install npm"

  - bash: |
      npm i -g @vue/cli-service
    displayName: "install vue-cli-service"

```

## publish-build-artifacts.yml
publish-build-artifacts.yml 는 build 결과물(artifacts)을 publish 하기 위한 template 입니다.
```yml
parameters:
  - name: contests
    type: string
    default: "dist/**"
  - name: artifactName
    type: string
    default: "GateXcannerWeb"

steps: 
# build 결과물 경로인 parameters.contests 를 Build.ArtifactStagingDirectory 로 복사합니다.
- task: CopyFiles@2
  inputs:
    Contents: ${{ parameters.contests }}
    TargetFolder: "$(Build.ArtifactStagingDirectory)"
  displayName: "copy dist dir to artifactstagingdirectory"

- task: CopyFiles@2
  inputs:
    sourceFolder: "$(Build.SourcesDirectory)"
    contents: "package.json"
    targetFolder: $(Build.ArtifactStagingDirectory)
  displayName: "copy package.json"

# Build.ArtifactStagingDirectory 를 artifact로 publish 합니다.
- task: PublishBuildArtifacts@1
  inputs:
    PathtoPublish: "$(Build.ArtifactStagingDirectory)"
    ArtifactName: ${{ parameters.artifactName }}
  displayName: "publish artifact: ${{ parameters.artifactName }}"

- bash: |
    cd $(Build.ArtifactStagingDirectory)
    echo $(Build.ArtifactStagingDirectory)
  displayName: change current directory

```
---
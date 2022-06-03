---
layout: post
title: Standalone application for Raspberry Pi
description: The design, implementation and continuous integration of a standalone application to be run on a Raspberry Pi. 
author-id: "galera"
categories: [raspberry-pi, python, ci/cd]
tags: [linux,nodejs,raspberry-pi,devops,docker]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/standalone-raspberry/featured-image.jpg"
thumbnail: "assets/img/posts/standalone-raspberry/featured-image.jpg"
image: "assets/img/posts/standalone-raspberry/featured-image.jpg"
---
<p>I'm building a small application to give treats to my dog in a remote manner. </p>

<p>I setup a Raspberry Pi with a very basic HTTP server connected to a servo motor that will open or close the deposit where the treats are stored. </p>

<p>In this article I'll explain all the challenges I found to make this application standalone.</p>
<p><!--more--></p>

<h2>Requirements</h2>

- Accessible via web
- To have a camera, a button to give treats and a button to play a sound
- Easily installable in a Raspberry Pi: no need to install trillions of libraries
- Production ready: even though this is a personal project, I want the app to be 100% test covered and to have a full CI/CD cycle

<h2>Solutions</h2>

First of all, I did some small investigations and tackle every requirement in a separate way. This way I manage to found scripts that:

- <a href="https://picamera.readthedocs.io/en/release-1.13/recipes2.html#web-streaming">Create a MJPEG stream out of the Raspberry Pi camera</a>
- <a href="https://raspberrypi.stackexchange.com/questions/22708/is-there-some-trick-to-getting-aplay-audio-output-working">Play a sound from the disk</a>
- <a href="https://peppe8o.com/how-to-remote-control-a-servo-motor-sg90-using-raspberry-pi-zero-w-with-python-and-arduino/">Interact with a servo motor</a>

<h3>Backend</h3>

Once the parts are working independently, I made a python project with a very basic HTTP server based on `BaseHTTPRequestHandler` that receive request to the stream, to interact with the servo and to play a sound.

The interesting thing here was to be able to develop this project without using the Raspberry Pi. This is challenging because the required libraries are hardware specific to the Raspberry. But I manage to mock the camera and the servo libraries by using `unittest` python package

```python
from unittest.mock import MagicMock, patch


def mock_rpi_gpio():
    MockRPi = MagicMock()
    modules = {
        "RPi": MockRPi,
        "RPi.GPIO": MockRPi.GPIO,
    }
    patcher = patch.dict("sys.modules", modules)
    patcher.start()


def mock_pi_camera():
    picamera = MagicMock()
    modules = {"picamera": picamera, "picamera.PiCamera": picamera.PiCamera}
    patcher = patch.dict("sys.modules", modules)
    patcher.start()


mock_rpi_gpio()
mock_pi_camera()
```

`unittest` module allows you to define a `conftest.py` file that will be executed as a configuration step for you unit tests. Having done that, we can have tests that covers all the required functionality, even without installing the required libraries:

```python
from unittest.mock import call

import RPi.GPIO as mockGPIO

from dogfeeder.servo import Servo


def test_initialize_closed_servo():
    Servo()
    mockGPIO.setmode.assert_called_once_with(mockGPIO.BCM)
    mockGPIO.setup.assert_called_once_with(Servo.SERVO_PIN, mockGPIO.OUT)
    mockGPIO.PWM.assert_called_once_with(Servo.SERVO_PIN, 50)
    mock_pwm = mockGPIO.PWM()
    mock_pwm.start.assert_called_once_with(Servo.CLOSED)
```
<h3>Frontend</h3>

The implementation of the frontend is super simple. I used React to create 3 components:

- CallButton: the button that plays an audio file
- DispenseTreat: the button that interacts with the servo
- WebcamContainer: the img that prints the MJPEG stream out of the Pi Camera.

When any button is pressed and API call to backend is done in the background.

Nothing really fancy to see here.

<h3>CI/CD</h3>

When all the logic is done and the tests are passing, I decided that I wanted to go full professional and create a CI/CD pipeline for the project. In order to do that, I used gitlab.com

This has been the most challenging piece of the project. I wanted to create a standalone application so the installation process is keep to the minimum bar. In order to do so, I created a docker image with all the required dependencies to be used by Gitlab pipeline.

<h4>Docker image</h4>

```
FROM balenalib/raspberrypi3-python:3.7-buster
RUN apt update && apt upgrade
RUN apt install build-essential binutils zlib1g-dev
RUN apt install python3-picamera python3-rpi.gpio
RUN pip3 install pyinstaller pytest pytest-cov flake8 requests
```

It's based on <a href="https://github.com/balena-io-library/base-images/tree/master/balena-base-images/device-base/raspberrypi3">balenalib/raspberrypi3-python</a> Docker image, that simulates even the hardware and processor architecture of the Raspberry Pi 3. The docker image also contains all the libraries required to work (picamera, gpio, ...) and tools for the CI/CD (pytest, flake8).

`pyinstaller` is installed in order to generate the executable file of the backend


<h4>Pipeline</h4>

The pipeline contains four stages:
- test: unit tests of the backend and frontend
- release: to generate semantic versioned tags of the project
- build: to generate the standalone executable file of the backend and the web site for the frontend. Thanks to <a href="https://threedots.tech/post/automatic-semantic-versioning-in-gitlab-ci/">https://threedots.tech/post/automatic-semantic-versioning-in-gitlab-ci/</a>
- publish: I decided to store the generated artifacts within the Gitlab Package Registry

```yaml
stages:
  - test
  - release
  - build
  - publish

test-backend:
  image: registry.gitlab.com/adrian.galera/dogfeeder/python-ci
  stage: test      
  script: 
    - cd dogfeeder-backend
    - pytest --cov --cov-fail-under=100
  only:
    - master
    - branches

test-frontend:
  image: node:12.13-alpine
  stage: test
  script:
    - cd dogfeeder-web
    - npm ci --cache .npm --prefer-offline
    - npm test
  cache:
    key: "node-modules"
    paths:
      - .npm/
  only:
    - master
    - branches

release:
  image: python:3.7-stretch
  stage: release
  before_script:
    # Allow gitlab runner push code to gitlab.com
    # see: https://threedots.tech/post/automatic-semantic-versioning-in-gitlab-ci/
    - mkdir -p ~/.ssh && chmod 700 ~/.ssh
    - ssh-keyscan gitlab.com >> ~/.ssh/known_hosts && chmod 644 ~/.ssh/known_hosts
    - eval $(ssh-agent -s)
    - ssh-add <(echo "$SSH_PRIVATE_KEY")
    - pip install semver
  script:
    - python3 gen-semver
  only:
    - master
  when: manual

build-backend:
  image: registry.gitlab.com/adrian.galera/dogfeeder/python-ci
  stage: build
  script: 
    - cd dogfeeder-backend
    - pyinstaller dogfeeder/main.py -F --name dogfeeder-server
  only:
    - tags        
  artifacts:
    paths:
     - "dogfeeder-backend/dist/dogfeeder-server"

build-frontend:
  image: node:12.13-alpine
  stage: build
  script:
    - cd dogfeeder-web
    - npm ci --cache .npm --prefer-offline
    - npm run build
    - npm run zip
  cache:
    key: "node-modules"
    paths:
      - .npm/
  artifacts:
    paths:
     - "dogfeeder-web/dogfeeder-web_.zip"
  only:
    - tags  

publish:
  image: curlimages/curl:latest
  stage: publish
  script:
   - VERSION=${CI_COMMIT_REF_NAME}
   - 'curl --header "JOB-TOKEN: $CI_JOB_TOKEN" --upload-file dogfeeder-backend/dist/dogfeeder-server "${CI_API_V4_URL}/projects/${CI_PROJECT_ID}/packages/generic/dogfeeder/${VERSION}/dogfeeder-server"'
   - 'curl --header "JOB-TOKEN: $CI_JOB_TOKEN" --upload-file dogfeeder-web/dogfeeder-web_.zip "${CI_API_V4_URL}/projects/${CI_PROJECT_ID}/packages/generic/dogfeeder/${VERSION}/dogfeeder-web.zip"'  
  only:
    - tags
```

<h2>Installation</h2>

Now that the packages are stored in Gitlab, the installation is super simple. I created a script that downloads the artifacts from Gitlab and unzip the web into a running nginx and replace the executable file that will be picked up from a `supervisorctl` process:

```bash
VERSION=$1
TOKEN=${GITLAB_ACCESS_TOKEN}

if [ -z "$1" ]; then
    echo "You need to provide PACKAGE_VERSION as argument: sudo ./install-dogfeeder.sh <PACKAGE_VERSION>"
    exit 1
fi

if [ -z "$TOKEN" ]; then
    echo "You need to set GITLAB_ACCESS_TOKEN environment variable"
    exit 1
fi

wget --header "PRIVATE-TOKEN: ${TOKEN}" "https://gitlab.com/api/v4/projects/24187261/packages/generic/dogfeeder/${VERSION}/dogfeeder-server" -O /tmp/dogfeeder-server-${VERSION}
wget --header "PRIVATE-TOKEN: ${TOKEN}" "https://gitlab.com/api/v4/projects/24187261/packages/generic/dogfeeder/${VERSION}/dogfeeder-web.zip" -O /tmp/dogfeeder-web-${VERSION}.zip

unzip -o /tmp/dogfeeder-web-${VERSION}.zip -d /var/www/html/.
# Kill the process and supervisorctl will start it again:
ps -eaf | grep "dogfeeder-server" | grep -v grep | awk '{ print $2 }' | xargs kill -9 && cp /tmp/dogfeeder-server-${VERSION} /home/pi/.local/bin/dogfeeder-server
chmod +x /home/pi/.local/bin/dogfeeder-server
```
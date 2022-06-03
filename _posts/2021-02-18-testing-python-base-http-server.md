---
layout: post
title: Testing python BaseHttpServer
description: How do I manage to put BaseHttpServer and mock dependencies
author-id: "galera"
categories: [python, testing]
tags: [python, testing, mocking]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/test-python-basehttpserver/featured-image.jpg"
thumbnail: "assets/img/posts/test-python-basehttpserver/featured-image.jpg"
image: "assets/img/posts/test-python-basehttpserver/featured-image.jpg"
---
<p>While the development of <a href="https://www.agalera.eu/standalone-app-raspberry-pi">https://www.agalera.eu/standalone-app-raspberry-pi/</a> I needed to use python's BaseHttpServer and inject some dependencies into it.</p>

<p>It turns out, there's no easy way of doing that. Moreover, I wanted to achieve 100% code coverage testing, so I should found a way of testing that code.</p>

<p><!--more--></p>

Here's the code I need to test:

```python
import socketserver
from http import server


class DogFeederServer(server.BaseHTTPRequestHandler):
    def __init__(self, camera_output, call_dog, servo, *args, **kwargs):
        self.camera_output = camera_output
        self.call_dog = call_dog
        self.servo = servo
        # BaseHTTPRequestHandler calls do_GET **inside** __init__ !!!
        # So we have to call super().__init__ after setting attributes.
        super().__init__(*args, **kwargs)

    def do_GET(self):
        if self.path == "/stream.mjpg":
            self.send_response(200)
            # do some magic with HTTP Streaming
        else:
            self.send_error(404)
        self.end_headers()

    def do_POST(self):
        if self.path == "/api/call":
            if self.call_dog():
                self.send_response(200)
            else:
                self.send_response(500)
        elif self.path == "/api/treat":
            self.servo.open_and_close()
            self.send_response(200)
        else:
            self.send_error(404)
        self.end_headers()


class StreamingServer(socketserver.ThreadingMixIn, server.HTTPServer):
    allow_reuse_address = True
    daemon_threads = True

```
As you can see, the code is really simple.

The problem comes when you realise there are no easy way of calling the constructor of the server and pass the dependencies

## Passing dependencies on the constructor

Hopefully I discovered this StackOverflow post where someone has experience the same issue: <a href="https://stackoverflow.com/questions/21631799/how-can-i-pass-parameters-to-a-requesthandler">https://stackoverflow.com/questions/21631799/how-can-i-pass-parameters-to-a-requesthandler</a>

I really like the approach of the "partial" application: we pass the arguments before and once the app is created with the arguments, is passed to the server:

```python
address = ("", 8000)
handler = partial(
    DogFeederServer,
    camera_output,
    call_dog,
    servo,
)
server = StreamingServer(address, handler)
server.serve_forever()
```

Once we have the "partial" approach, we could easily provide mocks for the dependencies in the tests

## Test the server

The only way of testing the base HTTP server I found is to create some sort of "integration testing": provide mocks to the server but actually start the HTTP server. To test the whole logic, we could use `requests` library to do the HTTP calls:

```python
import socket
from functools import partial
from threading import Thread
from unittest import TestCase
from unittest.mock import MagicMock

import requests

from dogfeeder.server import DogFeederServer, StreamingServer


class ServerTest(TestCase):
    def setUp(self):
        super(ServerTest, self).setUp()
        self.get_free_port()
        self.camera_output_mock = MagicMock()
        self.call_dog_mock = MagicMock()
        self.servo_mock = MagicMock()
        address = ("", self.mock_server_port)
        handler = partial(
            DogFeederServer,
            self.camera_output_mock,
            self.call_dog_mock,
            self.servo_mock,
        )
        self.mock_server = StreamingServer(address, handler)

        # Start running mock server in a separate thread.
        # Daemon threads automatically shut down when the main process exits.
        self.mock_server_thread = Thread(target=self.mock_server.serve_forever)
        self.mock_server_thread.setDaemon(True)
        self.mock_server_thread.start()

    def test_servo_open_close(self):
        url = f"http://localhost:{self.mock_server_port}/api/treat"
        response = requests.post(url)
        self.servo_mock.open_and_close.assert_called_once()
        assert response.status_code == 200


    def test_invalid_path(self):
        url = f"http://localhost:{self.mock_server_port}/unknown"
        response = requests.post(url)
        assert response.status_code == 404
        response = requests.get(url)
        assert response.status_code == 404

    def tearDown(self):
        super(ServerTest, self).tearDown()

    def get_free_port(self):
        s = socket.socket(socket.AF_INET, type=socket.SOCK_STREAM)
        s.bind(("localhost", 0))
        __, port = s.getsockname()
        s.close()
        self.mock_server_port = port
```

The key here is to start a daemon thread (that will die when the test ends) to start the HTTP server
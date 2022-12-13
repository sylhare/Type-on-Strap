---
layout: post
title: Python cli application
description: In this article I describe how to create a python cli application using the click library and install it using pipx
author-id: "galera"
categories: [python, cli, click]
tags: [python, cli, click]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/python-cli-app/featured.jpg"
thumbnail: "assets/img/posts/python-cli-app/featured.jpg"
image: "assets/img/posts/python-cli-app/featured.jpg"
---
I have been developing an Python application to be able to be run as a CLI tool. In this article I speak about how I implement that using click library.

<p><!--more--></p>

Since I want to create a nice CLI tool, I search some libraries and found the click library: <a href="https://pypi.org/project/click/">https://pypi.org/project/click/</a>. 

This tool allows the developer to create a beatiful CLI interaction with minimal code.

Each operation you can perform with the application is called Command. In thise case we want to perform two operations: `hello` and `bye`, so we'll create two commands for that.

To define a command, we need to do two things:
- Create the source and annotate it with `@click.group()`:

```python
@click.group()
def bye_source():
    pass
```

- Create the command method and annotate it with the source created before

```python
@bye_source.command()
@user_option_required
@verbose_option
def bye(**opts):
    user_name = opts["user"]
    say_bye(user_name)


def say_bye(user_name):
    debug(f"Saying bye to {user_name}")
    info(f"Bye {user_name}")
```

I decide to separate the `bye` in two methods because this way testing the actual logic is way easier. If we have it only in the `bye`, we need to use the `click` library in testing as well.

Once we have created the commands, we must define an entrypoint which imports them and setup click library:

```python
import click
from examplecli.commands.hello import hello_source
from examplecli.commands.bye import bye_source

WELCOME_MESSAGE = """
Welcome to Example CLI!
"""

def start():
    cli = click.CommandCollection(
        sources=[hello_source, bye_source], help=WELCOME_MESSAGE)
    cli()


if __name__ == '__main__':
    start()
```

Now we can run the program directly from the console:

```bash
python3 examplecli/entrypoint.py      
Usage: entrypoint.py [OPTIONS] COMMAND [ARGS]...

  Welcome to Example CLI!

Options:
  --help  Show this message and exit.

Commands:
  bye
  hello
```

If we want to package the application, we need to do one more thing. In the `setup.py` file we must specify the console entrypoiny:

```python
#!/usr/bin/env python

from setuptools import setup, find_packages
from os import environ

with open('requirements.txt') as fp:
    install_requires = fp.read()

with open('requirements_test.txt') as fp:
    tests_require = fp.read()

setup(name='example-cli',
      version=environ.get('EXAMPLE_CLI_VERSION', '0.0.1'),
      description='Example CLI',
      author='Adrian Galera',
      author_email='',
      python_requires='>=3.6.*',
      packages=find_packages(),
      install_requires=install_requires,
      tests_require=tests_require,
      entry_points={
          'console_scripts': [
              'example-cli = examplecli.entrypoint:start',
          ]
      }
      )
```

Now we can build the app and install it using pipx:

```bash
EXAMPLE_CLI_VERSION=$VERSION python3 setup.py sdist bdist_wheel
pipx install dist/example_cli-0.0.1+local*-py3-none-any.whl --force
```

pipx will create the binary cli tool. Now we can run it as a standalone application:

```bash
example-cli hello --user test
2022-12-13T15:40:23.029Z | Hello test
```
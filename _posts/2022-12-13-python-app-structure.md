---
layout: post
title: Python app structure
description: Description of how you should organise a Python executable application to publish and generate an executable
author-id: "galera"
categories: [python, architecture]
tags: [python, architecture]
excerpt_separator: <!--more-->
feature-img: "assets/img/posts/python-app-structure/featured.jpg"
thumbnail: "assets/img/posts/python-app-structure/featured.jpg"
image: "assets/img/posts/python-app-structure/featured.jpg"
---

Recently I was starting a Python application project from scratch and I had some issues understanding the correct project structure. Here's what I have learnt.

<p><!--more--></p>

## Project structure

Here you can find the project structure of a python application implementing a cli tool named `examplecli`. You can find that app here: <a href="https://github.com/adriangalera/examplecli">https://github.com/adriangalera/examplecli</a>

```
├── Makefile
├── README.md
├── examplecli
│   ├── __init__.py
│   ├── commands
│   │   ├── __init__.py
│   │   ├── bye.py
│   │   ├── hello.py
│   │   └── options.py
│   ├── common
│   │   ├── __init__.py
│   │   └── logging.py
│   └── entrypoint.py
├── requirements.txt
├── requirements_test.txt
├── scripts
│   ├── local-build.sh
│   └── set-module-in-path.sh
├── setup.py
└── tests
    ├── __init__.py
    └── commands
        ├── __init__.py
        ├── test_bye.py
        ├── test_hello.py
        └── test_options.py
```

Let's list the most important stuff:

- Makefile: contains a series of instructions on how to perform common tasks: clean, test, lint, coverage and building locally
- examplecli: main module of the application. Contains all the application code organized as well in sub-modules.
- README.md: Readme file that contains some documentation and help materials
- requirements.txt: the pip libraries the app needs to be able to execute
- requirements_test.txt: the pip libraries to run the app tests
- scripts: useful tools for local development
- setup.py: we'll discuss in a following section about this file
- tests: folder that contains the tests for the app

## Virtual environment and dependencies

In order to have a clean environment, it is very recommended to use virtual environments. This fancy feature will isolate the dependencies needed for every application.

In order to boostrap the virtual environment, you should create it (if not created) and activate it:

```
python -m venv .venv
source .venv/bin/activate
```

Once the virtual environemnt is setup, you can install the dependencies in:

```
pip install -r requirements.txt
pip install -r requirements_test.txt
```

This will store the dependencies under the `.venv` folder

## Python modules

You should organise the python source code around the idea of modules. Modules are basically a folder with an empty file named `__init__.py` and the source code.

*it is very important to add the .py extension, otherwise it's not recognized as a module*

Different parts of the application can import the modules, e.g.:

```python
import click

from examplecli.common.logging import debug, info
from examplecli.commands.options import verbose_option, user_option_required


@click.group()
def hello_source():
    pass


@hello_source.command()
@verbose_option
@user_option_required
def hello(**opts):
    user_name = opts["user"]
    say_hello(user_name)


def say_hello(user_name):
    debug(f"Saying hello to {user_name}")
    info(f"Hello {user_name}")
```
This file is importing the methods `debug` and `info` from `logging` file in module `examplecli.common`

## setup.py

This file is only needed if we want to package the application. Without it, we are still able to run the application by calling the entrypoint:

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

However, for this application, we want to package it into a wheel file and install it with `pipx`. In order to do that, we need the setup.py file. Let's see what's inside of this file:

```python
#!/usr/bin/env python

from setuptools import setup, find_packages
from os import environ

# reads the requirements and stores them into a variable
with open('requirements.txt') as fp:
    install_requires = fp.read()

# reads the test requirements and stores them into a variable
with open('requirements_test.txt') as fp:
    tests_require = fp.read()

setup(name='example-cli', # name of the package
      version=environ.get('EXAMPLE_CLI_VERSION', '0.0.1'), # reads the variable from a environment variable
      description='Example CLI', # provides a description of the package
      author='Adrian Galera', # provides the author of the package
      author_email='',
      python_requires='>=3.6.*', # details the version compatibility.
      packages=find_packages(), # the find_packages method scans the folder for modules and sub-modules
      install_requires=install_requires,
      tests_require=tests_require,
      entry_points={
          'console_scripts': [
              'example-cli = examplecli.entrypoint:start', # required for click framework to find the starting point
          ]
      }
      )
```

In order to generate the wheel file, the user should run the following command:

```bash
EXAMPLE_CLI_VERSION=$VERSION python3 setup.py sdist bdist_wheel
```
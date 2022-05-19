---
layout: post
title: Dark Mode
tags: [Katex, Mermaid, Markdown]
---

More colors with less light. Click the **half-moon** most top-right button to turn the lights ON/OFF.
Here is a bit of everything, so you can check how the theme look, have fun! 👌


# Headers
## Level 2
### Level 3
#### Level 4
##### Level 5
###### Level 6
 
# [Headers with links](http://localhost)
## [Level 2](http://localhost)
### [Level 3](http://localhost)
#### [Level 4](http://localhost)
##### [Level 5](http://localhost) 
###### [Level 6](http://localhost)

## Code highlight
Mode specific code highlighting themes. [Kramdown](https://kramdown.gettalong.org/) which is responsible for the color highlighting may be more limited than your IDE.

```python
#!/usr/bin/env python
"""
Test file for syntax
"""
# TODO: Use dark mode
from sys import os

def foo(bar): 
    try:
        print(bar)
    except NameError:
        print("Variable bar is not defined")


class Bar(object): 
    def __init__(self):
        foo(1)
        self.octal = '\04'
        self.text = """Example \t\n"""
    
    def __exit__(self, *args):
        print('exit\u1111\xFF')
        pass
    
    @staticmethod
    def example():
        assert (1.0 and 2L) or True
        return { "example": [(1,), (r'raw', u'unicode')]}
```

## Tables

| hex | dec | oct |
| -   | -   | -   |
| 0   | 0   | 0   |
| 5   | 5   | 5   |
| A   | 10  | 12  |
| F   | 16  | 20  |
| F5  | 21  | 25  |

## KaTeX

Some KaTeX diagrams to check in dark mode:

$$
\begin{CD}
A @>a>> B \\
@VbVV @AAcA \\
C @= D
\end{CD}
$$

$$\utilde{AB}$$

## Mermaid

<div class="mermaid">
flowchart TB
    c1-->a2
    subgraph one
    a1-->a2
    end
    subgraph two
    b1-->b2
    end
    subgraph three
    c1-->c2
    end
</div>

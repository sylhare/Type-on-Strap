---
published: true
tag:
  - environment
  - tools
---

## docker for emacs, rtags
suitable for run on remote server, code, autocomplete with rtags (quicker than vscode when you need complete huge code base but sometime rtags indexing incorrect)  

Feature on emacs : magit, vterm, multiedit, sneaker, plenty of things with org mode.

rtags usage: 
+ [rtags](https://github.com/Andersbakken/rtags)
+ [how to compile database for rtags](https://sarcasm.github.io/notes/dev/compilation-database.html)

Command :
``` 
    docker pull longkl/rtags 
    env TERM=xterm-256color emacs -nw
```

FAQ:
1. How stop emacs when it hang in docker ?  
-> at host terminal : sudo pkill -SIGUSR2 emacs

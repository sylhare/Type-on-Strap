---
published: true
tag:
  - dev-log
  - tools
---

# Table of Contents

1.  [where did I start ? what problems are ?](#orgcc86f88)
2.  [Solution](#orgde7d8fb)
3.  [Note](#orge865e86)

<a id="orgcc86f88"></a>

# where did I start ? what problems are ?

-   c/c++ programmer (embedded)
-   Work with huge code base which located on a server (strong cpu, ram, disk, &#x2026;)
-   Try to navigate source code at local machine with moderm text editor, ide.
    -   Pull a part of code from remote, build. VScode, sublime text seem slow when try to indexing definition, reference
    -   Using text editor, ide to navigate huge source doesn&rsquo;t seem a good idea.
    -   don&rsquo;t have enough header for libraries. trying to get all header to local in case of huge source also isn&rsquo;t good idea. I tried that with header (7xx MB), vscode is tired with the header :( .
-   Maping source code to local from remote (NFS) also painful
-   Without support from ide, text editor for flychecking
    -   re-built when syntax errors, don&rsquo;t realize potential bugs which can be alerted from compiler


<a id="orgde7d8fb"></a>

# Solution

-   Use tool for indexing huge source code
    -   rtags, cscope, ctags, &#x2026;
-   Use flychecking which save me much time.
    -   syntax
    -   warning from compiler (init, cast, miss header, &#x2026;.). this reduce time for sanity tasks.

-   navigate code on built server instead of local machine
    -   local machine weak than built server
    -   local machine take time when try to sync code from the repotories
    -   if built environment isn&rsquo;t docker, try to navigate code in local might get incorect result

-   with emacs
    -   [rtags](<https://github.com/Andersbakken/rtags>)
    -   flychecking with rtags


<a id="orge865e86"></a>

# Note

-   a file should be independence. Otherwise flychecking might meet many problems
    -   header is independent
    -   source code is independent


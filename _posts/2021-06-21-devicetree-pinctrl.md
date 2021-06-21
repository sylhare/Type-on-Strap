---
published: true
tag:
  - embedded linux
---
## Concept, Example about device tree and pinctrl module in Linux

1. Concept  
A Device Tree (DT) is an easy-to-read hardware description file, with a JSON-likeformatting style, which is a simple tree structure where devices are represented by nodeswith their properties.  
A DT is enabled in the kernel by setting the CONFIG_OF option to Y

header:  
```
#include <linux/of.h>
#include <linux/of_device.h>
```



## Reference:  
* [GPIO API](https://lwn.net/Articles/532714/)






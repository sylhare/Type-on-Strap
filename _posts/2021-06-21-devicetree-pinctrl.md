---
published: true
tag:
  - embedded linux
---
## Concept, Example about device tree and pinctrl module in Linux

### Device Tree
1. Concept  
A Device Tree (DT) is an easy-to-read hardware description file, with a JSON-likeformatting style, which is a simple tree structure where devices are represented by nodeswith their properties.  

2. Mechanic
A DT is enabled in the kernel by setting the CONFIG_OF option to Y

header:  
```
#include <linux/of.h>
#include <linux/of_device.h>
```
some data types used in DTs:

```
nodename@reg{ // @reg is optional   
	string-property = "a string";   //string
	string-list = "red fish", "blue fish";  // list of strings 
	one-int-property = <197>; //integer 32 bits unsigned integer
	boolean-property = true;
	status = okay;
}

Alias, lables, and phandle
aliases {    
	ethernet0 = &fec;    
	gpio0 = &gpio1;    
	gpio1 = &gpio2;    
	mmc0 = &usdhc1;    
	[...]
};
gpio1: gpio@0209c000 {    
	compatible = "fsl,imx6q-gpio", "fsl,imx35-gpio";    
	[...]};
node_label: nodename@reg {    
	[...];    
     // &gpio1 is converted to a phandle so that it refers to the gpio1
	gpios = <&gpio1 7 GPIO_ACTIVE_HIGH>;
};
```
runtime representation of a DT in /proc/device-tree ( need CONFIG_PROC_DEVICETREE=y )

3. Representing and addressing devices










## Reference:  
* [GPIO API](https://lwn.net/Articles/532714/)






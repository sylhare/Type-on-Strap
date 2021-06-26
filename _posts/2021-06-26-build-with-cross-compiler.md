---
published: true
publish: true
tag: embedded linux
---
## Example, Note

To reduce time for build a program and saving environment for development later. I decided to use cross-compiler to build binary for Beaglebone Black.

### Example build DTC
1. Problem with original dtc  
currently when I try to use dtc -I fs /proc/device-tree to get all information about device tree in BBB. I met segmentation fault error, after check on the internet I saw this bug is fixed in new dtc version but BBB's OS don't support new version.  

2. Build dtc from host machine with cross-compiler  
+ Download [toolchain](https://developer.arm.com/tools-and-software/open-source-software/developer-tools/gnu-toolchain/gnu-a/downloads)  
+ Download dtc on github https://github.com/dgibson/dtc  
+ Build with following README docs and change some code, config to code can be built.  

```
command line :  
make  SETUP_PREFIX=~/working/device-driver/pre_install PREFIX=~/working/device-driver/pre_install  
NO_PYTHON=1  CC=/opt/gcc-arm-linux/bin/arm-linux-gnueabihf-gcc CFLAGS='-Wl,-rpath=/opt/gcc-arm-  
linux/arm-linux-gnueabihf/libc/ --sysroot=/opt/gcc-arm-linux/arm-linux-gnueabihf/libc/ -fPIC' install

SETUP_PREFIX and NO_PYTHON are flags which I follow README file  

PREFIX is what directory I wan't install to.  

CC is cross-compiler  
CLFAGS -rpath is library for ld processing searching share library  

--sysroot change root directory due to I use lib, bin, include of BBB(debian, arm 32) instead of host machine (x86_64)



```

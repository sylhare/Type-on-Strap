---
published: true
tag: embedded linux
---
## Theory, API, Example  

## Platform Device
hey are **a part of SoC**, they can't be removed, are **non-discoverable**, and are also called platform devices.  

**platform_driver_register**() registers and puts the driver into a list ofdrivers maintained by the kernel, so that its probe() function can be called ondemand whenever a new match occurs.  

With **platform_driver_probe**(), the kernel immediately runs the match loop,checks if there is a platform device with the matching name, and then calls thedriver's probe() if a match occurred, meaning that the device is present,[defintion here](https://elixir.bootlin.com/linux/v4.1/source/drivers/base/platform.c#L611).  

```
module_platform_driver(struct platform_driver) for platform drivers,dedicated to devices that do not sit on conventional physical buses (we just usedit before)
module_spi_driver(struct spi_driver) for SPI drivers
module_i2c_driver(struct i2c_driver) for I2C drivers
module_pci_driver(struct pci_driver) for PCI drivers
module_usb_driver(struct usb_driver) for USB drivers
module_mdio_driver(struct mdio_driver) for mdio
```

```
struct platform_device {   
	const char *name;   
    u32 id;   
    struct device dev;   
    u32 num_resources;   //number element in resource array
    struct resource *resource; //array
    }
```
### Providing data for driver
1. Device provisioning - the old and deprecated way.  
this method is to be used with the kernel version that does not support a device tree.  
**Resources**  

```
struct resource {        
	resource_size_t start;        
	resource_size_t end;        
	const char *name;        
	unsigned long flags;  
}
```

Resource type 

```
#define IORESOURCE_IO  0x00000100  /* PCI/ISA I/O ports */
#define IORESOURCE_MEM 0x00000200  /* Memory regions */
#define IORESOURCE_REG 0x00000300  /* Register offsets */
#define IORESOURCE_IRQ 0x00000400  /* IRQ line */
#define IORESOURCE_DMA 0x00000800  /* DMA channels */
#define IORESOURCE_BUS 0x00001000  /* Bus *

```

defined resource like this 
```
* Our resource array 
static struct resource needed_resources[] = {   
[0] = {        
    .start = JZ4740_UDC_BASE_ADDR,         
    .end   = JZ4740_UDC_BASE_ADDR + 0x10000 - 1,         
    .flags = IORESOURCE_MEM,         
    .name  = "mem1",   
},   
[1] = {         
	.start = JZ4740_UDC_BASE_ADDR2,
    [....],
```

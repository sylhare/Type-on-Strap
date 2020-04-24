---
layout: post
title: Blogging with title
tags: [Test, Markdown]
---


Search should be working even for complicated escape symbols
```
sed -i 's/\"hostname\"\:.*$/\"hostname\"\: \"'$IPADDR'\"\,/g' open-falcon/agent/config/cfg.json
```

## Some tech stuff blog

Because you might put code in your blog post and you want to make sure 
it will look good in here.
And that the search function is working!

### Java

java example

```java
public class Demo {
  private static final String CONSTANT = "String";
  private Object o;
  /**
   * Creates a new demo.
   * @param o The object to demonstrate.
   */
  public Demo(Object o) {
    this.o = o;
    String s = CONSTANT + "Other";
    int i = 123;
  }
  public static void main(String[] args) {
    Demo demo = new Demo();
  }
}
```

### HTML

html example

```html
HTML
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">
<!-- Sample comment -->
<HTML>
<head>
<title>IntelliJ IDEA</title>
</head>
<body>
  <h1>IntelliJ IDEA</h1>
  <p><br><b><IMG border=0 height=12 src="images/hg.gif" width=18 >
    Hello&nbsp;World! &#x00B7; &Alpha; </b><br><br>
</body>
</html>

```

### Ruby

ruby example

```ruby
require "test"
CONSTANT = 777

# Sample comment

class Module::Class
  include Testcase

  render :action => 'foo'
  def foo(parameter)
    @parameter = parameter
  end

  local_var = eval <<-"FOO";\
  printIndex "Hello world!"
  And now this is heredoc!
  printIndex "Hello world again!"
  FOO
  foo("#{$GLOBAL_TIME >> $`} is \Z sample \"string\"" * 777);
  if ($1 =~ /sample regular expression/ni) then
  begin
    puts %W(sample words), CONSTANT, :fooo;
    do_something :action => "action"
  end
  expect{counter[0]}.to_be eq 1
  1.upto(@@n) do |index| printIndex 'Hello' + index end
  \\\\\\\\\
  end
end
```

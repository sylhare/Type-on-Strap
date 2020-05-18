---
layout: post
title: Generate a table of content
tags: [toc.js, kramdown, Markdown, Customization]
author: sylhare
excerpt_separator: <!--more-->
---

# Using Kramdown GFM <!--more-->

<!-- To be placed at the beginning of the post, it is where the table of content will be generated -->
* TOC
{:toc}

## Basic Usage


You need to put this at the beginning of the page where you want the table of content to be displayed

```html
* TOC
{:toc}
```

It will then render the markdown and html titles (lines that begins with `#` or using the `<h1></h1>` tages)

# Using toc.js

Demo display of [jekyll-table-of-contents](https://github.com/ghiculescu/jekyll-table-of-contents) by ghiculescu.

<!-- To be placed at the beginning of the post, it is where the table of content will be generated -->
<div id="toc"></div>

## Customize with toc.js

[toc.js](https://github.com/ghiculescu/jekyll-table-of-contents) stands for table of content, it is a js plugin that generates automatically a table of content of a post.

### Use with this jekyll template

If you want to customize the theme it is up to you, you can add the `toc.js` file into the `asset > js` and add it into the `page.html` layout with:

```html
<script src="{{ "/assets/js/toc.js" | relative_url }}" ></script>
```
Then you can use it as it is said on the repository.

## Basic Usage

The script requires jQuery. First, reference toc.js in templates where you would like to add the table of content. Then, create an HTML element wherever you want your table of contents to appear:

```html
<div id="toc"></div>
```

Then you put your post with titles and all like:

```apiblueprint
## Title
## Mid title 1
This is text on page one
## Mid title 2
This is text for page two
### Sub title 2.a
Some more text
```

Then at the end of your post, you call the .toc() function using:

```html
<script type="text/javascript">
$(document).ready(function() {
    $('#toc').toc();
});
</script>
```

## How it would look like

![image](https://user-images.githubusercontent.com/20642750/39189661-c22099f2-47a0-11e8-826e-2ec3ef4cc4f4.png)

<script>
// toc.js 
// Copied here for the example, can be placed in assets/js for real use in your template.
// https://github.com/ghiculescu/jekyll-table-of-contents
(function($){
  $.fn.toc = function(options) {
    var defaults = {
      noBackToTopLinks: false,
      title: '<i>Jump to...</i>',
      minimumHeaders: 3,
      headers: 'h1, h2, h3, h4, h5, h6',
      listType: 'ol', // values: [ol|ul]
      showEffect: 'show', // values: [show|slideDown|fadeIn|none]
      showSpeed: 'slow', // set to 0 to deactivate effect
      classes: { list: '',
                 item: ''
               }
    },
    settings = $.extend(defaults, options);

    function fixedEncodeURIComponent (str) {
      return encodeURIComponent(str).replace(/[!'()*]/g, function(c) {
        return '%' + c.charCodeAt(0).toString(16);
      });
    }

    function createLink (header) {
      var innerText = (header.textContent === undefined) ? header.innerText : header.textContent;
      return "<a href='#" + fixedEncodeURIComponent(header.id) + "'>" + innerText + "</a>";
    }

    var headers = $(settings.headers).filter(function() {
      // get all headers with an ID
      var previousSiblingName = $(this).prev().attr( "name" );
      if (!this.id && previousSiblingName) {
        this.id = $(this).attr( "id", previousSiblingName.replace(/\./g, "-") );
      }
      return this.id;
    }), output = $(this);
    if (!headers.length || headers.length < settings.minimumHeaders || !output.length) {
      $(this).hide();
      return;
    }

    if (0 === settings.showSpeed) {
      settings.showEffect = 'none';
    }

    var render = {
      show: function() { output.hide().html(html).show(settings.showSpeed); },
      slideDown: function() { output.hide().html(html).slideDown(settings.showSpeed); },
      fadeIn: function() { output.hide().html(html).fadeIn(settings.showSpeed); },
      none: function() { output.html(html); }
    };

    var get_level = function(ele) { return parseInt(ele.nodeName.replace("H", ""), 10); };
    var highest_level = headers.map(function(_, ele) { return get_level(ele); }).get().sort()[0];
    var return_to_top = '<i class="icon-arrow-up back-to-top"> </i>';

    var level = get_level(headers[0]),
      this_level,
      html = settings.title + " <" +settings.listType + " class=\"" + settings.classes.list +"\">";
    headers.on('click', function() {
      if (!settings.noBackToTopLinks) {
        window.location.hash = this.id;
      }
    })
    .addClass('clickable-header')
    .each(function(_, header) {
      this_level = get_level(header);
      if (!settings.noBackToTopLinks && this_level === highest_level) {
        $(header).addClass('top-level-header').after(return_to_top);
      }
      if (this_level === level) // same level as before; same indenting
        html += "<li class=\"" + settings.classes.item + "\">" + createLink(header);
      else if (this_level <= level){ // higher level than before; end parent ol
        for(var i = this_level; i < level; i++) {
          html += "</li></"+settings.listType+">"
        }
        html += "<li class=\"" + settings.classes.item + "\">" + createLink(header);
      }
      else if (this_level > level) { // lower level than before; expand the previous to contain a ol
        for(i = this_level; i > level; i--) {
          html += "<" + settings.listType + " class=\"" + settings.classes.list +"\">" +
                  "<li class=\"" + settings.classes.item + "\">"
        }
        html += createLink(header);
      }
      level = this_level; // update for the next one
    });
    html += "</"+settings.listType+">";
    if (!settings.noBackToTopLinks) {
      $(document).on('click', '.back-to-top', function() {
        $(window).scrollTop(0);
        window.location.hash = '';
      });
    }

    render[settings.showEffect]();
  };
})(jQuery);
</script>

<!-- To be copied at the end of the post to render the table of content -->
<script type="text/javascript">
$(document).ready(function() {
    $('#toc').toc();
});
</script>

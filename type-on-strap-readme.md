# Type on Strap 🎨

[![Build](https://github.com/sylhare/Type-on-Strap/actions/workflows/jekyll-build.yml/badge.svg)](https://github.com/sylhare/Type-on-Strap/actions/workflows/jekyll-build.yml)
[![Gem Version](https://badge.fury.io/rb/type-on-strap.svg)](https://badge.fury.io/rb/type-on-strap)
[![Docker Pulls](https://img.shields.io/docker/pulls/sylhare/type-on-strap)](https://hub.docker.com/r/sylhare/type-on-strap)

[![Default Type on Strap blog](https://github.com/Sylhare/Type-on-Strap/blob/master/assets/img/screenshot.png?raw=true)](https://sylhare.github.io/Type-on-Strap/)

A free and open-source [Jekyll](https://jekyllrb.com) theme. 
Based on Rohan Chandra [type-theme](https://github.com/rohanchandra/type-theme) packed with extra features and easily customizable:

* Responsive design on all devices (🖥, 💻, 📱, ...)
* Portfolio 🗂, Gallery 🖼 pages for your projects
* Multi comments 💬 options  
* Tags compatibility 🏷
* Handle _Bootstrap_'ed pages: [Get Bootstrap](http://getbootstrap.com/)
* 🔎 Search feature: [Simple-Jekyll-Search](https://github.com/christian-fei/Simple-Jekyll-Search)
* Math Rendering : [KateX](https://github.com/Khan/KaTeX)
* Diagram Rendering: [Mermaid-js](https://github.com/mermaid-js/mermaid)
* 🖋 Nice fonts: [Font Awesome](https://fontawesome.com/), [Source Sans Pro](https://fonts.google.com/specimen/Source+Sans+Pro), [Pacifico](https://fonts.google.com/specimen/Pacifico?selection.family=Pacifico) 
* Seo Tags: [Jekyll-seo-tag](https://github.com/jekyll/jekyll-seo-tag)
* 🛠 Syntax Highlighting: Easily customisable [Base16](https://github.com/chriskempson/base16)
* 💡 Light and dark theme supported
* Find free of rights images on [pexels](https://www.pexels.com/)

> [Demo Site](https://sylhare.github.io/Type-on-Strap/) 

## Usage

### As a ruby gem 💎

Check out this tutorial: [Use as Ruby Gem](#use-as-ruby-gem-)

[![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#/https://github.com/sylhare/Type-On-Strap)

### As a github page 📋

1. Fork and clone the [Type on Strap repo](https://github.com/sylhare/Type-On-Strap): `git clone https://github.com/Sylhare/Type-on-Strap.git`
2. Install [Jekyll](https://jekyllrb.com/docs/installation/): `gem install jekyll`, check [#1](https://github.com/Sylhare/Type-on-Strap/issues/1) if you have a problem.
3. Install the theme's dependencies: `bundle install`
4. Customize the theme
	- GitHub Page: [update `_config.yml`](#site-configuration)
5. Run the Jekyll server: `bundle exec jekyll serve`

## Structure

Here are the main files of the template

```shell
./Type-on-Strap
├── _includes	               # Theme includes
├── _layouts                   # Theme layouts (see below for details)
├── _portfolio	               # Collection of articles for the portfolio page
├── _posts                     # Blog posts
├── _sass                      # Sass partials (compiled into css at runtime)
├── assets
|  ├── js	               # JS compiled for distribution + raw sources
|  ├── css                     # CSS compiled for distribution
|  ├── fonts		       # Font-Awesome, and other fonts
|  └── img		       # Images used for the template
├── pages
|   ├── 404.md		       # To be displayed when url is wrong
|   ├── about.md               # About example page
|   ├── gallery.md             # Gallery page for your photos
|   ├── portfolio.md	       # Portfolio page for your projects
|   ├── search.md	       # Search page
|   └── tags.md                # The tag page
├── _config.yml                # sample configuration
├── _data.yml
|  ├── authors.yml             # Update the post authors configurations 
|  ├── language.yml            # Localization configuration
|  ├── biblio.yml              # To create a reference bibliography
|  ├── social.yml              # Social configurations to share posts (RSS, shares, ...)
|  └── icons.yml               # Footer icons (Twitter, Github, Stackoverflow, ...)
└── index.html                 # sample home page (blog page paginated)
```
	
## Configure Type on Strap 🛠

Open `_config.yml` in a text editor to change most of the blog's settings.

If a variable in this document is marked as "optional", disable the feature by removing all text from the variable. 

### Site configuration

#### Base url

Configure Jekyll as your own blog or with a "baseurl" in `_config.yml`:

Jekyll website *without* a "baseurl" (such as a **GitHub Pages website** with your username as the repository name):

```yml
baseurl: ""
url: "https://username.github.io"
```

Jekyll website *with* "baseurl" (like the Type on Strap [demo](https://sylhare.github.io/Type-on-Strap/) page):

```yml
baseurl: "/sub-directory"
url: "https://username.github.io"
```

#### Jekyll blog configuration 

And here is the basic information you will need in your `_config.yml` for it to work properly:

```yaml
# BLOG CONFIGURATION
post_navigation: true
paginate: 10
paginate_path: "blog/page:num"
plugins: [jekyll-paginate, jekyll-seo-tag, jekyll-feed]
```

To configure the blog part and default plugins. Those plugins are validated by GitHub page.

#### Meta and Branding

_Meta variables_ hold basic information about your Jekyll site, which will be used throughout the site 
and as _meta properties_ that are used for search engines, browsers, and the site's RSS feed.

Change these variables in `_config.yml`:

```yml
title: My Jekyll Blog                 # Name of website
avatar: assets/img/avatar.png         # Path of avatar image, to be displayed in the theme's header
description: My blog posts            # Short description, primarily used by search engines
favicon: assets/favicon.ico           # Icon displayed in the tab
color_theme: auto                     # color theme auto, dark or light
```

You can also customize the seo tags default option following the jekyll-seo-tag plugin [documentation](http://jekyll.github.io/jekyll-seo-tag/advanced-usage/).
The color theme can be set to dark or light (customize it in _variables.scss_). 
Using _auto_ you'll have a tiny icon in the navbar allowing the use to manually switch from dark to light theme.

### Theme customization 🎨

#### Footer and Header text

Customize your site header/footer with these variables in `_config.yml`:

```yml
header_text: Welcome to my Jekyll blog
footer_text: Copyright 2017
```

If you don't want anything, replace the value by `" "`.

#### Header's image

The header's image (tested with 2480x1280) can be set as one image with `header_feature_image`
but can also be responsive:

```yml
header_feature_image: assets/img/header/my-header-image.png
header_feature_image_responsive: true
```

By setting `header_feature_image_responsive` to true, it will look for images 
with suffix `-small` (620x320) and `-medium` (1240x640) to display on smaller screen.

#### Localisation string

Localization string is a way to quickly change the template language for text like *Next Post* or *Follow on*, ...
You can find all the properties in `_data/language.yml`.

By default, it is in English, but you can easily add your own language.

### Google Analytics

To enable Google Analytics (GA4), add your [Measurement ID](https://support.google.com/analytics/answer/12270356?hl=en&sjid=1593376271608310401-NA) 
to `_config.yml` like so:

```yml
google_analytics: G-XXXXXXXXXX
```

It will use the [Google Tag Manager](https://support.google.com/analytics/answer/10220869?hl=en&ref_topic=9355633&sjid=1593376271608310401-NA)

### Comments 💬

#### Disqus

If you have a [Disqus](https://disqus.com/) account, you can show a comments section below each post.

To enable Disqus comments, add your [Disqus shortname](https://help.disqus.com/customer/portal/articles/466208) 
to your project's `_config.yml` file:

```yml
comments:
  disqus_shortname: my_disqus_shortname
```

#### Cusdis

[Cusdis](https://cusdis.com/) is an open-source alternative to Disqus.
You can read more about it in the [documentation](https://cusdis.com/doc#/)

To enable it, set your Cusdis name in `_config.yml`:

```yaml
comments:
  cusdis_app_id: my_data-app-id                                     
```

#### Utterances

[Utterances](https://utteranc.es) is another open source alternative linked to one's GitHub account.
It stores the comments as GitHub issues on a repository for each page.

Install the utterance [app](https://github.com/apps/utterances) to your repo.
After installing, add your info in the `_config.yml`:

```yaml
comments:
  utterances:              # Enable by filling below information. For more info, go to https://utteranc.es
    repo:                  # your public comments repository (e.g. owner/repo)
    issue-term:            # Issue term (e.g. "comment" consider issues with this word in the title as comments)
    theme:                 # OPTIONAL: Take the `color_theme` by default, or set a custom one like github-dark-orange
    label:                 # OPTIONAL: Adds an issue label in the issue
```

### Math typesetting with KateX

When KateX is set in `_config.yml`:

```yml
katex: true # to enable it
```

You can then wrap math expressions with `$$` signs in your posts and make sure you have set the `katex` variable 
in `_config.yml` to `true` for math typesetting.

For inline math typesetting, type your math expression on the *same line* as your content. For example:

```latex
Type math within a sentence $$2x^2 + x + c$$ to display inline
```

For display math typesetting, type your math expression on a *new line*. For example:

```latex
$$
  \bar{y} = {1 \over n} \sum_{i = 1}^{n}y_i
$$
```

You can find a cheat sheet of the compatible LaTex symbols [online](https://artofproblemsolving.com/wiki/index.php/LaTeX:Symbols).

### Diagrams with Mermaid

Enable the [mermaid-js](https://github.com/mermaid-js/mermaid) diagram rendering by setting mermaid to true in the `_config.yml`.
This will load and init the [mermaid.min.js](https://mermaid-js.github.io/mermaid/getting-started/n00b-gettingStarted.html#4-calling-mermaid-from-a-relative-link).

```yml
mermaid: default # Enable mermaid-js for diagrams, use theme: base, forest, dark, default, neutral
```

Find all the help you need on the official [mermaid documentation](https://mermaid-js.github.io/mermaid/).
You can create with ease diagrams. Add your mermaid script inside two mermaid divs (default Kramdown does not yet support mermaid).
With the `class="mermaid"` inside the `<div>`:

```html
<div class="mermaid">
sequenceDiagram
    Alice->>John: Hello John, how are you?
    John-->>Alice: Great!
</div>
```

### Social icons

In `_data/social.yml` you can customize the social icons that will be displayed in the post to share your post.
You can also enable RSS.
The site icons come from [Font Awesome](https://fontawesome.com/).

In `_data/icons.yml` you can set the footer icon that will appear at the bottom of the page.
They will redirect the user on your profile on to other platforms like Twitter, GitHub and so many more!

### Cookie consent

You can add a cookie consent with a disclaimer if you use Google Analytics while respecting the [GDPR](https://en.wikipedia.org/wiki/General_Data_Protection_Regulation).
Set to true, there will be a banner at the bottom of the page with the disclaimer, and an _approve_ button.
Once the user clicks on "Approve" the cookies will be created for Google Analytics.

#### Share in article

The share icons are the one at the bottom of the blog page if enabled.
They will on click redirect you to the logo's platform to share the article.

#### Footer

Display icons in the footer. 
All icon variables should be your username enclosed in quotes (e.g. "username") in `_data/icons.yml`.

You can update the RSS settings in `_data/social` to change the default feed path (generated by [jekyll-feel](https://github.com/jekyll/jekyll-feed)).
To enable the share icons at the bottom of each article set to true the one you'd like under `share` in the `_data/social.yml` file.

### Personalize your Blog Posts 📝

When writing a post, be sure in jekyll to:
 - Put it in the `_posts` folder
 - Name it with the date first like `2019-08-21-This-is-my-blog-post.md`

Please refer to the [Jekyll docs for writing posts](https://jekyllrb.com/docs/posts/). 

#### Layout: Post

These are the basic features you can use with the `post` layout, in the comment the `Opt` means that
it is optional.

```yml

---
layout: post
title: Hello World                                # Title of the page
hide_title: true                                  # [Opt] Hide the title when displaying the post, but shown in lists of posts
feature-img: "assets/img/sample.png"              # [Opt] Add a feature-image to the post
thumbnail: "assets/img/thumbnails/sample-th.png"  # [Opt] Add a thumbnail image on blog view
color: rgb(80,140,22)                             # [Opt] Add the specified colour as feature image, and change link colors in post
position: 1                                       # [Opt] Set position on the menu navigation bar
tags: [sample, markdown, html]                    # [Opt] Add tags to the page
---
```

With `thumbnail`, you can add a smaller image than the `feature-img`. 
If you don't have a thumbnail, you can still use the same image as the feature one. Or use the gulp task to create it.

If you don't use a feature image, but `color`, the transparent background is set comes from `lineart.png`. 
You can edit it in the config file (`_config.yml > color_image`). If you want another one, put it in `/assets/img` as well. 

For position, if not set on all pages, it will be by alphabetical order without `position` then by `position` order.
If two pages have the same position number, the order is decided by alphabetical order on the page title.

There's also `bootstrap: true` which is not mandatory and only useful if you want to add HTML content in your page that
requires [bootstrap](http://getbootstrap.com/).
It will respect the page and theme layout, mind the padding on the sides.

#### Post excerpt

The [excerpt](https://jekyllrb.com/docs/posts/#post-excerpts) are the first lines of an article that is displayed on the blog page. 
The length of the excerpt has a default of around `250` characters or can be manually set in the post using:

in `conflig.yml`:

```yml
excerpt: true
```

Then in your post, add the `excerpt separator`:

```yml

---
layout: post
title: Sample Page
excerpt_separator: <!--more-->
---

some text in the excerpt
<!--more-->
... rest of the text not shown in the excerpt ...
```

The html is stripped out of the excerpt, so it only displays text.

#### Image aligner

To easily add align images side by side in your article using the `aligner.html` include:

```ruby
{% include aligner.html images="path/to/img1.png,path/to/img2.png,path/to/img3.png" column=3 %}
```

Use it in any markdown file. There are two fields in the _include_ you need to look into:
  - _images_: Takes a string separated with `,` of all the images' path. 
    - It by default look into `assets/img/` so give the path from there.
  - _column_: (OPTIONAL) Set the number of column you want your imaged displayed in.
    - default is 2 columns
    - `column=3` set 3 columns
    - `column="auto"` makes as many columns as images

#### Code highlight

Like all CSS variables in the theme, you can edit the color of the code highlight in `_sass > base > _variables.scss`.
The code highlighting works with [base16](https://github.com/chriskempson/base16-html-previews/tree/master/css) you can find existing example 
of your favourite highlight color scheme on this format.

## Feature pages and layouts 

All feature pages besides the "home" one are stored in the `page` folder, 
they will appear in the navigation bar unless you set `Hide: true` in the front matter. 

Here are the documentation for the other feature pages that can be added through `_config.yml`. 

Non-standard features are documented below.

### Layout: Default

This layout includes the head, navigation bar and footer around your content. 
Unless you are making a custom layout you won't need it.

### Layout: Home 🏡

This page is used as the home page of the template (in the `index.html`). It displays the list of articles in `_posts`.
You can use this layout in another page (adding a title to it will make it appear in the navigation bar).

The recommended width and height for the home picture is width:`2484px;` and height:`1280px` 
which are the dimensions of the actual picture for it to be rolling down as you scroll the page.

If your posts are not displaying ensure that you have added the line `paginate: 5` to `_config.yml`.

### Layout: Page 📄

The page layout has a bit more features explained here.

```yml

---
layout: page
title: "About" 
subtitle: "This is a subtitle"   
feature-img: "assets/img/sample.png" 
permalink: /about/                   # Set a permalink your your page
hide: true                           # Prevent the page title to appear in the navbar
icon: "fa-search"                    # Will Display only the fontawesome icon (here: fa-search) and not the title
tags: [sample, markdown, html]
---
```

The hide only hides your page from the navigation bar, it is, however, still generated and can be accessed through its link. 

### Feature: Portfolio 🗂

Portfolio is a feature page that will take all the markdown/html files in the `_portfolio` folder to create a 3-columns image portfolio matrix.

To use the portfolio, simply create a `portfolio.md` with this information inside:

```yml

--- 
layout: page
title : Portfolio 
---

{% include default/portfolio.html %}
```

#### Portfolio posts

You can format the portfolio posts in the `_portfolio` folder using the `post layout`. 
Here is a little explanation on some of the possible features you can use.

If you decide to use a date, please be sure to use one that can be parsed such as `yyyy-mm-dd`. 
You can see more format examples in the demo posts that are available for the theme:

```yml

---
layout: post
title: Circus				       # Title of the portfolio post
feature-img: "assets/img/portfolio/cake.png"   # Will display the image in the post
img: "assets/img/portfolio/cake.png"           # Will display the image in the portfolio page
date: 2019-07-25		 	       # Not mandatory, however needs to be in date format to display the date
---
```

#### Portfolio in gem

Make sure your `_config.yml` contains the following if you are using the theme as a gem:

```yml
# PORTFOLIO
collections:
  portfolio:
    output: true
    permalink: /:collection/:name
```    

This creates the collection for Jekyll, so it can find and display your portfolio posts.

### Feature: Gallery 🖼

You can create a gallery using [Masonry JS](https://masonry.desandro.com/) which will placing the pictures at the optimal position 
based on available vertical space. 
You need to specify the `gallery_path` which will be used to find the pictures to render. 
It will take all the pictures under that directory. Then use the `include` to add it in your page. 

```yml

---
layout: page
title: Gallery
gallery: "assets/img/pexels"
---

{% include default/gallery.html gallery_path=page.gallery %}
```

### Feature: Search 🔍

The search feature is based on [Simple-Jekyll-search](https://github.com/christian-fei/Simple-Jekyll-Search) 
there is a `search.liquid` file that will create a list of all the site posts, pages and portfolios. 
Then there's a script displaying the formatted results in the _search page_.

To exclude contents from the search add the `exclude: true` option in the markdown header. 
By default, all posts, pages, and collections are available in the search.
Hide the search page from the navigation bar with the `hide: true` option. 
You can remove the icon by removing `icon`:

```yml

---
layout: search
title: Search
icon: "search"
---
```

### Feature: Tags 🏷

Tags should be placed between `[]` in your post metadata. Separate each tag with a comma. 
Tags are recommended for posts and portfolio items.

For example:

```yml

---
layout: post
title: Markdown and HTML
tags: [sample, markdown, html]
---
```

> Tags are case-sensitive `Tag_nAme` ≠ `tag_name`

All the tags will be listed on the "tags" page with a link toward the pages or posts.
The Tag page can be hidden with the `hide` option. You can remove the icon by removing `icon` (like for the search page).

## Advanced

### Liquid tags

Jekyll works with [liquid](https://shopify.github.io/liquid/) tags usually represented by:

```
{{ liquid.tag | filter }}
```

These are useful to render your jekyll files. 
You can learn more about them on [shopify's doc](https://help.shopify.com/themes/liquid/basics)

### Gulp toolbox

#### Requirements

Before you need to have *node* and `npm` installed:

- Windows: https://nodejs.org/
- Ubuntu/Debian: `apt-get install nodejs npm libgl1 libxi6`
- Fedora (dnf) / RHEL/CentOS (yum): `dnf install node npm libglvnd-glx libXi`

Then you need to install [`gulp-cli`](https://gulpjs.com/) and its dependencies:

```bash
cd assets/
sudo npm install gulp-cli -g
npm install
```

#### Minimizing and optimizing: css, js and images

You can run the default task that will compress the js, css and images and create the thumbnails for the supported image
formats:

```bash
cd assets/
gulp default
gulp thumbnails-all # to create all of the images thumbnails
gulp thumbnails     # to create thumbnails for the feature-img/ only
# tip: run a git status to see the changes
git status
```

You can find more about the gulp tasks in the [gulpfile.js](assets/gulpfile.js).

#### Create a post

To create a `.md` file in the *_posts/* section with the jekyll format of today's date.
Use this command with the title you'd like to create the very basic post.

```bash
gulp post -n 'title of the post'
```

A file will be created following the format `yyyy-mm-dd-title-of-the-post.md` with default post attributes inside.
Nothing will happen if the file exists already.

### Use as Ruby Gem 💎

You can use Type-on-strap as a [gem](https://rubygems.org/gems/type-on-strap). 

Using the [Ruby Gem Method](https://sylhare.github.io/2021/03/25/Run-type-on-strap-jekyll-theme-locally.html).
Add this line to your Jekyll site's Gemfile (or create one):

```ruby
gem "type-on-strap"
```

Add this line to your Jekyll site's `_config.yml` file:

```yml
theme: type-on-strap
```

Then run Bundler to install the theme gem and dependencies:

```bash
bundle install
```

Then you can start adding content like:
  - Add a `index.html` file
  - Add the feature page you want. (ex: as it is already in `pages`)
  - Add posts in `_posts` and `_portfolio` to be displayed

### Remote Theme

Now you can use any theme gem with GitHub pages with [29/11/2017 GitHub Pages Broadcast](https://github.com/blog/2464-use-any-theme-with-github-pages).
For that remove all `theme:` attributes from `_config.yml` and add instead:

```yml
remote_theme: sylhare/Type-on-Strap 
```

## License

This theme is licensed under the [MIT License (MIT)](/LICENSE)

- Pictures from [Pexels](https://www.pexels.com/) are under Creative Commons Zero (CC0) license
- Fonts are licensed under the [SIL Open Font License (OFL)](https://scripts.sil.org/cms/scripts/page.php?site_id=nrsi&id=OFL)

# Type on Strap

![Default Type on Strap blog](https://raw.githubusercontent.com/Sylhare/Type-on-Strap/dev/screenshot.png)

A free and open-source [Jekyll](https://jekyllrb.com) theme. Based on Rohan Chandra [type-theme](https://github.com/rohanchandra/type-theme) with a few new features:

- Responsive design
- Include bootstrap and Jekyll search
- Portfolio, tags, search page layout
- Automatic generation of seo tag (for search engines)
- Free of rights images from [pexels](https://www.pexels.com/)

## Table of Contents

1. [Usage](https://github.com/Sylhare/Type-on-Strap#Usage)
2. [Struture](https://github.com/Sylhare/Type-on-Strap#structure)
3. [Configure Type on Strap](https://github.com/Sylhare/Type-on-Strap#configure-type-on-strap)
4. [Layout](https://github.com/Sylhare/Type-on-Strap#layout)
5. [Feature pages](https://github.com/Sylhare/Type-on-Strap#feature-pages)
6. [License](https://github.com/Sylhare/Type-on-Strap#license)

## Usage

1. Fork and clone the [Type on Strap repo](https://github.com/sylhare/Type-On-Strap): `git clone https://github.com/Sylhare/Type-on-Strap.git`
2. Install [Jekyll](https://jekyllrb.com/docs/installation/): `gem install jekyll` check [#1](https://github.com/Sylhare/Type-on-Strap/issues/1) if you have a problem.
3. Install the theme's dependencies: `bundle install`
4. Customize the theme (see below)
5. Run the Jekyll server: `jekyll serve`

## Structure

Here are the main files of the template

```bash
jekyll-theme-basically-basic
├── _portofolio	               # collection of article to be populated in the portfolio page
├── _includes	               # theme includes
├── _layouts                   # theme layouts (see below for details)
├── _sass                      # Sass partials 
├── assets
|  ├── js	               # theme javascript, Katex, jquery, bootstrap, jekyll search, 
|  ├── css                     # isolated Bootstrap, font-awesome, katex and main css
|  ├── fonts		       # Font-Awesome, Glyphicon, and other fonts
|  └── img		       # Images used for the template
├── pages
|   ├── 404.md		       # To be displayed when url is wrong
|   ├── about.md               # About example page
|   ├── portfolio.html	       # Portfolio bootstrapped page
|   ├── search.html	       # Search page
|   └── search.json            # Specify the search target (page, post, collection)
├── _config.yml                # sample configuration
└── index.html                 # sample home page (blog page paginated)
```
	
## Configure Type on Strap

Open `_config.yml` in a text editor to change most of the blog's settings.

If a variable in this document is marked as "optional", disable the feature by removing all text from the variable. 


### Site configuration
Configure Jekyll as your own blog or with a subpath in in `_config.yml`:

Jekyll website *without* a subpath (such as a GitHub Pages website for a given username):

```yml
  baseurl: ""
  url: "https://username.github.io"
```

Jekyll website *with* subpath (like the Type Theme demo page):

```yml
  baseurl: "/sub-directory"
  url: "https://username.github.io/"
```

Please configure this  before using the theme.

### Meta and Branding

Meta variables hold basic information about your Jekyll site which will be used throughout the site and as meta properties for search engines, browsers, and the site's RSS feed.

Change these variables in `_config.yml`:

```yml
  theme_settings:
    title: My Jekyll Blog                 # Name of website
    avatar: assets/img/triangular.svg     # Path of avatar image, to be displayed in the theme's header
    gravatar: f98....6bfc                 # MD5 hash of your email address
    description: My blog posts            # Short description, primarily used by search engines
```

### Customizing text

#### Footer and Header's text

Customize your site header/footer with these variables in `_config.yml`:

```yml
  theme_settings:
    header_text: Welcome to my Jekyll blog
    header_text_feature_image: assets/img/sample3.png
    footer_text: Copyright 2017
```

#### Localisation string

Change localization string variables in `_config.yml`.

English text used in the theme has been grouped  so you can quickly translate the theme or change labels to suit your needs.

```yml
  theme_settings:
     str_follow_on: "Follow on"
     str_rss_follow: "Follow RSS feed"
     str_email: "Email"
     str_next_post: "Next post"
     str_previous_post: "Previous post"
     str_next_page: "Next"
     str_previous_page: "Prev"
     str_continue_reading: "Continue reading"
     str_javascript_required_disqus: "Please enable JavaScript to view comments."
```


### Other features

### Footer's icons

Display the site's icon from [Font Awesome](https://fortawesome.github.io/Font-Awesome/) in the footer. All icon variables should be your username enclosed in quotes (e.g. "username") in `_config.yml`, except for the following variables:

```yml
  theme_settings:
     rss: true
     email_address: type@example.com
     linkedin: ttps://www.linkedin.com/in/FirstLast
     stack_exchange: https://stackoverflow.com/users/0000/first-last
```

### Comments (via Disqus)

Optionally, if you have a [Disqus](https://disqus.com/) account, you can show a 
comments section below each post.

To enable Disqus comments, add your [Disqus shortname](https://help.disqus.com/customer/portal/articles/466208) to your project's `_config.yml` file:

```yml
  theme_settings:
     disqus_shortname: my_disqus_shortname
```

### Google Analytics

To enable Google Analytics, add your [tracking ID](https://support.google.com/analytics/answer/1032385) 
to `_config.yml` like so:

```yml
  theme_settings:
     google_analytics: UA-NNNNNNNN-N
```

### Math typesetting

When KateX is set in `_config.yml`:

```yml
  theme_settings:
     katex: true # to Enable it
```

You can then wrap math expressions with `$$` signs in your posts and make sure you have set the `katex` variable in `_config.yml` to `true` for math typesetting.

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

## Layout
Please refer to the [Jekyll docs for writing posts](https://jekyllrb.com/docs/posts/). Non-standard features are documented below.

### Layout: Post

This are the basic features you can use with the  `post` layout.

```yml
---
layout: post

# Title of the page
title: Hello World

# A subtitle can be displayed below your title
subtitle: "This is a subtitle"

# Add a feature-image to the post
feature-img: "assets/img/sample.png"

tags: [sample, markdown, html]
---
```

### Layout: Page

The page layout have a bit more features explained here.

```yml
---
layout: page
title: "About" 
subtitle: "This is a subtitle"   

# Add a feature-image to the post
feature-img: "assets/img/sample.png"

tags: [sample, markdown, html]

# Set a permalink your your page
permalink: /about.html  

# to prevent the page from showing up in the header's navigation bar (visitors can still visit the URL through other means).
hide: true  
---
```

### Layout: Bootstrap

This is the page layout modified to have bootstrap activated to format your content accordingly with the theme.

```yml
--- 
layout: bootstrap
---
```

### Layout: Default

This layout includes the head, navigation bar and footer around your content.

## Feature pages

All feature pages are stored in the `page` folder, they will appear in the navigation bar unless you set `Hide: true` in the front matter.

### Portfolio

Portfolio is a feature bootstrapped page that will take all the markdown/html files in the `_portfolio` folder to create a 3x3 image portfolio matrix.

### Search

The search feature is based on [Simple-Jekyll-search](https://github.com/christian-fei/Simple-Jekyll-Search) there is a `search.json` file that will create a list of all of the site posts, pages and portfolios. 

Then there's a `search.js` displaying the formated results entered in the `search.html` page. 

### Tags

Post tags should be placed between `[]` in your post metadata. Seperate each tag with a comma.

For example:

```yml
---
layout: post
title: Markdown and HTML
tags: [sample, markdown, html]
---
```

All the tags will be listed in `tags.html` with a link toward the pages or posts.


## License

[The MIT License (MIT)](https://github.com/rohanchandra/type-theme/blob/master/LICENSE)

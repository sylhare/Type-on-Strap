#!/bin/bash
JEKYLL_ENV=production bundle exec jekyll build
# Add generated robots.txt:

# Minify assets
cd assets || exit
gulp default

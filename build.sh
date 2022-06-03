#!/bin/bash
JEKYLL_ENV=production bundle exec jekyll build --incremental

# Minify assets
cd assets || exit
gulp default

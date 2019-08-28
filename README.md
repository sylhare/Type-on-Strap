# SAIL Blog

## For post authors 
To write a new post - 
1. install Jekyll
2. Run bundle install in the root dir 
3. Try running serve_jekyll, and go to the link printed in console in browser
4. If all looks good, go to \_posts and follow examples of prior posts in writing
a new one in the right format. To convert from Google doc to markdown, you can use [this add-on](https://gsuite.google.com/marketplace/app/docs_to_markdown/700168918607), download as docx and use pandoc, or do it manually.
5. Once you've written a full draft, send it on to whoever from SAIL blog you are in contact with (one of the editors listed here http://ai.stanford.edu/blog/about/) and it'll go through some editing, and then be released.

## For authors 

The way things work is that we have a 'source' branch with all the markdown and jekyll files, and the master branch has the compiled HTML. This master branch is closed to /afs/.cs/group/ai/www/blog/ and is how we update the site's contents. 

So to sync with online version, run in terminal 
1. bundle exec jekyll clean
2. export JEKYLL_ENV=production
3. jekyll build, 
4. 'octopress deploy'
5. go to /afs/.cs/group/ai/www/blog/ and pull latest from master

If you don't have octopress, install it with `gem install octopress`

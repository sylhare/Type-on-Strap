# SAIL Blog

## For post authors 
To write a new post - 
1. Get a draft of your post in Google doc, and email the lead editors ([Andrey Kurenkov and Michelle Lee](http://ai.stanford.edu/blog/about/)) to get an editor assigned to the draft.
2. Once you have a draft that is finalized, you need to create a pull request with markdown and images of your post. First, fork the repo, clone the fork and pull source branch
3. Install Jekyll
4. Run bundle install in the root dir 
5. Try running serve_jekyll, and go to the link printed in console in browser
6. If all looks good, go to \_posts and follow examples of prior posts in writing
a new one in the right format. To convert from Google doc to markdown, you can use pandoc [this add-on](https://gsuite.google.com/marketplace/app/docs_to_markdown/700168918607), download as docx and use pandoc, or do it manually.
7. Once it all looks good, submit a pull request and email your editor to let them know, and it'll be posted on the site fairly promptly.

## For editors 

The way things work is that we have a 'source' branch with all the markdown and jekyll files, and the master branch has the compiled HTML. This master branch is closed to /afs/.cs/group/ai/www/blog/ and is how we update the site's contents. 

So to sync with online version, run in terminal 
1. bundle exec jekyll clean
2. export JEKYLL_ENV=production
3. jekyll build, 
4. 'octopress deploy'
5. go to /afs/.cs/group/ai/www/blog/ and pull latest from master

If you don't have octopress, install it with `gem install octopress`

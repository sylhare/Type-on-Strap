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

### Formatting Tips

For embeddings figures, we generally recommend using {% figure %} wrapping and our supported CSS classes. See examples below: 

Single image that is 3/4th of the width of text:

    {% figure %}
    <img class="postimage_75" src="{{ site.baseurl }}/assets/img/posts/2019-10-28-learning-from-partners/image13.png"/>
    {% endfigure %}

Two images side by side:

    {% figure %}
    <img class="postimagehalf" src="{{ site.baseurl }}/assets/img/posts/2019-10-28-learning-from-partners/image16.png"/>
    <img class="postimagehalf" src="{{ site.baseurl }}/assets/img/posts/2019-10-28-learning-from-partners/image1.png"/>
    {% endfigure %}

Three videos side by side, with a caption:

    {% figure %}
    <video autoplay loop muted playsinline class="postimagethird">
      <source src="{{ site.baseurl }}/assets/img/posts/2019-11-08-roboturk/robonet.mp4" type="video/mp4">
    </video>
    <video autoplay loop muted playsinline class="postimagethird">
      <source src="{{ site.baseurl }}/assets/img/posts/2019-11-08-roboturk/collecting_2.mp4" type="video/mp4">
    </video>
    <video autoplay loop muted playsinline class="postimagethird">
      <source src="{{ site.baseurl }}/assets/img/posts/2019-11-08-roboturk/collecting_3_crop.mp4" type="video/mp4">
    </video>
    <figcaption>
    Prior data collection methodologies include autonomous data collection (left), human supervision with web interfaces (middle), and human teleoperation with motion interfaces (right).
    </figcaption>
    {% endfigure %}
   
Two rows of three images side by side, with a sub-caption for each image and a caption for overall figure:

    {% figure %}
    <figure class="postfigurethird" >
      <img src="{{ site.baseurl }}/assets/img/posts/2019-11-26-robonet/align_tshirt.gif" class="postimage_unpadded"/>
      <figcaption>
      Kuka can align shirts next to the others
      </figcaption>
    </figure>
    <figure class="postfigurethird" >
      <img src="{{ site.baseurl }}/assets/img/posts/2019-11-26-robonet/color_stripe_front.gif" class="postimage_unpadded"/>
      <figcaption>
      Baxter can sweep the table with cloth
      </figcaption>
    </figure>
    <figure class="postfigurethird">
      <img src="{{ site.baseurl }}/assets/img/posts/2019-11-26-robonet/marker_marker.gif" class="postimage_unpadded"/>
      <figcaption>
      Franka can grasp and reposition the markers
      </figcaption>
    </figure>
    <figure class="postfigurethird" >
      <img src="{{ site.baseurl }}/assets/img/posts/2019-11-26-robonet/move_plate.gif" class="postimage_unpadded"/>
      <figcaption>
      Kuka can move the plate to the edge of the table
      </figcaption>
    </figure>
    <figure class="postfigurethird" >
      <img src="{{ site.baseurl }}/assets/img/posts/2019-11-26-robonet/socks_right.gif" class="postimage_unpadded"/>
      <figcaption>
      Baxter can pick up and reposition socks 
      </figcaption>
    </figure>
    <figure class="postfigurethird">
      <img src="{{ site.baseurl }}/assets/img/posts/2019-11-26-robonet/towel_stack.gif" class="postimage_unpadded" />
      <figcaption>
      Franka can stack the towel on the pile
      </figcaption>
    </figure>
    <figcaption>
    Here you can see examples of visual foresight fine-tuned to perform basic control tasks in three entirely different environments. For the experiments, the target robot and environment was subtracted from RoboNet during pre-training. Fine-tuning was accomplished with data collected in one afternoon.
    </figcaption>
    {% endfigure %}


The full list of CSS classes is this:

    .postimage {
      width: 100%;
      height: auto;
      margin: auto;
    }

    .postimage_unpadded {
      width: 100%;
      height: auto;
      margin: auto;
      padding: 0;
    }

    .postimagesmall {
      display:block;
      margin-left:auto;
      margin-right:auto;
      height: auto;
    }

    .postimageactual {
      display:block;
      margin-left:auto;
      margin-right:auto;
      height: auto;
    }

    .postimagesmaller {
      width: 75%;
      display:block;
      margin-left:auto;
      margin-right:auto;
      height: auto;
    }

    .postimage_75 {
      width: 75%;
      display:block;
      height: auto;
      margin-left:auto;
      margin-right:auto;
    }

    .postimage_50 {
      width: 50%;
      display:block;
      height: auto;
      margin-left:auto;
      margin-right:auto;
    }

    .postimage_25 {
      width: 25%;
      display:block;
      height: auto;
      margin-left:auto;
      margin-right:auto;
    }

    .postimagehalf {
      width: 49%;
      height: auto;
      margin-left:auto;
      margin-right:auto;
      padding: 1%;
    }

    .postimagethird {
      width: 30%;
      height: auto;
      margin-left:auto;
      margin-right:auto;
      padding: 1%;
    }

    .postimagequarter {
      width: 24%;
      height: auto;
      margin-left:auto;
      margin-right:auto;
    }
    
    .postfigurehalf {
      width: 49%;
      height: auto;
      margin-left:auto;
      margin-right:auto;
      padding: 1%;
      display: inline-block;
    }

    .postfigurethird {
      width: 30%;
      height: auto;
      margin-left:auto;
      margin-right:auto;
      padding: 1%;
      display: inline-block;
    }

    .postfigurequarter {
      width: 24%;
      height: auto;
      margin-left:auto;
      margin-right:auto;
      display: inline-block;
    }

  
### Citations

We support [bigfoot](http://www.bigfootjs.com/) pop-up citations and recommend using them. In text, use [^<name>] as follows:

    - **Autonomous Data Collection**: Many data collection mechanisms and algorithms such as Self-Supervised Learning[^SSL][^robonet] 

Then at bottom of file:

    [^robonet]: Dasari, S., Ebert, F., Tian, S., Nair, S., Bucher, B., Schmeckpeper, K., ... & Finn, C. (2019). RoboNet: Large-Scale Multi-Robot Learning. arXiv preprint arXiv:1910.11215.
    [^SSL]: Levine, S., Pastor, P., Krizhevsky, A., & Quillen, D. (2016, October). Learning hand-eye coordination for robotic grasping with large-scale data collection. In International Symposium on Experimental Robotics (pp. 173-184). Springer, Cham.

## For editors 

The way things work is that we have a 'source' branch with all the markdown and jekyll files, and the master branch has the compiled HTML. This master branch is cloned to /afs/.cs/group/ai/www/blog/ and is how we update the site's contents. 

So to sync with online version, run in terminal 
1. bundle exec jekyll clean
2. export JEKYLL_ENV=production
3. bundle exec jekyll build 
4. 'octopress deploy'
5. go to /afs/.cs/group/ai/www/blog/ and pull latest from master

If you don't have octopress, install it with `gem install octopress`. You don't technically need it, you can also just copy the \_site contents after build and commit them to master manually, but octopress deploy is a little shortcut that makes this simpler.

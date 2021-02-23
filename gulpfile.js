/*
 * Here are all of the gulp tasks you can use to help manage your blog
 * Use `npm install` to install all the dependencies located in package.json
 * If you have an issue with sharp, try: `npm rebuild`.
 * Then `gulp default` to minimize css and images.
 */
const gulp = require('gulp');
const concat = require('gulp-concat');
const uglify = require('gulp-uglify');
const pipeline = require('readable-stream').pipeline;
const imagemin = require('gulp-imagemin');
const less = require('gulp-less');
const cleanCSS = require('gulp-clean-css');
const replace = require('gulp-replace');
const webp = require('gulp-webp');
const responsive = require('gulp-responsive');
const fs = require('fs');

//TODO: Use Grunt if Gulp isn't enough

// Create an empty post with today's date
// usage: gulp post -n <title of the post>
gulp.task('post', function (callback) {
  let args = process.argv;
  let title = args[args.length - 1];
  let filename = new Date().toLocaleDateString('en-CA') + '-' + title.replaceAll(' ', '-') + '.md';
  let content = '---\n' +
    'layout: post\n' +
    'title: ' + title + '\n' +
    //'feature-img: "assets/img/"\n' +
    //'thumbnail: "assets/img/thumbnails/"\n' +
    'tags: []\n' +
    '---';
  console.log('[' + new Date().toLocaleTimeString('en-CA', {hour12: false}) + '] File created: _posts/' + filename);
  fs.writeFile(__dirname + '/_posts/' + filename, content, callback);
});

// Minify JS
gulp.task('js', function() {
  return pipeline(
     gulp.src('assets/_js/*.js'),
     concat('main.min.js'),
     uglify({output: {comments: 'some'}}), //will preserve multi-line comments w/ @preserve, @license or @cc_on
     gulp.dest('assets/js')
  );
});

//FIXME: Updated Bootstrap to 4.6
// Isolate Bootstrap
gulp.task('isolate', function() {
  return pipeline(
     gulp.src('assets/_css/bootstrap-iso.less'),
     less({strictMath: 'on'}),
     replace('.bootstrap-iso html', ''),
     replace('.bootstrap-iso body', ''),
     gulp.dest('assets/css/')
  );
});

// Minify Bootstrap CSS
gulp.task('css', function() {
  return pipeline(
     gulp.src('assets/css/bootstrap-iso.css'),
     cleanCSS(),
     concat('bootstrap-iso.min.css'),
     gulp.dest('assets/css/')
  );
});

// Optimize IMGs
gulp.task("img", function() {
  return pipeline(
     gulp.src('assets/_img/**/*.{png,svg,jpg,webp,jpeg,gif}'),
     imagemin(),
     gulp.dest('assets/img/')
  );
});

// Alternative using "sharp" in case "imagemin" does not work.
// Supported formats: heic, heif, jpeg, jpg, png, raw, tiff, webp
gulp.task('sharp_img', function() {
  let settings = {
    quality: 85,
    progressive: true,
    compressionLevel: 6,
  };
  
  return pipeline(
    gulp.src('assets/_img/**/*.{png,jpg,webp,jpeg}'),
    responsive({
      '**/*.*': settings,
      '*.*': settings,
    }),
    gulp.dest('assets/img')
  );
});

// Convert IMGs to WEBP
gulp.task('webp', function() {
  return pipeline(
    gulp.src('assets/_img/**/*.{png,svg,jpg,jpeg,gif}'),
    webp({
      quality: 85,
      preset: 'photo',
      method: 6
    }),
    gulp.dest('assets/img')
  );
});

// Generate thumbnails
gulp.task('thumbnails', function() {
  let settings = {
    width: '50%', //FIXME: Relative size of a non-absolute
    //format: 'jpeg', // convert to jpeg format
  };

  return pipeline(
    gulp.src('assets/_img/feature-img/*'),
    responsive({
      '**/*.*': settings,
      '*.*': settings,
    }),
    gulp.dest('assets/img/thumbnails/feature-img')
  );
});

gulp.task('thumbnails-all', function () {
  let settings = {
      width: '50%', //FIXME: Relative size of a non-absolute
    //format: 'jpeg', // convert to jpeg format
  };

  return pipeline( 
      gulp.src('assets/_img/*.{png,jpg,webp,jpeg}'),
      responsive({'*.*': settings}),
      gulp.dest('assets/img/thumbnails')
  );
});

// Tasks
gulp.task("isolate-bootstrap-css", gulp.series('isolate', 'css'));
gulp.task("default", gulp.series(gulp.parallel('js', 'css', 'img')));


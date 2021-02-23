/*
 * Here are all of the gulp tasks you can use to help manage your blog
 * Use `npm install` to install all the dependencies located in package.json
 * If you have an issue with sharp, try: `npm rebuild`.
 * Then `gulp default` to minimize css and images.
 */
const gulp = require('gulp');
const concat = require('gulp-concat');
const uglify = require('gulp-uglify');
const pump = require('pump');
const imagemin = require('gulp-imagemin');
const less = require('gulp-less');
const cleanCSS = require('gulp-clean-css');
const replace = require('gulp-replace');
const webp = require('gulp-webp');
const responsive = require('gulp-responsive');
const fs = require('fs');
const changed = require('gulp-changed');

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
gulp.task('js', function(cb) {
  pum([
     gulp.src('assets/_js/*.js'),
     concat('main.min.js'),
     uglify({output: {comments: 'some'}}), //will preserve multi-line comments w/ @preserve, @license or @cc_on
     gulp.dest('assets/js')
  ],
  cb
  );
});

//FIXME: Updated Bootstrap to 4.6
// Isolate Bootstrap
gulp.task('bsIsolate', function(cb) {
  pump([
     gulp.src('assets/_css/bootstrap-iso.less'),
     less({strictMath: 'on'}),
     replace('.bootstrap-iso html', ''),
     replace('.bootstrap-iso body', ''),
     gulp.dest('assets/_css/')
  ],
  cb
  );
});

// Minify Bootstrap CSS
gulp.task('bsMinify', function(cb) {
  pump([
     gulp.src('assets/_css/bootstrap-iso.css'),
     cleanCSS(),
     concat('bootstrap-iso.min.css'),
     gulp.dest('assets/css/')
  ],
  cb
  );
});

// Resize IMGs
const paths = {
    featured: {
        src: 'assets/_img/featured/*.{gif,jpg,jpeg,png,svg}',
        dest: 'assets/img/featured/'
    }
}

gulp.task("featured", function(cb) {
  pump([
     gulp.src(paths.featured.src),
     changed(paths.featured.dest),
     responsive({'*': {width: 1920}}),
     gulp.dest(paths.featured.dest)
  ],
  cb
  );
})

// Minify IMGs
gulp.task("imgMinify", function(cb) {
  pump([
     gulp.src('assets/img/featured/*.{gif,jpg,jpeg,png,svg}'),
     imagemin({verbose: true}),
     gulp.dest('assets/img/featured')
  ],
  cb
  );
});

// Convert IMGs to WEBP
gulp.task('webp', function(cb) {
  pump([
    gulp.src('assets/_img/**/*.{png,svg,jpg,jpeg,gif}'),
    webp({quality: 85, preset: 'photo', method: 6}),
    gulp.dest('assets/img')
  ],
  cb
  );
});

// Generate thumbnails
gulp.task('thumbnails', function(cb) {
  let settings = {
    width: '50%', //FIXME: Relative size of a non-absolute
    //format: 'jpeg', // convert to jpeg format
  };

  pump([
    gulp.src('assets/_img/feature-img/*'),
    responsive({
      '**/*.*': settings,
      '*.*': settings,
    }),
    gulp.dest('assets/img/thumbnails/feature-img')
  ],
  cb
  );
});

gulp.task('thumbnails-all', function () {
  let settings = {
      width: '50%', //FIXME: Relative size of a non-absolute
    //format: 'jpeg', // convert to jpeg format
  };

  pump([ 
      gulp.src('assets/_img/*.{png,jpg,webp,jpeg}'),
      responsive({'*.*': settings}),
      gulp.dest('assets/img/thumbnails')
  ],
  cb
  );
});

// Tasks
gulp.task("bootstrap", gulp.series('bsIsolate', 'bsMinify'));
//gulp.task("img", gulp.series('imgResize', 'imgMinify'));
gulp.task("default", gulp.series(gulp.parallel('js', 'bootstrap')));


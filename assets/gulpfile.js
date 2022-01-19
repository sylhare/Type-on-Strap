/*
 * Here are all of the gulp tasks you can use to help manage your blog
 * Use `npm install` to install all the dependencies located in package.json
 * If you have an issue with sharp, try: `npm rebuild`.
 * Then `gulp default` to minimize css and images.
 */
const gulp = require('gulp');
const concat = require('gulp-concat');
const uglify = require('gulp-uglify');
const imagemin = require('gulp-imagemin');
const less = require('gulp-less');
const cleanCSS = require('gulp-clean-css');
const replace = require('gulp-replace');
const webp = require('gulp-webp');
const responsive = require('gulp-responsive');
const fs = require('fs');

// Use it gulp post -n <title of the post>
gulp.task('post', function (callback) {
  let args = process.argv;
  let title = args[args.length - 1];
  let filename = new Date().toLocaleDateString('en-CA') + '-' + title.replace(/ /g, '-') + '.md';
  let content = '---\n' +
    'layout: post\n' +
    'title: ' + title + '\n' +
    //'feature-img: "assets/img/"\n' +
    //'thumbnail: "assets/img/thumbnails/"\n' +
    'tags: []\n' +
    '---';
  console.log('[' + new Date().toLocaleTimeString('en-CA', {hour12: false}) + '] File created: _posts/' + filename);
  fs.writeFile(__dirname + '/../_posts/' + filename, content, callback);
});

gulp.task('js', function minijs() {
  return gulp.src(['js/partials/**.js'])
    .pipe(concat('main.min.js'))
    .pipe(uglify())
    .on('error', (err) => {
      console.log(err.toString())
    })
    .pipe(gulp.dest("js/"))
});

gulp.task("img", function imging() {
  return gulp.src('img/**/*.{png,svg,jpg,webp,jpeg,gif}')
    .pipe(imagemin())
    .on('error', (err) => {
      console.log(err.toString())
    })
    .pipe(gulp.dest('img/'))
});

// Alternative using "sharp" in case "imagemin" does not work, supported formats: heic, heif, jpeg, jpg, png, raw, tiff, webp
gulp.task('sharp_img', function () {
  let settings = {
    quality: 85,
    progressive: true,
    compressionLevel: 6,
  };

  return gulp.src('img/**/*.{png,jpg,webp,jpeg}')
    .pipe(responsive({
      '**/*.*': settings,
      '*.*': settings,
    }))
    .pipe(gulp.dest('img'))
});

gulp.task('thumbnails', function () {
  let settings = {
    width: '50%',
    //format: 'jpeg', // convert to jpeg format
  };

  return gulp.src('img/feature-img/*')
    .pipe(responsive({
      '**/*.*': settings,
      '*.*': settings,
    }))
    .pipe(gulp.dest('img/thumbnails/feature-img'))
});


gulp.task('thumbnails-all', function () {
  let settings = {
    width: '50%',
    //format: 'jpeg', // convert to jpeg format
  };

  return gulp.src('img/*.{png,jpg,webp,jpeg}')
      .pipe(responsive({'*.*': settings}))
      .pipe(gulp.dest('img/thumbnails')) &&
    gulp.src('img/!(thumbnails)/*.{png,jpg,webp,jpeg}')
      .pipe(responsive({'**/*.*': settings}))
      .pipe(gulp.dest('img/thumbnails'))
});

gulp.task('webp', () =>
  gulp.src('img/**/*.{png,svg,jpg,jpeg,gif}')
    .pipe(webp({
      quality: 85,
      preset: 'photo',
      method: 6
    }))
    .pipe(gulp.dest('img'))
);

gulp.task('css', function minicss() {
  return gulp.src('css/vendor/bootstrap-iso.css')
    .pipe(cleanCSS())
    .on('error', (err) => {
      console.log(err.toString())
    })
    .pipe(concat('bootstrap-iso.min.css'))
    .pipe(gulp.dest('css/vendor/'));
});

gulp.task('isolate', function isolateBootstrap() {
  return gulp.src('css/bootstrap-iso.less')
    .pipe(less({strictMath: 'on'}))
    .pipe(replace('.bootstrap-iso html', ''))
    .pipe(replace('.bootstrap-iso body', ''))
    .pipe(gulp.dest('css/vendor/'));
});

gulp.task("isolate-bootstrap-css", gulp.series('isolate', 'css'));
gulp.task("default", gulp.series(gulp.parallel('js', 'css', 'img')));

/*
 * Here are all of the gulp tasks you can use to help manage your blog
 * Use `npm install` to install all the dependencies located in package.json 
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
const responsive = require('gulp-responsive'); // Supported formats: heic, heif, jpeg, jpg, png, raw, tiff, webp

gulp.task('js', function minijs() {
  return gulp.src(['js/partials/**.js'])
    .pipe(concat('main.min.js'))
    .pipe(uglify())
    .on('error', (err) => { console.log(err.toString()) })
    .pipe(gulp.dest("js/"))
});

gulp.task("img", function imging() {
  return gulp.src('img/**/*.{png,svg,jpg,webp,jpeg,gif}')
    .pipe(imagemin())
    .on('error', (err) => { console.log(err.toString()) })
    .pipe(gulp.dest('img/'))
});

// Alternative using "sharp" in case "imagemin" does not work
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

  return gulp.src('img/**/*.{png,jpg,webp,jpeg}')
    .pipe(responsive({
      '**/*.*': settings,
      '*.*': settings,
    }))
    .pipe(gulp.dest('thumbnails'))
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
    .on('error', (err) => { console.log(err.toString()) })
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

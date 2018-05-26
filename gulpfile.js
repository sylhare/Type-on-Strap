const gulp = require('gulp');
const concat = require('gulp-concat');
const uglify = require('gulp-uglify');
const imagemin = require('gulp-imagemin');
const gutil = require('gulp-util');
const shell = require("gulp-shell")

// Concatenates and uglifies global JS files
gulp.task('js', function () {
    gulp.src(['assets/js/partials/**.js'])
        .pipe(concat('main.min.js'))
        .pipe(uglify())
        .on('error', (err) => {
            gutil.log(gutil.colors.red('[Error]'), err.toString());
        })
        .pipe(gulp.dest("assets/js/"))
});

gulp.task("img", function () {
    gulp.src('assets/img/**/*.{png,svg,jpg,gif}')
        .pipe(imagemin())
        .on('error', (err) => {
            gutil.log(gutil.colors.red('[Error]'), err.toString());
        })
        .pipe(gulp.dest('assets/img/'))
});

gulp.task("serve", shell.task([
  "python -m webbrowser 'http://localhost:4000/Type-on-Strap/' && bundle exec jekyll serve --watch"
]));

gulp.task('default', ['js', 'img']);

/* 
 * Use `npm install` to install all the dependencies located in package.json 
 * Then `gulp default` to minimize css and images.
 */
const gulp = require('gulp');
const concat = require('gulp-concat');
const uglify = require('gulp-uglify');
const imagemin = require('gulp-imagemin');
const shell = require('gulp-shell');
const less = require('gulp-less');
const cssmin = require('gulp-cssmin')
const replace = require('gulp-replace');

gulp.task('js', function minijs() {
    return gulp.src(['js/partials/**.js'])
        .pipe(concat('main.min.js'))
        .pipe(uglify())
        .on('error', (err) => {
            console.log(err.toString());
        })
        .pipe(gulp.dest("js/"))
});

gulp.task("img", function imging() {
    return gulp.src('img/**/*.{png,svg,jpg,gif}')
        .pipe(imagemin())
        .on('error', (err) => {
            console.log(err.toString());
        })
        .pipe(gulp.dest('img/'))
});

gulp.task('css', function minicss() {
    return gulp.src('css/vendor/bootstrap-iso.css')
        .pipe(cssmin())
        .on('error', (err) => {
            console.log(err.toString());
        })
        .pipe(concat('bootstrap-iso.min.css'))
        .pipe(gulp.dest('css/vendor/'));
})

gulp.task("isolate-bootstrap-css", gulp.series('css', function isolating() {
    return gulp.src('css/bootstrap-iso.less')
        .pipe(less())
        .pipe(replace('.bootstrap-iso html', ''))
        .pipe(replace('.bootstrap-iso body', ''))
        .pipe(gulp.dest('css/vendor/'));
}));

gulp.task("serve", function serving(done) {
    console.log('... not working at the moment try \ncd .. && bundle exec jekyll serve --watch\n',
                'then go to \nhttp://localhost:4000/Type-on-Strap/');
    shell.task([
        "python -m webbrowser 'http://localhost:4000/Type-on-Strap/'; cd .. && bundle exec jekyll serve --watch"
    ]);
    done();
});

gulp.task("default", gulp.series(gulp.parallel('js', 'css', 'img')));

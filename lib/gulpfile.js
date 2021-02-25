/*
 * Here are all of the gulp tasks you can use to help manage your blog
 * Use `npm install` to install all the dependencies located in package.json
 * If you have an issue with sharp, try: `npm rebuild`.
 * Then `gulp default` to minimize css and images.
 */
const { src, dest, series, parallel } = require('gulp');
const pump = require('pump');
const concat = require('gulp-concat');
const uglify = require('gulp-uglify');
const imagemin = require('gulp-imagemin');
const less = require('gulp-less');
const cleanCSS = require('gulp-clean-css');
const replace = require('gulp-replace');
const webp = require('gulp-webp');
const responsive = require('gulp-responsive');
const fs = require('fs');
const changed = require('gulp-changed');

// Paths
const postsPath =  '../_posts/'
const assetsPath = '../assets/'
const paths = {
    js: {
        src: assetsPath+'_js/*.js',
        dest: assetsPath+'js/'
    },
    css: {
        src: assetsPath+'_css/',
        dest: assetsPath+'css/'
    },
    avatar: {
        src: assetsPath+'_img/avatar.{gif,jpg,jpeg,png,svg}',
        dest: assetsPath+'img/'
    },
    featured: {
        src: assetsPath+'_img/featured/*.{gif,jpg,jpeg,png,svg}',
        dest: assetsPath+'img/featured/'
    },
    thumbnails: {
        src: assetsPath+'img/featured/*.{gif,jpg,jpeg,png,svg,webp}',
        dest: assetsPath+'img/thumbnails/'
    },
    posts: {
        src: assetsPath+'_img/posts/*.{gif,jpg,jpeg,png,svg}',
        dest: assetsPath+'img/posts/'
    },
    portfolio: {
        src: assetsPath+'_img/portfolio/*.{gif,jpg,jpeg,png,svg}',
        dest: assetsPath+'img/portfolio/'
    },
     webp: {
        src: assetsPath+'img/featured/*.{jpg,jpeg,png}',
        dest: assetsPath+'img/featured/webp/'
    },
}

// Create an empty post with today's date
// usage: gulp post -n <title of the post>
const post = function(callback) {
  let args = process.argv;
  let title = args[args.length - 1];
  let filename = new Date().toLocaleDateString('en-CA') + '-' + title.replaceAll(' ', '-') + '.md';
  let content = '---\n' +
    'layout: post\n' +
    'title: ' + title + '\n' +
    '#author:\n' +
    'color: rgb(80,140,22)\n' +
    '#feature-img: assets/img/template.jpg\n' +
    '#thumbnail: assets/img/thumbnails/template.jpg\n' +
    'tags: []\n' +
    '#permalink:\n' +    
    'hide_title: false\n' +
    'hidden: false\n' +
    'bootstrap: false\n' + 
    '#excerpt_separator: <!--more-->\n' +
    '---'+
    '\nWrite here...';
  console.log('[' + new Date().toLocaleTimeString('en-CA', {hour12: false}) + '] File created: _posts/' + filename);
  fs.writeFile(postsPath + filename, content, callback);
}

// Minify JS
const js = function(cb) {
  pump([
     src(paths.js.src),
     concat('main.min.js'),
     uglify({output: {comments: 'some'}}), //will preserve multi-line comments w/ @preserve, @license or @cc_on
     dest(paths.js.dest)
  ],
  cb()
  );
}

// TODO: Updated Bootstrap to 4.6
// Isolate Bootstrap
const bsIsolate = function(cb) {
  pump([
     src(paths.css.src+'bootstrap-iso.less'),
     less({strictMath: 'on'}),
     replace('.bootstrap-iso html', ''),
     replace('.bootstrap-iso body', ''),
     dest(paths.css.src)
  ],
  cb()
  );
}

// Minify Bootstrap CSS
const bsMinify = function(cb) {
  pump([
     src(paths.css.src+'bootstrap-iso.css'),
     cleanCSS(),
     concat('bootstrap-iso.min.css'),
     dest(paths.css.dest)
  ],
  cb()
  );
}

// Resize and Minify IMG
const imgAvatar = function(cb) {
  pump([
     src(paths.avatar.src),
     changed(paths.avatar.dest),
     responsive({'*': {width: 512}}),
     imagemin({verbose: true}),
     dest(paths.avatar.dest)
  ],
  cb()
  );
}

const imgFeatured = function(cb) {
  pump([
     src(paths.featured.src),
     changed(paths.featured.dest),
     responsive({'*': {width: 1920}}),
     imagemin({verbose: true}),
     dest(paths.featured.dest)
  ],
  cb()
  );
}

const imgThumbnails = function(cb) {
  pump([
     src(paths.thumbnails.src),
     changed(paths.thumbnails.dest),
     responsive({'*': {width: 1024}}),
     dest(paths.thumbnails.dest)
  ],
  cb()
  );
}

const imgPosts = function(cb) {
  pump([
     src(paths.posts.src),
     changed(paths.posts.dest),
     //responsive({'*': {width: 900}}),
     imagemin({verbose: true}),
     dest(paths.posts.dest)
  ],
  cb()
  );
}

const imgPortfolio = function(cb) {
  pump([
     src(paths.portfolio.src),
     changed(paths.portfolio.dest),
     responsive({'*': {width: 900}}),
     imagemin({verbose: true}),
     dest(paths.portfolio.dest)
  ],
  cb()
  );
}

const imgWebp = function(cb) {
  pump([
     src(paths.webp.src),
     changed(paths.webp.dest),
     webp({quality: 85, preset: 'photo', method: 6}),
     dest(paths.webp.dest)
  ],
  cb()
  );
}

// Tasks
exports.post = post;
exports.js = js;
exports.bootstrap = series(bsIsolate, bsMinify);

exports.avatar = imgAvatar;
exports.featured = imgFeatured;
exports.thumbs = imgThumbnails;
exports.posts = imgPosts;
exports.portfolio = imgPortfolio;
exports.webp = imgWebp;
exports.img = parallel(imgAvatar, imgPortfolio, series(imgFeatured, imgThumbnails));

// FIXME: Gulp will generate thumbnails only on 2nd run
exports.default = parallel(js, series(bsIsolate, bsMinify), imgAvatar, imgPortfolio, imgPosts, series(imgFeatured, imgThumbnails));


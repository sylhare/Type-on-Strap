'use strict';
const fs = require('fs');
const path = require('path');

function formatDate(date) {
  return date.toLocaleDateString('en-CA');
}

function formatTime(date) {
  return date.toLocaleTimeString('en-CA', { hour12: false });
}

function createFilename(title, date) {
  return formatDate(date) + '-' + title.replace(/ /g, '-') + '.md';
}

function createContent(title) {
  return '---\n' +
    'layout: post\n' +
    'title: ' + title + '\n' +
    //'feature-img: "assets/img/"\n' +
    //'thumbnail: "assets/img/thumbnails/"\n' +
    'tags: []\n' +
    '---';
}

function createPost(title, postsDir, date = new Date()) {
  const filename = createFilename(title, date);
  const filepath = path.join(postsDir, filename);
  const content = createContent(title);
  fs.writeFileSync(filepath, content);
  console.log('[' + formatTime(date) + '] File created: _posts/' + filename);
  return { filename, filepath, content };
}

module.exports = { formatDate, formatTime, createFilename, createContent, createPost };

if (require.main === module) {
  const title = process.argv.slice(2).join(' ');
  if (!title) {
    console.error('Usage: npm run post "Post Title"');
    process.exit(1);
  }
  const postsDir = path.join(process.cwd(), '_posts');
  createPost(title, postsDir);
}

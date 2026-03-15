'use strict';
const fs = require('fs');
const path = require('path');
const { globSync } = require('glob');
const esbuild = require('esbuild');
const less = require('less');
const CleanCSS = require('clean-css');

function getJsPartials(dir) {
  return globSync(path.join(dir, '*.js')).sort();
}

function concatFiles(filePaths) {
  return filePaths.map(f => fs.readFileSync(f, 'utf8')).join('\n');
}

async function buildJs(partialFiles, outFile) {
  const code = concatFiles(partialFiles);
  const result = await esbuild.transform(code, { minify: true });
  fs.writeFileSync(outFile, result.code);
}

async function compileLess(inputPath, outputPath) {
  const input = fs.readFileSync(inputPath, 'utf8');
  const result = await less.render(input, { filename: inputPath, strictMath: true });
  let css = result.css.replace(/\.bootstrap-iso (?:html|body)/g, '');
  fs.writeFileSync(outputPath, css);
  return css;
}

function minifyCSS(css, outputPath) {
  const result = new CleanCSS().minify(css);
  fs.writeFileSync(outputPath, result.styles);
}

module.exports = { getJsPartials, concatFiles, buildJs, compileLess, minifyCSS };

if (require.main === module) {
  const cwd = process.cwd();
  const arg = process.argv[2];

  async function runJs() {
    const partialsDir = path.join(cwd, 'assets/js/partials');
    const partials = getJsPartials(partialsDir);
    const commentsFile = path.join(cwd, 'assets/js/comments-lazy-load.js');
    await Promise.all([
      buildJs(partials, path.join(cwd, 'assets/js/main.min.js'))
        .then(() => console.log('Built assets/js/main.min.js')),
      buildJs([commentsFile], path.join(cwd, 'assets/js/comments-lazy-load.min.js'))
        .then(() => console.log('Built assets/js/comments-lazy-load.min.js')),
    ]);
  }

  async function runCss() {
    const lessIn = path.join(cwd, 'assets/css/bootstrap-iso.less');
    const cssOut = path.join(cwd, 'assets/css/vendor/bootstrap-iso.css');
    const css = await compileLess(lessIn, cssOut);
    console.log('Compiled assets/css/vendor/bootstrap-iso.css');

    minifyCSS(css, path.join(cwd, 'assets/css/vendor/bootstrap-iso.min.css'));
    console.log('Built assets/css/vendor/bootstrap-iso.min.css');
  }

  (async () => {
    switch (arg) {
      case 'js':  await runJs();  break;
      case 'css': await runCss(); break;
      default:
        await runJs();
        await runCss();
    }
  })().catch(err => {
    console.error(err);
    process.exit(1);
  });
}

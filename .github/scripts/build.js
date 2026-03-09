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
  let css = result.css;
  css = css.replace(/\.bootstrap-iso html/g, '');
  css = css.replace(/\.bootstrap-iso body/g, '');
  fs.writeFileSync(outputPath, css);
}

function minifyCSS(inputPath, outputPath) {
  const input = fs.readFileSync(inputPath, 'utf8');
  const result = new CleanCSS().minify(input);
  fs.writeFileSync(outputPath, result.styles);
}

module.exports = { getJsPartials, concatFiles, buildJs, compileLess, minifyCSS };

if (require.main === module) {
  const cwd = process.cwd();
  const arg = process.argv[2];

  async function runJs() {
    const partialsDir = path.join(cwd, 'assets/js/partials');
    const partials = getJsPartials(partialsDir);
    await buildJs(partials, path.join(cwd, 'assets/js/main.min.js'));
    console.log('Built assets/js/main.min.js');

    const commentsFile = path.join(cwd, 'assets/js/comments-lazy-load.js');
    await buildJs([commentsFile], path.join(cwd, 'assets/js/comments-lazy-load.min.js'));
    console.log('Built assets/js/comments-lazy-load.min.js');
  }

  async function runCss() {
    const lessIn = path.join(cwd, 'assets/css/bootstrap-iso.less');
    const cssOut = path.join(cwd, 'assets/css/vendor/bootstrap-iso.css');
    await compileLess(lessIn, cssOut);
    console.log('Compiled assets/css/vendor/bootstrap-iso.css');

    minifyCSS(cssOut, path.join(cwd, 'assets/css/vendor/bootstrap-iso.min.css'));
    console.log('Built assets/css/vendor/bootstrap-iso.min.css');
  }

  (async () => {
    if (arg === 'js') {
      await runJs();
    } else if (arg === 'css') {
      await runCss();
    } else {
      await runJs();
      await runCss();
    }
  })().catch(err => {
    console.error(err);
    process.exit(1);
  });
}

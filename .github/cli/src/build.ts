import * as fs from 'node:fs';
import * as path from 'node:path';
import { globSync } from 'glob';
import * as esbuild from 'esbuild';
import * as less from 'less';
import CleanCSS from 'clean-css';
import { logger } from './utils/logger';

export function getJsPartials(dir: string): string[] {
  return globSync(path.join(dir, '*.js')).sort();
}

export function concatFiles(filePaths: string[]): string {
  return filePaths.map(f => fs.readFileSync(f, 'utf8')).join('\n');
}

export async function buildJs(partialFiles: string[], outFile: string): Promise<void> {
  const code = concatFiles(partialFiles);
  const result = await esbuild.transform(code, { minify: true });
  fs.writeFileSync(outFile, result.code);
}

export async function compileLess(inputPath: string, outputPath: string): Promise<string> {
  const input = fs.readFileSync(inputPath, 'utf8');
  const result = await less.render(input, { filename: inputPath, strictMath: true } as Less.Options);
  let css = result.css.replace(/\.bootstrap-iso (?:html|body)/g, '');
  fs.writeFileSync(outputPath, css);
  return css;
}

export function minifyCSS(css: string, outputPath: string): void {
  const result = new CleanCSS().minify(css);
  fs.writeFileSync(outputPath, result.styles);
}

if (require.main === module) {
  const cwd = process.cwd();
  const arg = process.argv[2];

  async function runJs(): Promise<void> {
    const partialsDir = path.join(cwd, 'assets/js/partials');
    const partials = getJsPartials(partialsDir);
    const commentsFile = path.join(cwd, 'assets/js/comments-lazy-load.js');
    await Promise.all([
      buildJs(partials, path.join(cwd, 'assets/js/main.min.js'))
        .then(() => logger.info('Built assets/js/main.min.js')),
      buildJs([commentsFile], path.join(cwd, 'assets/js/comments-lazy-load.min.js'))
        .then(() => logger.info('Built assets/js/comments-lazy-load.min.js')),
    ]);
  }

  async function runCss(): Promise<void> {
    const lessIn = path.join(cwd, 'assets/css/bootstrap-iso.less');
    const cssOut = path.join(cwd, 'assets/css/vendor/bootstrap-iso.css');
    const css = await compileLess(lessIn, cssOut);
    logger.info('Compiled assets/css/vendor/bootstrap-iso.css');

    minifyCSS(css, path.join(cwd, 'assets/css/vendor/bootstrap-iso.min.css'));
    logger.info('Built assets/css/vendor/bootstrap-iso.min.css');
  }

  (async () => {
    switch (arg) {
      case 'js':
        await runJs();
        break;
      case 'css':
        await runCss();
        break;
      default:
        await Promise.all([runJs(), runCss()]);
    }
  })().catch(err => {
    logger.error(String(err));
    process.exit(1);
  });
}

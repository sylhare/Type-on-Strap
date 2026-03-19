import fs from 'node:fs';
import path from 'node:path';
import os from 'node:os';
import { buildJs, compileLess, minifyCSS, getJsPartials } from '../../src/build';

const ROOT = path.resolve(__dirname, '../../../../');
const TMP = fs.mkdtempSync(path.join(os.tmpdir(), 'tos-build-'));

afterAll(() => {
  fs.rmSync(TMP, { recursive: true, force: true });
});

describe('build:js integration', () => {
  const outFile = path.join(TMP, 'main.min.js');

  beforeAll(async () => {
    const partials = getJsPartials(path.join(ROOT, 'assets/js/partials'));
    await buildJs(partials, outFile);
  });

  test('creates output file', () => {
    expect(fs.existsSync(outFile)).toEqual(true);
  });

  test('output is non-empty', () => {
    expect(fs.statSync(outFile).size).toBeGreaterThan(0);
  });

  test('output is minified (single line)', () => {
    const content = fs.readFileSync(outFile, 'utf8').trim();
    expect(content.split('\n').length).toEqual(1);
  });

  test('output is smaller than concatenated source', () => {
    const partials = getJsPartials(path.join(ROOT, 'assets/js/partials'));
    const srcSize = partials.reduce((sum, f) => sum + fs.statSync(f).size, 0);
    expect(fs.statSync(outFile).size).toBeLessThan(srcSize);
  });
});

describe('build:css integration', () => {
  const lessIn = path.join(ROOT, 'assets/css/bootstrap-iso.less');
  const cssOut = path.join(TMP, 'bootstrap-iso.css');
  const minOut = path.join(TMP, 'bootstrap-iso.min.css');

  beforeAll(async () => {
    const css = await compileLess(lessIn, cssOut);
    minifyCSS(css, minOut);
  });

  test('compileLess creates the CSS output file', () => {
    expect(fs.existsSync(cssOut)).toEqual(true);
  });

  test('compiled CSS does not contain .bootstrap-iso html', () => {
    const css = fs.readFileSync(cssOut, 'utf8');
    expect(css).not.toContain('.bootstrap-iso html');
  });

  test('compiled CSS does not contain .bootstrap-iso body', () => {
    const css = fs.readFileSync(cssOut, 'utf8');
    expect(css).not.toContain('.bootstrap-iso body');
  });

  test('compiled CSS still scopes rules under .bootstrap-iso', () => {
    const css = fs.readFileSync(cssOut, 'utf8');
    expect(css).toContain('.bootstrap-iso');
  });

  test('minified CSS file is created', () => {
    expect(fs.existsSync(minOut)).toEqual(true);
  });

  test('minified CSS is smaller than compiled CSS', () => {
    expect(fs.statSync(minOut).size).toBeLessThan(fs.statSync(cssOut).size);
  });
});

import fs from 'node:fs';
import * as path from 'node:path';
import * as os from 'node:os';
import { buildJs, getJsPartials } from '../../src/build';

const ROOT = path.resolve(__dirname, '../../../../');
const TMP = fs.mkdtempSync(path.join(os.tmpdir(), 'tos-build-'));

afterAll(() => {
  fs.rmSync(TMP, { recursive: true, force: true });
});

describe('build:js integration', () => {
  const outFile = path.join(TMP, 'main.min.js');
  let partials: string[];

  beforeAll(async () => {
    partials = getJsPartials(path.join(ROOT, 'assets/js/partials'));
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
    const srcSize = partials.reduce((sum, f) => sum + fs.statSync(f).size, 0);
    expect(fs.statSync(outFile).size).toBeLessThan(srcSize);
  });
});

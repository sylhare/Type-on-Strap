import fs from 'node:fs';
import path from 'node:path';

jest.mock('esbuild', () => ({ transform: jest.fn() }));
jest.mock('less', () => ({ render: jest.fn() }));
jest.mock('clean-css', () =>
  jest.fn().mockImplementation(() => ({ minify: jest.fn().mockReturnValue({ styles: '.a{}' }) }))
);
jest.mock('glob', () => ({ globSync: jest.fn() }));

import * as esbuild from 'esbuild';
import less from 'less';
import CleanCSS from 'clean-css';
import { globSync } from 'glob';
import { getJsPartials, concatFiles, buildJs, compileLess, minifyCSS } from '../src/build';

const mockEsbuild = esbuild as jest.Mocked<typeof esbuild>;
const mockLess = less as jest.Mocked<typeof less>;
const mockGlobSync = globSync as jest.Mock;
const MockCleanCSS = CleanCSS as jest.Mock;

describe('build.ts', () => {
  describe('getJsPartials()', () => {
    test('returns a sorted file list from the partials directory', () => {
      mockGlobSync.mockReturnValue(['/dir/c.js', '/dir/a.js', '/dir/b.js']);
      const result = getJsPartials('/dir');
      expect(result).toEqual(['/dir/a.js', '/dir/b.js', '/dir/c.js']);
    });

    test('passes the correct glob pattern', () => {
      mockGlobSync.mockReturnValue([]);
      getJsPartials('/assets/js/partials');
      expect(mockGlobSync).toHaveBeenCalledWith(
        path.join('/assets/js/partials', '*.js')
      );
    });

    test('returns empty array when no files found', () => {
      mockGlobSync.mockReturnValue([]);
      expect(getJsPartials('/empty')).toEqual([]);
    });
  });

  describe('concatFiles()', () => {
    test('joins file contents with newline', () => {
      jest.spyOn(fs, 'readFileSync')
        .mockReturnValueOnce('content A')
        .mockReturnValueOnce('content B');
      const result = concatFiles(['/a.js', '/b.js']);
      expect(result).toEqual('content A\ncontent B');
    });

    test('reads each file as utf8', () => {
      const spy = jest.spyOn(fs, 'readFileSync').mockReturnValue('');
      concatFiles(['/a.js']);
      expect(spy).toHaveBeenCalledWith('/a.js', 'utf8');
    });

    test('returns empty string for empty file list', () => {
      expect(concatFiles([])).toEqual('');
    });
  });

  describe('buildJs()', () => {
    test('transforms concatenated source with minify enabled', async () => {
      jest.spyOn(fs, 'readFileSync').mockReturnValue('function foo() {}');
      jest.spyOn(fs, 'writeFileSync').mockImplementation(() => {});
      mockEsbuild.transform.mockResolvedValue({ code: 'function foo(){}' } as any);

      await buildJs(['/partials/a.js'], '/out/main.min.js');

      expect(mockEsbuild.transform).toHaveBeenCalledWith(
        expect.any(String),
        { minify: true }
      );
    });
  });

  describe('compileLess()', () => {
    test('passes filename to less.render for import resolution', async () => {
      jest.spyOn(fs, 'readFileSync').mockReturnValue('');
      jest.spyOn(fs, 'writeFileSync').mockImplementation(() => {});
      mockLess.render.mockResolvedValue({ css: '' } as any);

      await compileLess('/path/to/input.less', '/output.css');

      expect(mockLess.render).toHaveBeenCalledWith(
        expect.any(String),
        expect.objectContaining({ filename: '/path/to/input.less' })
      );
    });
  });

  describe('minifyCSS()', () => {
    test('writes minified output to the given path', () => {
      jest.spyOn(fs, 'writeFileSync').mockImplementation(() => {});
      MockCleanCSS.mockImplementation(() => ({ minify: jest.fn().mockReturnValue({ styles: '.a{}' }) }));

      minifyCSS('.a { color: red; }', '/out/file.min.css');

      expect(fs.writeFileSync).toHaveBeenCalledWith('/out/file.min.css', '.a{}');
    });
  });
});

/**
 * @fileoverview Unit tests for build.js
 * @jest-environment node
 */

'use strict';
const fs = require('fs');
const path = require('path');

jest.mock('esbuild', () => ({ transform: jest.fn() }), { virtual: true });
jest.mock('less', () => ({ render: jest.fn() }), { virtual: true });
jest.mock('clean-css', () => jest.fn(), { virtual: true });
jest.mock('glob', () => ({ globSync: jest.fn() }), { virtual: true });

const esbuild = require('esbuild');
const less = require('less');
const CleanCSS = require('clean-css');
const { globSync } = require('glob');

const { getJsPartials, concatFiles, buildJs, compileLess, minifyCSS } =
  require('../../scripts/build');

describe('build.js', () => {
  describe('getJsPartials()', () => {
    test('returns a sorted file list from the partials directory', () => {
      globSync.mockReturnValue(['/dir/c.js', '/dir/a.js', '/dir/b.js']);
      const result = getJsPartials('/dir');
      expect(result).toEqual(['/dir/a.js', '/dir/b.js', '/dir/c.js']);
    });

    test('passes the correct glob pattern', () => {
      globSync.mockReturnValue([]);
      getJsPartials('/assets/js/partials');
      expect(globSync).toHaveBeenCalledWith(
        path.join('/assets/js/partials', '*.js')
      );
    });

    test('returns empty array when no files found', () => {
      globSync.mockReturnValue([]);
      expect(getJsPartials('/empty')).toEqual([]);
    });
  });

  describe('concatFiles()', () => {
    test('joins file contents with newline', () => {
      jest.spyOn(fs, 'readFileSync')
        .mockReturnValueOnce('content A')
        .mockReturnValueOnce('content B');
      const result = concatFiles(['/a.js', '/b.js']);
      expect(result).toBe('content A\ncontent B');
    });

    test('reads each file as utf8', () => {
      const spy = jest.spyOn(fs, 'readFileSync').mockReturnValue('');
      concatFiles(['/a.js']);
      expect(spy).toHaveBeenCalledWith('/a.js', 'utf8');
    });

    test('returns empty string for empty file list', () => {
      expect(concatFiles([])).toBe('');
    });
  });

  describe('buildJs()', () => {
    test('transforms concatenated source with minify enabled', async () => {
      jest.spyOn(fs, 'readFileSync').mockReturnValue('function foo() {}');
      jest.spyOn(fs, 'writeFileSync').mockImplementation(() => {});
      esbuild.transform.mockResolvedValue({ code: 'function foo(){}' });

      await buildJs(['/partials/a.js'], '/out/main.min.js');

      expect(esbuild.transform).toHaveBeenCalledWith(
        expect.any(String),
        { minify: true }
      );
    });

  });

  describe('compileLess()', () => {
    test('passes filename to less.render for import resolution', async () => {
      jest.spyOn(fs, 'readFileSync').mockReturnValue('');
      jest.spyOn(fs, 'writeFileSync').mockImplementation(() => {});
      less.render.mockResolvedValue({ css: '' });

      await compileLess('/path/to/input.less', '/output.css');

      expect(less.render).toHaveBeenCalledWith(
        expect.any(String),
        expect.objectContaining({ filename: '/path/to/input.less' })
      );
    });
  });

});


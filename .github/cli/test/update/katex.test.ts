jest.mock('../../src/utils/fs');

import { resolveVersion, generateScss } from '../../src/update/katex';
import { readVendorVersion } from '../../src/utils/fs';

const mockReadVendorVersion = readVendorVersion as jest.MockedFunction<typeof readVendorVersion>;

describe('update/katex', () => {
  describe('resolveVersion()', () => {
    test('returns the provided argument', () => {
      expect(resolveVersion('0.16.0')).toEqual('0.16.0');
    });

    test('reads version from vendor config when no argument provided', () => {
      mockReadVendorVersion.mockReturnValue('0.16.38');
      expect(resolveVersion()).toEqual('0.16.38');
    });
  });

  describe('generateScss()', () => {
    test('adds the katex font path header', () => {
      const result = generateScss('');
      expect(result).toContain('$katex-font-path: "../../assets/fonts/katex"');
    });

    test('removes stylelint-disable comments', () => {
      const result = generateScss('/* stylelint-disable some-rule */\n.foo {}');
      expect(result).not.toContain('stylelint-disable some-rule');
    });

    test('removes TTF format entries', () => {
      const css = 'src: url("fonts/KaTeX_Main-Regular.woff2") format("woff2"), url(KaTeX_Main-Regular.ttf) format(\'truetype\')';
      const result = generateScss(css);
      expect(result).not.toContain('.ttf');
    });

    test('replaces unquoted font URLs with SCSS variable', () => {
      const result = generateScss('src: url(fonts/KaTeX_Main-Regular.woff2)');
      expect(result).toContain('url("#{$katex-font-path}/KaTeX_Main-Regular.woff2")');
    });

    test('replaces quoted font URLs with SCSS variable', () => {
      const result = generateScss('src: url("fonts/KaTeX_Main-Regular.woff2")');
      expect(result).toContain('url("#{$katex-font-path}/KaTeX_Main-Regular.woff2")');
    });
  });
});

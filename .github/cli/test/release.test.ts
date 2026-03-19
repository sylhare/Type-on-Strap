import fs from 'node:fs';
import {
  parseVersion,
  bumpVersion,
  updateGemspec,
  updateDefaultHtml,
  updatePackageJson,
  updatePackageLockJson,
  updateGemBuildWorkflow,
  release,
} from '../src/release';

describe('release.ts', () => {
  describe('parseVersion()', () => {
    test('parses a valid semver string', () => {
      expect(parseVersion('2.5.0')).toEqual([2, 5, 0]);
    });

    test('parses major.minor.patch correctly', () => {
      expect(parseVersion('1.23.456')).toEqual([1, 23, 456]);
    });

    test('throws on invalid version format', () => {
      expect(() => parseVersion('2.5')).toThrow('Invalid version format: 2.5');
    });
  });

  describe('bumpVersion()', () => {
    test('bumps patch by default', () => {
      expect(bumpVersion('2.5.0', 'patch')).toEqual('2.5.1');
    });

    test('bumps minor and resets patch', () => {
      expect(bumpVersion('2.5.3', 'minor')).toEqual('2.6.0');
    });

    test('bumps major and resets minor and patch', () => {
      expect(bumpVersion('2.5.3', 'major')).toEqual('3.0.0');
    });

    test('accepts an explicit version string', () => {
      expect(bumpVersion('2.5.0', '3.1.4')).toEqual('3.1.4');
    });

    test('defaults to patch for unknown bump type', () => {
      expect(bumpVersion('2.5.0', 'unknown')).toEqual('2.5.1');
    });
  });

  describe('updateGemspec()', () => {
    test('replaces version in spec.version line', () => {
      jest.spyOn(fs, 'readFileSync').mockReturnValue('spec.version = "2.5.0"\n' as any);
      const write = jest.spyOn(fs, 'writeFileSync').mockImplementation(() => {});

      updateGemspec('/fake/type-on-strap.gemspec', '2.6.0');

      expect(write).toHaveBeenCalledWith('/fake/type-on-strap.gemspec', 'spec.version = "2.6.0"\n');
    });
  });

  describe('updateDefaultHtml()', () => {
    test('replaces version in theme comment', () => {
      const content = '    Type on Strap jekyll theme v2.5.0\n';
      jest.spyOn(fs, 'readFileSync').mockReturnValue(content as any);
      const write = jest.spyOn(fs, 'writeFileSync').mockImplementation(() => {});

      updateDefaultHtml('/fake/_layouts/default.html', '2.6.0');

      expect(write).toHaveBeenCalledWith(
        '/fake/_layouts/default.html',
        '    Type on Strap jekyll theme v2.6.0\n'
      );
    });
  });

  describe('updatePackageJson()', () => {
    test('updates version field in package.json', () => {
      const pkg = JSON.stringify({ name: 'test', version: '2.5.0' }, null, 2);
      jest.spyOn(fs, 'readFileSync').mockReturnValue(pkg as any);
      jest.spyOn(fs, 'writeFileSync').mockImplementation(() => {});

      updatePackageJson('/fake/package.json', '2.6.0');

      expect(fs.writeFileSync).toHaveBeenCalledWith(
        '/fake/package.json',
        expect.stringContaining('"version": "2.6.0"')
      );
    });
  });

  describe('updatePackageLockJson()', () => {
    test('updates root version and packages[""] version', () => {
      const lock = JSON.stringify({ version: '2.5.0', packages: { '': { version: '2.5.0' } } }, null, 2);
      jest.spyOn(fs, 'readFileSync').mockReturnValue(lock as any);
      jest.spyOn(fs, 'writeFileSync').mockImplementation(() => {});

      updatePackageLockJson('/fake/package-lock.json', '2.6.0');

      const writtenContent = (fs.writeFileSync as jest.Mock).mock.calls
        .find((c: unknown[]) => c[0] === '/fake/package-lock.json')?.[1] as string;
      const written = JSON.parse(writtenContent);
      expect(written.version).toEqual('2.6.0');
      expect(written.packages[''].version).toEqual('2.6.0');
    });
  });

  describe('updateGemBuildWorkflow()', () => {
    test('replaces version in gem install command', () => {
      const content = 'gem install type-on-strap --version "2.4.11" --source "https://rubygems.pkg.github.com/sylhare"\n';
      jest.spyOn(fs, 'readFileSync').mockReturnValue(content as any);
      const write = jest.spyOn(fs, 'writeFileSync').mockImplementation(() => {});

      updateGemBuildWorkflow('/fake/.github/workflows/gem-build.yml', '2.6.0');

      expect(write).toHaveBeenCalledWith(
        '/fake/.github/workflows/gem-build.yml',
        'gem install type-on-strap --version "2.6.0" --source "https://rubygems.pkg.github.com/sylhare"\n'
      );
    });
  });

  describe('release()', () => {
    test('calls all updaters and logs the new version', () => {
      jest.spyOn(fs, 'readFileSync').mockReturnValue('{}' as any);
      jest.spyOn(fs, 'writeFileSync').mockImplementation(() => {});
      const log = jest.spyOn(console, 'log').mockImplementation(() => {});

      release('2.6.0', '/fake');

      expect(log).toHaveBeenCalledWith('Version bumped to 2.6.0');
      expect(fs.writeFileSync).toHaveBeenCalledWith(expect.stringContaining('type-on-strap.gemspec'), expect.any(String));
      expect(fs.writeFileSync).toHaveBeenCalledWith(expect.stringContaining('default.html'), expect.any(String));
      expect(fs.writeFileSync).toHaveBeenCalledWith(expect.stringContaining('package.json'), expect.any(String));
      expect(fs.writeFileSync).toHaveBeenCalledWith(expect.stringContaining('package-lock.json'), expect.any(String));
      expect(fs.writeFileSync).toHaveBeenCalledWith(expect.stringContaining('gem-build.yml'), expect.any(String));
    });
  });
});

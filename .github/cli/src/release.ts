import * as fs from 'node:fs';
import * as path from 'node:path';
import { updateJsonFile, updateVersionInFile } from './utils/fs';
import { logger } from './utils/logger';

export type BumpType = 'major' | 'minor' | 'patch';

const SEMVER = /^(\d+)\.(\d+)\.(\d+)$/;

export function parseVersion(version: string): [number, number, number] {
  const match = version.match(SEMVER);
  if (!match) throw new Error(`Invalid version format: ${version}`);
  return [parseInt(match[1]), parseInt(match[2]), parseInt(match[3])];
}

export function bumpVersion(current: string, bump: BumpType | string): string {
  if (SEMVER.test(bump)) return bump;
  const [major, minor, patch] = parseVersion(current);
  switch (bump) {
    case 'major':
      return `${major + 1}.0.0`;
    case 'minor':
      return `${major}.${minor + 1}.0`;
    case 'patch':
    default:
      return `${major}.${minor}.${patch + 1}`;
  }
}

export function updateGemspec(filePath: string, newVersion: string): void {
  updateVersionInFile(filePath, /(spec\.version\s*=\s*")[^"]+(")/, `$1${newVersion}$2`);
}

export function updateDefaultHtml(filePath: string, newVersion: string): void {
  updateVersionInFile(filePath, /(Type on Strap jekyll theme v)\d+\.\d+\.\d+/, `$1${newVersion}`);
}

export function updatePackageJson(filePath: string, newVersion: string): void {
  updateJsonFile(filePath, (pkg: { version: string }) => {
    pkg.version = newVersion;
  });
}

export function updatePackageLockJson(filePath: string, newVersion: string): void {
  updateJsonFile(filePath, (lock: { version: string; packages?: Record<string, { version: string }> }) => {
    lock.version = newVersion;
    if (lock.packages?.['']) lock.packages[''].version = newVersion;
  });
}

export function release(newVersion: string, root: string): void {
  updateGemspec(path.join(root, 'type-on-strap.gemspec'), newVersion);
  updateDefaultHtml(path.join(root, '_layouts/default.html'), newVersion);
  updatePackageJson(path.join(root, 'package.json'), newVersion);
  updatePackageLockJson(path.join(root, 'package-lock.json'), newVersion);
  logger.success(`Version bumped to ${newVersion}`);
}

if (require.main === module) {
  const arg = process.argv[2] || 'patch';
  const gemspec = path.join(process.cwd(), 'type-on-strap.gemspec');
  const current = fs.readFileSync(gemspec, 'utf8').match(/spec\.version\s*=\s*"([^"]+)"/)?.[1];
  if (!current) {
    logger.error('Could not read current version from type-on-strap.gemspec');
    process.exit(1);
  }
  const newVersion = bumpVersion(current, arg);
  release(newVersion, process.cwd());
}

import path from 'node:path';
import { logger } from '../utils/logger';
import { readVendorVersion } from '../utils/fs';
import { ValidationResult } from './types';
import { validateFile, runAsMain } from './common';

const PROJECT_ROOT = path.resolve(__dirname, '../../../..');
const VENDOR_CONFIG = path.join(PROJECT_ROOT, 'vendor.config.json');

export async function validate(): Promise<ValidationResult> {
  const version = readVendorVersion(VENDOR_CONFIG, 'masonry');
  const failures: string[] = [];

  logger.header(`Masonry Validation (v${version})`);

  const ok = await validateFile(
    'masonry.pkgd.min.js',
    path.join(PROJECT_ROOT, 'assets/js/vendor/masonry.pkgd.min.js'),
    `https://unpkg.com/masonry-layout@${version}/dist/masonry.pkgd.min.js`
  );
  if (!ok) failures.push('masonry.pkgd.min.js');

  return { passed: failures.length === 0, failures };
}

runAsMain(module, 'Masonry', validate);

import path from 'node:path';
import { logger } from '../utils/logger';
import { readVendorVersion } from '../utils/fs';
import { ValidationResult } from './types';
import { validateFile, runAsMain } from './common';

const PROJECT_ROOT = path.resolve(__dirname, '../../../..');
const VENDOR_CONFIG = path.join(PROJECT_ROOT, 'vendor.config.json');

export async function validate(): Promise<ValidationResult> {
  const version = readVendorVersion(VENDOR_CONFIG, 'imagesLoaded');
  const failures: string[] = [];

  logger.header(`imagesLoaded Validation (v${version})`);

  const ok = await validateFile(
    'imagesloaded.pkgd.min.js',
    path.join(PROJECT_ROOT, 'assets/js/vendor/imagesloaded.pkgd.min.js'),
    `https://unpkg.com/imagesloaded@${version}/imagesloaded.pkgd.min.js`
  );
  if (!ok) failures.push('imagesloaded.pkgd.min.js');

  return { passed: failures.length === 0, failures };
}

runAsMain(module, 'imagesLoaded', validate);

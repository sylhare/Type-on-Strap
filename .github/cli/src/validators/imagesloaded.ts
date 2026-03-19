import path from 'node:path';
import { ValidationResult } from './types';
import { validateSingleVendorFile, runAsMain } from './common';

const PROJECT_ROOT = path.resolve(__dirname, '../../../..');
const VENDOR_CONFIG = path.join(PROJECT_ROOT, 'vendor.config.json');

export async function validate(): Promise<ValidationResult> {
  return validateSingleVendorFile(
    VENDOR_CONFIG, 'imagesLoaded', 'imagesLoaded',
    'imagesloaded.pkgd.min.js',
    path.join(PROJECT_ROOT, 'assets/js/vendor/imagesloaded.pkgd.min.js'),
    v => `https://unpkg.com/imagesloaded@${v}/imagesloaded.pkgd.min.js`
  );
}

runAsMain(module, 'imagesLoaded', validate);

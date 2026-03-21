import * as path from 'node:path';
import { PROJECT_ROOT, ValidationResult, VENDOR_CONFIG } from '../types';
import { runAsMain, validateSingleVendorFile } from './common';

export async function validate(): Promise<ValidationResult> {
  return validateSingleVendorFile(
    VENDOR_CONFIG, 'imagesLoaded', 'imagesLoaded',
    'imagesloaded.pkgd.min.js',
    path.join(PROJECT_ROOT, 'assets/js/vendor/imagesloaded.pkgd.min.js'),
    v => `https://unpkg.com/imagesloaded@${v}/imagesloaded.pkgd.min.js`
  );
}

runAsMain(module, 'imagesLoaded', validate);

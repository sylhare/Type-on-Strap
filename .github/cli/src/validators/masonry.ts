import * as path from 'node:path';
import { PROJECT_ROOT, ValidationResult, VENDOR_CONFIG } from '../types';
import { runAsMain, validateSingleVendorFile } from './common';

export async function validate(): Promise<ValidationResult> {
  return validateSingleVendorFile(
    VENDOR_CONFIG, 'masonry', 'Masonry',
    'masonry.pkgd.min.js',
    path.join(PROJECT_ROOT, 'assets/js/vendor/masonry.pkgd.min.js'),
    v => `https://unpkg.com/masonry-layout@${v}/dist/masonry.pkgd.min.js`
  );
}

runAsMain(module, 'Masonry', validate);

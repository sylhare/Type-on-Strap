import path from 'node:path';
import { ValidationResult } from './types';
import { validateSingleVendorFile, runAsMain } from './common';

const PROJECT_ROOT = path.resolve(__dirname, '../../../..');
const VENDOR_CONFIG = path.join(PROJECT_ROOT, 'vendor.config.json');

export async function validate(): Promise<ValidationResult> {
  return validateSingleVendorFile(
    VENDOR_CONFIG, 'mermaid', 'Mermaid',
    'mermaid.min.js',
    path.join(PROJECT_ROOT, 'assets/js/vendor/mermaid.min.js'),
    v => `https://cdn.jsdelivr.net/npm/mermaid@${v}/dist/mermaid.min.js`
  );
}

runAsMain(module, 'Mermaid', validate);

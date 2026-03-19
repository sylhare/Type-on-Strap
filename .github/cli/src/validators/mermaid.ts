import * as path from 'node:path';
import { PROJECT_ROOT, ValidationResult, VENDOR_CONFIG } from '../types';
import { runAsMain, validateSingleVendorFile } from './common';

export async function validate(): Promise<ValidationResult> {
  return validateSingleVendorFile(
    VENDOR_CONFIG, 'mermaid', 'Mermaid',
    'mermaid.min.js',
    path.join(PROJECT_ROOT, 'assets/js/vendor/mermaid.min.js'),
    v => `https://cdn.jsdelivr.net/npm/mermaid@${v}/dist/mermaid.min.js`
  );
}

runAsMain(module, 'Mermaid', validate);

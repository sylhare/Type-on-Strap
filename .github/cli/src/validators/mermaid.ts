import path from 'node:path';
import { logger } from '../utils/logger';
import { readVendorVersion } from '../utils/fs';
import { ValidationResult } from './types';
import { validateFile, runAsMain } from './common';

const PROJECT_ROOT = path.resolve(__dirname, '../../../..');
const VENDOR_CONFIG = path.join(PROJECT_ROOT, 'vendor.config.json');

export async function validate(): Promise<ValidationResult> {
  const version = readVendorVersion(VENDOR_CONFIG, 'mermaid');
  const failures: string[] = [];

  logger.header(`Mermaid Validation (v${version})`);

  const ok = await validateFile(
    'mermaid.min.js',
    path.join(PROJECT_ROOT, 'assets/js/vendor/mermaid.min.js'),
    `https://cdn.jsdelivr.net/npm/mermaid@${version}/dist/mermaid.min.js`
  );
  if (!ok) failures.push('mermaid.min.js');

  return { passed: failures.length === 0, failures };
}

runAsMain(module, 'Mermaid', validate);

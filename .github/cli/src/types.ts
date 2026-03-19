import * as path from 'node:path';

export const PROJECT_ROOT = path.resolve(__dirname, '../../..');
export const VENDOR_CONFIG = path.join(PROJECT_ROOT, 'vendor.config.json');
export const HEAD_LIQUID = path.join(PROJECT_ROOT, '_includes/default/head.liquid');

export interface ValidationResult {
  passed: boolean;
  failures: string[];
}

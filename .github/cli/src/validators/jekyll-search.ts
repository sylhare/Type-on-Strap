import * as path from 'node:path';
import { PROJECT_ROOT, ValidationResult, VENDOR_CONFIG } from '../types';
import { runAsMain, validateSingleVendorFile } from './common';

const GITHUB_REPO = 'sylhare/Simple-Jekyll-Search';
const VENDOR_JS = path.join(PROJECT_ROOT, 'assets/js/vendor/simple-jekyll-search.min.js');

export async function validate(): Promise<ValidationResult> {
  return validateSingleVendorFile(
    VENDOR_CONFIG, 'simpleJekyllSearch', 'Simple-Jekyll-Search',
    'simple-jekyll-search.min.js',
    VENDOR_JS,
    v => `https://raw.githubusercontent.com/${GITHUB_REPO}/v${v}/dest/simple-jekyll-search.min.js`
  );
}

runAsMain(module, 'Simple-Jekyll-Search', validate);

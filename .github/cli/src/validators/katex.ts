import * as fs from 'node:fs';
import * as path from 'node:path';
import { fetchBuffer } from '../utils/http';
import { logger } from '../utils/logger';
import { readVendorVersion } from '../utils/fs';
import { PROJECT_ROOT, ValidationResult, VENDOR_CONFIG } from '../types';
import { runAsMain, validateFile } from './common';

async function validateCssVersion(scssFile: string, expectedVersion: string): Promise<boolean> {
  logger.info('Validating KaTeX CSS/SCSS version...');
  if (!fs.existsSync(scssFile)) {
    logger.error(`SCSS file not found: ${scssFile}`);
    return false;
  }
  const [localContent, cdnCssBuffer] = await Promise.all([
    Promise.resolve(fs.readFileSync(scssFile, 'utf8')),
    fetchBuffer(`https://cdn.jsdelivr.net/npm/katex@${expectedVersion}/dist/katex.css`),
  ]);
  const cdnCss = cdnCssBuffer.toString('utf8');

  const localMatch = localContent.match(/content: *"([0-9.]+)"/);
  const localVersion = localMatch?.[1];
  const cdnMatch = cdnCss.match(/content: *"([0-9.]+)"/);
  const cdnVersion = cdnMatch?.[1];

  logger.info(`  Local SCSS version:  ${localVersion}`);
  logger.info(`  CDN CSS version:     ${cdnVersion}`);
  logger.info(`  Expected version:    ${expectedVersion}`);

  if (localVersion === expectedVersion && cdnVersion === expectedVersion) {
    logger.success('CSS version matches!');
    return true;
  } else if (localVersion !== expectedVersion) {
    logger.error(`Local SCSS version mismatch! Expected ${expectedVersion}, found ${localVersion}`);
    return false;
  } else {
    logger.warn(`CDN version mismatch (CDN may have updated)`);
    return false;
  }
}

export async function validate(): Promise<ValidationResult> {
  const version = readVendorVersion(VENDOR_CONFIG, 'katex');
  const failures: string[] = [];

  logger.header(`KaTeX Validation (v${version})`);

  const jsMain = await validateFile(
    'KaTeX main library',
    path.join(PROJECT_ROOT, 'assets/js/vendor/katex.min.js'),
    `https://cdn.jsdelivr.net/npm/katex@${version}/dist/katex.min.js`
  );
  if (!jsMain) failures.push('katex.min.js');

  const jsAutoRender = await validateFile(
    'KaTeX auto-render',
    path.join(PROJECT_ROOT, 'assets/js/vendor/katex.auto-render.min.js'),
    `https://cdn.jsdelivr.net/npm/katex@${version}/dist/contrib/auto-render.min.js`
  );
  if (!jsAutoRender) failures.push('katex.auto-render.min.js');

  const scssVersion = await validateCssVersion(
    path.join(PROJECT_ROOT, '_sass/external/katex/katex.scss'),
    version
  );
  if (!scssVersion) failures.push('katex.scss version');

  return { passed: failures.length === 0, failures };
}

runAsMain(module, 'KaTeX', validate);

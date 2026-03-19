import fs from 'node:fs';
import path from 'node:path';
import { sha256File, sha256Buffer } from '../utils/hash';
import { fetchBuffer } from '../utils/http';
import { logger } from '../utils/logger';
import { ValidationResult } from './types';

const PROJECT_ROOT = path.resolve(__dirname, '../../../..');
const VENDOR_CONFIG = path.join(PROJECT_ROOT, 'vendor.config.json');

function getVersion(): string {
  const config = JSON.parse(fs.readFileSync(VENDOR_CONFIG, 'utf8')) as { katex: { version: string } };
  return config.katex.version;
}

async function validateFile(name: string, localPath: string, cdnUrl: string): Promise<boolean> {
  logger.info(`Validating ${name}...`);
  if (!fs.existsSync(localPath)) {
    logger.error(`Local file not found: ${localPath}`);
    return false;
  }
  const localHash = await sha256File(localPath);
  const remoteHash = sha256Buffer(await fetchBuffer(cdnUrl));
  logger.info(`  Local SHA256:  ${localHash}`);
  logger.info(`  Remote SHA256: ${remoteHash}`);
  if (localHash === remoteHash) {
    logger.success('Files match!');
    return true;
  } else {
    logger.error('Files DO NOT match!');
    return false;
  }
}

async function validateCssVersion(scssFile: string, expectedVersion: string): Promise<boolean> {
  logger.info('Validating KaTeX CSS/SCSS version...');
  if (!fs.existsSync(scssFile)) {
    logger.error(`SCSS file not found: ${scssFile}`);
    return false;
  }
  const localContent = fs.readFileSync(scssFile, 'utf8');
  const localMatch = localContent.match(/content: *"([0-9.]+)"/);
  const localVersion = localMatch?.[1];

  const cdnCss = (await fetchBuffer(`https://cdn.jsdelivr.net/npm/katex@${expectedVersion}/dist/katex.css`)).toString('utf8');
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
  const version = getVersion();
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

if (require.main === module) {
  validate().then(({ passed, failures }) => {
    if (passed) {
      logger.success('KaTeX validation passed!');
    } else {
      logger.error(`KaTeX validation failed! (${failures.join(', ')})`);
      process.exit(1);
    }
  }).catch(err => {
    logger.error((err as Error).message);
    process.exit(1);
  });
}

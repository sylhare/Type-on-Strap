import fs from 'node:fs';
import path from 'node:path';
import { sha256File, sha256Buffer } from '../utils/hash';
import { fetchBuffer } from '../utils/http';
import { logger } from '../utils/logger';
import { ValidationResult } from './types';

const PROJECT_ROOT = path.resolve(__dirname, '../../../..');
const VENDOR_CONFIG = path.join(PROJECT_ROOT, 'vendor.config.json');

function getVersion(): string {
  const config = JSON.parse(fs.readFileSync(VENDOR_CONFIG, 'utf8')) as { imagesLoaded: { version: string } };
  return config.imagesLoaded.version;
}

export async function validate(): Promise<ValidationResult> {
  const version = getVersion();
  const failures: string[] = [];

  logger.header(`imagesLoaded Validation (v${version})`);

  const localPath = path.join(PROJECT_ROOT, 'assets/js/vendor/imagesloaded.pkgd.min.js');
  const cdnUrl = `https://unpkg.com/imagesloaded@${version}/imagesloaded.pkgd.min.js`;

  logger.info('Validating imagesloaded.pkgd.min.js...');

  if (!fs.existsSync(localPath)) {
    logger.error(`Local file not found: ${localPath}`);
    failures.push('imagesloaded.pkgd.min.js');
    return { passed: false, failures };
  }

  const localHash = await sha256File(localPath);
  const remoteHash = sha256Buffer(await fetchBuffer(cdnUrl));

  logger.info(`  Local SHA256:  ${localHash}`);
  logger.info(`  Remote SHA256: ${remoteHash}`);

  if (localHash === remoteHash) {
    logger.success('Files match!');
  } else {
    logger.error('Files DO NOT match!');
    failures.push('imagesloaded.pkgd.min.js');
  }

  return { passed: failures.length === 0, failures };
}

if (require.main === module) {
  validate().then(({ passed, failures }) => {
    if (passed) {
      logger.success('imagesLoaded validation passed!');
    } else {
      logger.error(`imagesLoaded validation failed! (${failures.join(', ')})`);
      process.exit(1);
    }
  }).catch(err => {
    logger.error((err as Error).message);
    process.exit(1);
  });
}

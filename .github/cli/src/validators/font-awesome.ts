import fs from 'node:fs';
import path from 'node:path';
import { sha256File, sha256Buffer } from '../utils/hash';
import { fetchBuffer } from '../utils/http';
import { logger } from '../utils/logger';

const PROJECT_ROOT = path.resolve(__dirname, '../../../..');
const VENDOR_CONFIG = path.join(PROJECT_ROOT, 'vendor.config.json');

function getVersion(): string {
  const config = JSON.parse(fs.readFileSync(VENDOR_CONFIG, 'utf8')) as { fontAwesome: { version: string } };
  return config.fontAwesome.version;
}

export interface ValidationResult {
  passed: boolean;
  failures: string[];
}

interface FileEntry {
  name: string;
  localPath: string;
  remoteUrl: string;
}

async function validateFile(entry: FileEntry): Promise<string | null> {
  logger.info(`Validating ${entry.name}...`);
  if (!fs.existsSync(entry.localPath)) {
    logger.error(`Local file not found: ${entry.localPath}`);
    return entry.name;
  }
  const localHash = await sha256File(entry.localPath);
  const remoteHash = sha256Buffer(await fetchBuffer(entry.remoteUrl));
  logger.info(`  Local SHA256:   ${localHash}`);
  logger.info(`  Remote SHA256:  ${remoteHash}`);
  if (localHash === remoteHash) {
    logger.success('Files match!');
    return null;
  } else {
    logger.error('Files DO NOT match!');
    return entry.name;
  }
}

export async function validate(): Promise<ValidationResult> {
  const version = getVersion();
  const CDN_BASE = `https://cdnjs.cloudflare.com/ajax/libs/font-awesome/${version}/webfonts`;
  const GITHUB_SCSS_BASE = `https://raw.githubusercontent.com/FortAwesome/Font-Awesome/${version}/scss`;

  logger.header(`Font Awesome Validation (v${version})`);

  const fontFiles: FileEntry[] = ['fa-brands-400', 'fa-regular-400', 'fa-solid-900'].flatMap(font =>
    ['woff2', 'ttf'].map(ext => ({
      name: `${font}.${ext}`,
      localPath: path.join(PROJECT_ROOT, `assets/fonts/font-awesome/${font}.${ext}`),
      remoteUrl: `${CDN_BASE}/${font}.${ext}`,
    }))
  );

  const scssNames = [
    '_animated.scss', '_bordered-pulled.scss', '_core.scss', '_fixed-width.scss',
    '_functions.scss', '_icons.scss', '_list.scss', '_mixins.scss',
    '_rotated-flipped.scss', '_screen-reader.scss', '_shims.scss', '_sizing.scss',
    '_stacked.scss', '_variables.scss', 'brands.scss', 'fontawesome.scss',
    'regular.scss', 'solid.scss', 'v4-shims.scss',
  ];

  const scssFiles: FileEntry[] = scssNames.map(scss => ({
    name: scss,
    localPath: path.join(PROJECT_ROOT, `_sass/external/font-awesome/${scss}`),
    remoteUrl: `${GITHUB_SCSS_BASE}/${scss}`,
  }));

  logger.info('-- Font files --\n');
  const fontResults = await Promise.allSettled(fontFiles.map(validateFile));

  logger.info('\n-- SCSS files --\n');
  const scssResults = await Promise.allSettled(scssFiles.map(validateFile));

  const failures: string[] = [];
  for (const result of [...fontResults, ...scssResults]) {
    if (result.status === 'fulfilled' && result.value !== null) {
      failures.push(result.value);
    } else if (result.status === 'rejected') {
      failures.push('unknown (error during validation)');
    }
  }

  return { passed: failures.length === 0, failures };
}

if (require.main === module) {
  validate().then(({ passed, failures }) => {
    if (passed) {
      logger.success('Font Awesome validation passed!');
    } else {
      logger.error(`Font Awesome validation failed! (${failures.length} check(s) failed)`);
      process.exit(1);
    }
  }).catch(err => {
    logger.error((err as Error).message);
    process.exit(1);
  });
}

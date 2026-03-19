import path from 'node:path';
import { logger } from '../utils/logger';
import { readVendorVersion } from '../utils/fs';
import { ValidationResult } from './types';
import { validateFile, runAsMain } from './common';

const PROJECT_ROOT = path.resolve(__dirname, '../../../..');
const VENDOR_CONFIG = path.join(PROJECT_ROOT, 'vendor.config.json');

interface FileEntry {
  name: string;
  localPath: string;
  remoteUrl: string;
}

export async function validate(): Promise<ValidationResult> {
  const version = readVendorVersion(VENDOR_CONFIG, 'fontAwesome');
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

  const [fontResults, scssResults] = await Promise.all([
    Promise.allSettled(fontFiles.map(({ name, localPath, remoteUrl }) =>
      validateFile(name, localPath, remoteUrl).then(ok => ok ? null : name)
    )),
    Promise.allSettled(scssFiles.map(({ name, localPath, remoteUrl }) =>
      validateFile(name, localPath, remoteUrl).then(ok => ok ? null : name)
    )),
  ]);

  const failures: string[] = [];
  for (const result of [...fontResults, ...scssResults]) {
    if (result.status === 'fulfilled' && result.value !== null) {
      failures.push(result.value);
    } else if (result.status === 'rejected') {
      logger.error((result.reason as Error).message);
      failures.push('unknown (error during validation)');
    }
  }

  return { passed: failures.length === 0, failures };
}

runAsMain(module, 'Font Awesome', validate);

import * as path from 'node:path';
import { downloadFile, fetchJson } from '../utils/http';
import { updateVendorConfig, updateVersionInFile } from '../utils/fs';
import { logger } from '../utils/logger';
import { HEAD_LIQUID, PROJECT_ROOT, VENDOR_CONFIG } from '../types';

const VENDOR_JS = path.join(PROJECT_ROOT, 'assets/js/vendor/mermaid.min.js');

export async function fetchLatestVersion(): Promise<string> {
  logger.info('No version specified, querying npm registry for latest...');
  const data = await fetchJson<{ version: string }>('https://registry.npmjs.org/mermaid/latest');
  return data.version;
}

export async function resolveVersion(arg?: string): Promise<string> {
  if (arg) return arg;
  return fetchLatestVersion();
}

export async function updateMermaid(version: string): Promise<void> {
  logger.info(`\nUpdating Mermaid to v${version}...\n`);

  const cdnUrl = `https://cdn.jsdelivr.net/npm/mermaid@${version}/dist/mermaid.min.js`;
  logger.info(`Downloading ${cdnUrl}...`);
  await downloadFile(cdnUrl, VENDOR_JS);
  logger.info('  assets/js/vendor/mermaid.min.js');

  logger.info('\nUpdating version strings...');
  updateVendorConfig(VENDOR_CONFIG, 'mermaid', version);
  logger.info(`  vendor.config.json → mermaid.version="${version}"`);

  updateVersionInFile(
    HEAD_LIQUID,
    /<!-- Mermaid [\d.]+ -->/,
    `<!-- Mermaid ${version} -->`
  );
  logger.info(`  _includes/default/head.liquid → <!-- Mermaid ${version} -->`);

  logger.success('Mermaid update complete!');
  logger.info('   Run: npm run validate:mermaid');
}

if (require.main === module) {
  resolveVersion(process.argv[2])
    .then(version => updateMermaid(version))
    .catch(err => {
      logger.error((err as Error).message);
      process.exit(1);
    });
}

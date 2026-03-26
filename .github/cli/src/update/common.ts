import { downloadFile } from '../utils/http';
import { updateVendorConfig, updateVersionInFile } from '../utils/fs';
import { logger } from '../utils/logger';
import { HEAD_LIQUID, PROJECT_ROOT, VENDOR_CONFIG } from '../types';

export async function updateSingleVendorFile(
  displayName: string,
  vendorKey: string,
  vendorFilePath: string,
  validateScript: string,
  version: string,
  downloadUrl: (version: string) => string,
): Promise<void> {
  logger.info(`\nUpdating ${displayName} to v${version}...\n`);

  const url = downloadUrl(version);
  logger.info(`Downloading ${url}...`);
  await downloadFile(url, vendorFilePath);
  logger.info(`  ${vendorFilePath.replace(PROJECT_ROOT + '/', '')}`);

  logger.info('\nUpdating version strings...');
  updateVendorConfig(VENDOR_CONFIG, vendorKey, version);
  logger.info(`  vendor.config.json → ${vendorKey}.version="${version}"`);

  updateVersionInFile(
    HEAD_LIQUID,
    new RegExp(`<!-- ${displayName} [\\d.]+ -->`),
    `<!-- ${displayName} ${version} -->`
  );
  logger.info(`  _includes/default/head.liquid → <!-- ${displayName} ${version} -->`);

  logger.success(`${displayName} update complete!`);
  logger.info(`   Run: ${validateScript}`);
}

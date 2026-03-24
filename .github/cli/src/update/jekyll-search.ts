import * as path from 'node:path';
import { downloadFile, fetchJson } from '../utils/http';
import { updateVendorConfig, updateVersionInFile } from '../utils/fs';
import { logger } from '../utils/logger';
import { HEAD_LIQUID, PROJECT_ROOT, VENDOR_CONFIG } from '../types';

const GITHUB_REPO = 'sylhare/Simple-Jekyll-Search';
const VENDOR_JS = path.join(PROJECT_ROOT, 'assets/js/vendor/simple-jekyll-search.min.js');

interface GithubRelease {
  tag_name: string;
}

export async function fetchLatestVersion(): Promise<string> {
  logger.info('No version specified, querying GitHub for latest release...');
  const data = await fetchJson<GithubRelease>(`https://api.github.com/repos/${GITHUB_REPO}/releases/latest`);
  return data.tag_name.replace(/^v/, '');
}

export async function resolveVersion(arg?: string): Promise<string> {
  if (arg) return arg.replace(/^v/, '');
  return fetchLatestVersion();
}

export async function updateJekyllSearch(version: string): Promise<void> {
  logger.info(`\nUpdating Simple-Jekyll-Search to v${version}...\n`);

  const downloadUrl = `https://raw.githubusercontent.com/${GITHUB_REPO}/v${version}/dest/simple-jekyll-search.min.js`;
  logger.info(`Downloading ${downloadUrl}...`);
  await downloadFile(downloadUrl, VENDOR_JS);
  logger.info('  assets/js/vendor/simple-jekyll-search.min.js');

  logger.info('\nUpdating version strings...');
  updateVendorConfig(VENDOR_CONFIG, 'simpleJekyllSearch', version);
  logger.info(`  vendor.config.json → simpleJekyllSearch.version="${version}"`);

  updateVersionInFile(
    HEAD_LIQUID,
    /<!-- Simple Jekyll Search [\d.]+ -->/,
    `<!-- Simple Jekyll Search ${version} -->`
  );
  logger.info(`  _includes/default/head.liquid → <!-- Simple Jekyll Search ${version} -->`);

  logger.success('Simple-Jekyll-Search update complete!');
  logger.info('   Run: npm run validate:search');
}

if (require.main === module) {
  resolveVersion(process.argv[2])
    .then(version => updateJekyllSearch(version))
    .catch(err => {
      logger.error((err as Error).message);
      process.exit(1);
    });
}

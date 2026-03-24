import * as path from 'node:path';
import { fetchJson } from '../utils/http';
import { logger } from '../utils/logger';
import { PROJECT_ROOT } from '../types';
import { updateSingleVendorFile } from './common';

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
  await updateSingleVendorFile(
    'Simple-Jekyll-Search', 'simpleJekyllSearch', VENDOR_JS, 'npm run validate:search',
    version,
    v => `https://raw.githubusercontent.com/${GITHUB_REPO}/v${v}/dest/simple-jekyll-search.min.js`
  );
}

if (require.main === module) {
  resolveVersion(process.argv[2])
    .then(version => updateJekyllSearch(version))
    .catch(err => {
      logger.error((err as Error).message);
      process.exit(1);
    });
}

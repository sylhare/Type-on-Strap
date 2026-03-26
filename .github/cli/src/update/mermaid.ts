import * as path from 'node:path';
import { fetchJson } from '../utils/http';
import { logger } from '../utils/logger';
import { PROJECT_ROOT } from '../types';
import { updateSingleVendorFile } from './common';

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
  await updateSingleVendorFile(
    'Mermaid', 'mermaid', VENDOR_JS, 'npm run validate:mermaid',
    version,
    v => `https://cdn.jsdelivr.net/npm/mermaid@${v}/dist/mermaid.min.js`
  );
}

if (require.main === module) {
  resolveVersion(process.argv[2])
    .then(version => updateMermaid(version))
    .catch(err => {
      logger.error((err as Error).message);
      process.exit(1);
    });
}

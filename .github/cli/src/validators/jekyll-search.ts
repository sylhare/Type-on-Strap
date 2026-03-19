import fs from 'node:fs';
import path from 'node:path';
import { sha256File, sha256Buffer } from '../utils/hash';
import { fetchBuffer, fetchJson } from '../utils/http';
import { logger } from '../utils/logger';
import { ValidationResult } from './types';

const PROJECT_ROOT = path.resolve(__dirname, '../../../..');
const GITHUB_REPO = 'sylhare/Simple-Jekyll-Search';
const LOCAL_FILE = path.join(PROJECT_ROOT, 'assets/js/vendor/simple-jekyll-search.min.js');

interface GithubRelease {
  tag_name: string;
}

export async function validate(): Promise<ValidationResult> {
  const failures: string[] = [];

  logger.header('Simple-Jekyll-Search Validation');
  logger.info(`Repository: https://github.com/${GITHUB_REPO}\n`);

  if (!fs.existsSync(LOCAL_FILE)) {
    logger.error(`Local file not found: ${LOCAL_FILE}`);
    return { passed: false, failures: ['simple-jekyll-search.min.js'] };
  }

  const localContent = fs.readFileSync(LOCAL_FILE, 'utf8');
  const localLines = localContent.split('\n').slice(0, 3).join('\n');
  const localVersionMatch = localLines.match(/v(\d+\.\d+\.\d+)/);
  const localVersion = localVersionMatch ? `v${localVersionMatch[1]}` : null;

  logger.info(`Local version:  ${localVersion ?? 'Unable to detect'}`);

  let latestVersion: string | null = null;
  let downloadUrl: string | null = null;

  try {
    const release = await fetchJson<GithubRelease>(
      `https://api.github.com/repos/${GITHUB_REPO}/releases/latest`
    );
    latestVersion = release.tag_name;
    downloadUrl = `https://raw.githubusercontent.com/${GITHUB_REPO}/${latestVersion}/dest/simple-jekyll-search.min.js`;
    logger.info(`Latest release: ${latestVersion}`);
  } catch {
    logger.info('Latest release: Unable to fetch (network unavailable or rate limited)');
    logger.warn('Could not verify against latest release');
    logger.info(`  Check manually: https://github.com/${GITHUB_REPO}/releases`);
    return { passed: failures.length === 0, failures };
  }

  if (localVersion === latestVersion) {
    logger.info('✓ Version matches');
  } else {
    logger.warn('Version mismatch detected');
    logger.info(`   Local:  ${localVersion}`);
    logger.info(`   Latest: ${latestVersion}`);
    failures.push('version mismatch');
  }

  if (downloadUrl) {
    logger.info('\nContent check:');
    try {
      const remoteBuf = await fetchBuffer(downloadUrl);
      const localHash = await sha256File(LOCAL_FILE);
      const remoteHash = sha256Buffer(remoteBuf);

      if (localHash === remoteHash) {
        logger.info('✓ Content matches official release');
      } else {
        logger.error('Content differs from official release');
        logger.info(`   Local hash:  ${localHash.slice(0, 16)}...`);
        logger.info(`   Remote hash: ${remoteHash.slice(0, 16)}...`);
        failures.push('content mismatch');
      }
    } catch {
      logger.warn('Unable to download for comparison');
    }
  }

  return { passed: failures.length === 0, failures };
}

if (require.main === module) {
  validate().then(({ passed }) => {
    if (passed) {
      logger.success('Simple-Jekyll-Search validation passed!');
    } else {
      logger.error('Validation failed');
      process.exit(1);
    }
  }).catch(err => {
    logger.error((err as Error).message);
    process.exit(1);
  });
}

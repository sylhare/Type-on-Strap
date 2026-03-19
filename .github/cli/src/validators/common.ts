import { sha256File, sha256Buffer } from '../utils/hash';
import { fetchBuffer } from '../utils/http';
import { logger } from '../utils/logger';
import { readVendorVersion } from '../utils/fs';
import { ValidationResult } from '../types';

export async function validateFile(name: string, localPath: string, cdnUrl: string): Promise<boolean> {
  logger.info(`Validating ${name}...`);
  try {
    const [localHash, remoteHash] = await Promise.all([
      sha256File(localPath),
      fetchBuffer(cdnUrl).then(sha256Buffer),
    ]);
    logger.info(`  Local SHA256:  ${localHash}`);
    logger.info(`  Remote SHA256: ${remoteHash}`);
    if (localHash === remoteHash) {
      logger.success('Files match!');
      return true;
    } else {
      logger.error('Files DO NOT match!');
      return false;
    }
  } catch (err) {
    logger.error((err as Error).message);
    return false;
  }
}

export async function validateSingleVendorFile(
  vendorConfigPath: string,
  vendorKey: string,
  displayName: string,
  fileName: string,
  localPath: string,
  remoteUrl: (version: string) => string
): Promise<ValidationResult> {
  const version = readVendorVersion(vendorConfigPath, vendorKey);
  logger.header(`${displayName} Validation (v${version})`);
  const ok = await validateFile(fileName, localPath, remoteUrl(version));
  return { passed: ok, failures: ok ? [] : [fileName] };
}

export function runAsMain(callerModule: NodeModule, name: string, validate: () => Promise<ValidationResult>): void {
  if (require.main === callerModule) {
    validate().then(({ passed, failures }) => {
      if (passed) {
        logger.success(`${name} validation passed!`);
      } else {
        logger.error(`${name} validation failed! (${failures.join(', ')})`);
        process.exit(1);
      }
    }).catch(err => {
      logger.error((err as Error).message);
      process.exit(1);
    });
  }
}

import * as fs from 'node:fs';
import * as path from 'node:path';
import { globSync } from 'glob';
import * as esbuild from 'esbuild';
import { logger } from './utils/logger';

export function getJsPartials(dir: string): string[] {
  return globSync(path.join(dir, '*.js')).sort();
}

export function concatFiles(filePaths: string[]): string {
  return filePaths.map(f => fs.readFileSync(f, 'utf8')).join('\n');
}

export async function buildJs(partialFiles: string[], outFile: string): Promise<void> {
  const code = concatFiles(partialFiles);
  const result = await esbuild.transform(code, { minify: true });
  fs.writeFileSync(outFile, result.code);
}

if (require.main === module) {
  const cwd = process.cwd();

  async function runJs(): Promise<void> {
    const partialsDir = path.join(cwd, 'assets/js/partials');
    const partials = getJsPartials(partialsDir);
    const commentsFile = path.join(cwd, 'assets/js/comments-lazy-load.js');
    await Promise.all([
      buildJs(partials, path.join(cwd, 'assets/js/main.min.js'))
        .then(() => logger.info('Built assets/js/main.min.js')),
      buildJs([commentsFile], path.join(cwd, 'assets/js/comments-lazy-load.min.js'))
        .then(() => logger.info('Built assets/js/comments-lazy-load.min.js')),
    ]);
  }

  runJs().catch(err => {
    logger.error(String(err));
    process.exit(1);
  });
}

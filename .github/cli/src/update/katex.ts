import fs from 'node:fs';
import path from 'node:path';
import os from 'node:os';
import { execSync } from 'child_process';
import { globSync } from 'glob';
import { updateVersionInFile, updateVendorConfig, readVendorVersion } from '../utils/fs';
import { logger } from '../utils/logger';
import { PROJECT_ROOT, VENDOR_CONFIG, HEAD_LIQUID } from '../types';

const FONTS_DIR = path.join(PROJECT_ROOT, 'assets/fonts/katex');
const KATEX_SCSS = path.join(PROJECT_ROOT, '_sass/external/katex/katex.scss');
const KATEX_ENTRY = path.join(PROJECT_ROOT, '_sass/external/_katex.scss');

export function resolveVersion(arg?: string): string {
  if (arg) return arg;
  return readVendorVersion(VENDOR_CONFIG, 'katex');
}

function copyFile(src: string, dest: string): void {
  fs.mkdirSync(path.dirname(dest), { recursive: true });
  fs.copyFileSync(src, dest);
}

/**
 * Transforms `dist/katex.css` into a SCSS partial.
 * TTF entries are stripped manually: `USE_TTF=false` only applies to `katex.min.css`.
 */
export function generateScss(distCss: string): string {
  const header = '$katex-font-path: "../../assets/fonts/katex" !default;\n\n/* stylelint-disable font-family-no-missing-generic-family-keyword */\n';

  let css = distCss
    .replace(/\/\*\s*stylelint-disable[^*]*\*\/\n?/g, '')
    .replace(/,\s*url\(["']?[^)"']*\.ttf["']?\)\s*format\(["']truetype["']\)/g, '')
    .replace(/url\(fonts\/([^)]+)\)/g, 'url("#{$katex-font-path}/$1")')
    .replace(/url\(["']fonts\/([^"')]+)["']\)/g, 'url("#{$katex-font-path}/$1")');

  return header + css;
}

export async function updateKatex(version: string): Promise<void> {
  logger.info(`\nUpdating KaTeX to v${version}...\n`);

  const tmp = fs.mkdtempSync(path.join(os.tmpdir(), 'katex-'));
  logger.info(`Temp dir: ${tmp}`);

  try {
    logger.info(`Cloning KaTeX v${version}...`);
    execSync(`git clone --depth 1 --branch v${version} https://github.com/KaTeX/KaTeX.git ${tmp}`, { stdio: 'inherit' });

    logger.info('\nInstalling dependencies...');
    execSync('yarn install --frozen-lockfile', { cwd: tmp, stdio: 'inherit' });

    logger.info('\nBuilding KaTeX (USE_TTF=false)...');
    execSync('yarn build', {
      cwd: tmp,
      stdio: 'inherit',
      env: { ...process.env, USE_TTF: 'false' },
    });

    logger.info('\nCopying JS files...');
    copyFile(
      path.join(tmp, 'dist/katex.min.js'),
      path.join(PROJECT_ROOT, 'assets/js/vendor/katex.min.js')
    );
    logger.info('  assets/js/vendor/katex.min.js');
    copyFile(
      path.join(tmp, 'dist/contrib/auto-render.min.js'),
      path.join(PROJECT_ROOT, 'assets/js/vendor/katex.auto-render.min.js')
    );
    logger.info('  assets/js/vendor/katex.auto-render.min.js');

    logger.info('\nUpdating fonts...');
    const ttfFiles = globSync(path.join(FONTS_DIR, '*.ttf'));
    for (const f of ttfFiles) {
      fs.unlinkSync(f);
      logger.info(`  Deleted ${path.basename(f)}`);
    }
    const newFonts = globSync(path.join(tmp, 'dist/fonts/*.{woff2,woff}'));
    for (const src of newFonts) {
      const dest = path.join(FONTS_DIR, path.basename(src));
      fs.copyFileSync(src, dest);
    }
    logger.info(`  Copied ${newFonts.length} font files (woff2 + woff)`);

    logger.info('\nGenerating SCSS...');
    const distCss = fs.readFileSync(path.join(tmp, 'dist/katex.css'), 'utf8');
    const scss = generateScss(distCss);
    fs.writeFileSync(KATEX_SCSS, scss);
    logger.info('  _sass/external/katex/katex.scss');

    logger.info('\nUpdating version strings...');
    updateVendorConfig(VENDOR_CONFIG, 'katex', version);
    logger.info(`  vendor.config.json → katex.version="${version}"`);

    updateVersionInFile(
      KATEX_ENTRY,
      /KaTeX v[\d.]+/,
      `KaTeX v${version}`
    );
    logger.info(`  _sass/external/_katex.scss → KaTeX v${version}`);

    updateVersionInFile(
      HEAD_LIQUID,
      /<!-- KaTeX [\d.]+ -->/,
      `<!-- KaTeX ${version} -->`
    );
    logger.info(`  _includes/default/head.liquid → KaTeX ${version}`);

  } finally {
    fs.rmSync(tmp, { recursive: true, force: true });
    logger.info('\nTemp directory cleaned up.');
  }

  logger.success('KaTeX update complete!');
  logger.info('   Run: npm run validate:katex');
}

if (require.main === module) {
  let version: string;
  try {
    version = resolveVersion(process.argv[2]);
  } catch (err) {
    logger.error((err as Error).message);
    process.exit(1);
  }

  updateKatex(version).catch(err => {
    logger.error(String(err));
    process.exit(1);
  });
}

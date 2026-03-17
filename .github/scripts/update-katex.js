'use strict';
const fs = require('fs');
const path = require('path');
const os = require('os');
const { execSync } = require('child_process');
const { globSync } = require('glob');

const PROJECT_ROOT = path.resolve(__dirname, '../..');
const VALIDATE_SCRIPT = path.join(PROJECT_ROOT, '.github/scripts/validate-katex.sh');
const FONTS_DIR = path.join(PROJECT_ROOT, 'assets/fonts/katex');
const KATEX_SCSS = path.join(PROJECT_ROOT, '_sass/external/katex/katex.scss');
const KATEX_ENTRY = path.join(PROJECT_ROOT, '_sass/external/_katex.scss');

function resolveVersion(arg) {
  if (arg) return arg;
  const script = fs.readFileSync(VALIDATE_SCRIPT, 'utf8');
  const match = script.match(/^KATEX_VERSION="([^"]+)"/m);
  if (!match) throw new Error('Could not find KATEX_VERSION in validate-katex.sh');
  return match[1];
}

function copyFile(src, dest) {
  fs.mkdirSync(path.dirname(dest), { recursive: true });
  fs.copyFileSync(src, dest);
}

/**
 * Transforms `dist/katex.css` into a SCSS partial.
 * We use the unminified file because it is multi-line and safe to process with regex.
 * TTF entries are stripped manually: `USE_TTF=false` only applies to `katex.min.css`.
 */
function generateScss(distCss) {
  const header = '$katex-font-path: "../../assets/fonts/katex" !default;\n\n/* stylelint-disable font-family-no-missing-generic-family-keyword */\n';

  let css = distCss
    // Strip any existing stylelint-disable comment (avoid duplication)
    .replace(/\/\*\s*stylelint-disable[^*]*\*\/\n?/g, '')
    // Remove TTF src entries — handles both quoted and unquoted URLs
    .replace(/,\s*url\(["']?[^)"']*\.ttf["']?\)\s*format\(["']truetype["']\)/g, '')
    // Replace unquoted font path: url(fonts/file.ext) → url("#{$katex-font-path}/file.ext")
    .replace(/url\(fonts\/([^)]+)\)/g, 'url("#{$katex-font-path}/$1")')
    // Replace quoted font path: url("fonts/file.ext") → url("#{$katex-font-path}/file.ext")
    .replace(/url\(["']fonts\/([^"')]+)["']\)/g, 'url("#{$katex-font-path}/$1")');

  return header + css;
}

function updateVersionInFile(filePath, pattern, replacement) {
  const content = fs.readFileSync(filePath, 'utf8');
  if (!pattern.test(content)) throw new Error(`Pattern not found in ${filePath}`);
  const updated = content.replace(pattern, replacement);
  if (updated !== content) fs.writeFileSync(filePath, updated);
}

async function updateKatex(version) {
  console.log(`\nUpdating KaTeX to v${version}...\n`);

  // 1. Create temp dir and clone
  const tmp = fs.mkdtempSync(path.join(os.tmpdir(), 'katex-'));
  console.log(`Temp dir: ${tmp}`);

  try {
    // 2. Clone
    console.log(`Cloning KaTeX v${version}...`);
    execSync(`git clone --depth 1 --branch v${version} https://github.com/KaTeX/KaTeX.git ${tmp}`, { stdio: 'inherit' });

    // 3. Install
    console.log('\nInstalling dependencies...');
    execSync('yarn install --frozen-lockfile', { cwd: tmp, stdio: 'inherit' });

    // 4. Build without TTF
    console.log('\nBuilding KaTeX (USE_TTF=false)...');
    execSync('yarn build', {
      cwd: tmp,
      stdio: 'inherit',
      env: { ...process.env, USE_TTF: 'false' },
    });

    // 5. Copy JS files
    console.log('\nCopying JS files...');
    copyFile(
      path.join(tmp, 'dist/katex.min.js'),
      path.join(PROJECT_ROOT, 'assets/js/vendor/katex.min.js')
    );
    console.log('  assets/js/vendor/katex.min.js');
    copyFile(
      path.join(tmp, 'dist/contrib/auto-render.min.js'),
      path.join(PROJECT_ROOT, 'assets/js/vendor/katex.auto-render.min.js')
    );
    console.log('  assets/js/vendor/katex.auto-render.min.js');

    // 6. Update fonts — delete existing TTF, copy woff2 + woff
    console.log('\nUpdating fonts...');
    const ttfFiles = globSync(path.join(FONTS_DIR, '*.ttf'));
    for (const f of ttfFiles) {
      fs.unlinkSync(f);
      console.log(`  Deleted ${path.basename(f)}`);
    }
    const newFonts = globSync(path.join(tmp, 'dist/fonts/*.{woff2,woff}'));
    for (const src of newFonts) {
      const dest = path.join(FONTS_DIR, path.basename(src));
      fs.copyFileSync(src, dest);
    }
    console.log(`  Copied ${newFonts.length} font files (woff2 + woff)`);

    // 7. Generate SCSS from built CSS
    console.log('\nGenerating SCSS...');
    const distCss = fs.readFileSync(path.join(tmp, 'dist/katex.css'), 'utf8');
    const scss = generateScss(distCss);
    fs.writeFileSync(KATEX_SCSS, scss);
    console.log('  _sass/external/katex/katex.scss');

    // 8. Update version strings
    console.log('\nUpdating version strings...');
    updateVersionInFile(
      VALIDATE_SCRIPT,
      /^KATEX_VERSION="[^"]+"/m,
      `KATEX_VERSION="${version}"`
    );
    console.log(`  validate-katex.sh → KATEX_VERSION="${version}"`);

    updateVersionInFile(
      KATEX_ENTRY,
      /KaTeX v[\d.]+/,
      `KaTeX v${version}`
    );
    console.log(`  _sass/external/_katex.scss → KaTeX v${version}`);

    const HEAD_LIQUID = path.join(PROJECT_ROOT, '_includes/default/head.liquid');
    updateVersionInFile(
      HEAD_LIQUID,
      /<!-- KaTeX [\d.]+ -->/,
      `<!-- KaTeX ${version} -->`
    );
    console.log(`  _includes/default/head.liquid → KaTeX ${version}`);

  } finally {
    // 9. Cleanup
    fs.rmSync(tmp, { recursive: true, force: true });
    console.log('\nTemp directory cleaned up.');
  }

  console.log('\n✅ KaTeX update complete!');
  console.log('   Run: bash .github/scripts/validate-katex.sh');
}

module.exports = { updateKatex, resolveVersion, generateScss };

if (require.main === module) {
  const arg = process.argv[2];
  let version;
  try {
    version = resolveVersion(arg);
  } catch (err) {
    console.error(err.message);
    process.exit(1);
  }

  updateKatex(version).catch(err => {
    console.error(err);
    process.exit(1);
  });
}

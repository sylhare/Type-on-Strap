'use strict';
const fs = require('fs');
const path = require('path');
const https = require('https');

const PROJECT_ROOT = path.resolve(__dirname, '../..');
const VENDOR_JS = path.join(PROJECT_ROOT, 'assets/js/vendor/mermaid.min.js');
const VALIDATE_SCRIPT = path.join(PROJECT_ROOT, '.github/scripts/validate-mermaid.sh');
const HEAD_LIQUID = path.join(PROJECT_ROOT, '_includes/default/head.liquid');

function fetchLatestVersion() {
  return new Promise((resolve, reject) => {
    https.get('https://registry.npmjs.org/mermaid/latest', res => {
      let data = '';
      res.on('data', chunk => { data += chunk; });
      res.on('end', () => {
        try {
          resolve(JSON.parse(data).version);
        } catch (err) {
          reject(new Error('Failed to parse npm registry response'));
        }
      });
    }).on('error', reject);
  });
}

function downloadFile(url, dest) {
  return new Promise((resolve, reject) => {
    const file = fs.createWriteStream(dest);
    https.get(url, res => {
      if (res.statusCode === 301 || res.statusCode === 302) {
        file.close();
        return downloadFile(res.headers.location, dest).then(resolve).catch(reject);
      }
      if (res.statusCode !== 200) {
        file.close();
        return reject(new Error(`HTTP ${res.statusCode} for ${url}`));
      }
      res.pipe(file);
      file.on('finish', () => file.close(resolve));
    }).on('error', err => {
      fs.unlink(dest, () => {});
      reject(err);
    });
  });
}

function updateVersionInFile(filePath, pattern, replacement) {
  const content = fs.readFileSync(filePath, 'utf8');
  if (!pattern.test(content)) throw new Error(`Pattern not found in ${filePath}`);
  const updated = content.replace(pattern, replacement);
  if (updated !== content) fs.writeFileSync(filePath, updated);
}

async function resolveVersion(arg) {
  if (arg) return arg;
  console.log('No version specified, querying npm registry for latest...');
  return fetchLatestVersion();
}

async function updateMermaid(version) {
  console.log(`\nUpdating Mermaid to v${version}...\n`);

  // 1. Download mermaid.min.js from CDN
  const cdnUrl = `https://cdn.jsdelivr.net/npm/mermaid@${version}/dist/mermaid.min.js`;
  console.log(`Downloading ${cdnUrl}...`);
  await downloadFile(cdnUrl, VENDOR_JS);
  console.log('  assets/js/vendor/mermaid.min.js');

  // 2. Update version in validate-mermaid.sh
  console.log('\nUpdating version strings...');
  updateVersionInFile(
    VALIDATE_SCRIPT,
    /^MERMAID_VERSION="[^"]+"/m,
    `MERMAID_VERSION="${version}"`
  );
  console.log(`  validate-mermaid.sh → MERMAID_VERSION="${version}"`);

  // 3. Update version comment in head.liquid
  updateVersionInFile(
    HEAD_LIQUID,
    /<!-- Mermaid [\d.]+ -->/,
    `<!-- Mermaid ${version} -->`
  );
  console.log(`  _includes/default/head.liquid → <!-- Mermaid ${version} -->`);

  console.log('\n✅ Mermaid update complete!');
  console.log('   Run: bash .github/scripts/validate-mermaid.sh');
}

module.exports = { updateMermaid, resolveVersion };

if (require.main === module) {
  resolveVersion(process.argv[2])
    .then(version => updateMermaid(version))
    .catch(err => {
      console.error(err.message);
      process.exit(1);
    });
}

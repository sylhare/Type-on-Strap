import { logger } from '../utils/logger';
import { validate as validateKatex } from './katex';
import { validate as validateMermaid } from './mermaid';
import { validate as validateFontAwesome } from './font-awesome';
import { validate as validateMasonry } from './masonry';
import { validate as validateImagesLoaded } from './imagesloaded';
import { validate as validateJekyllSearch } from './jekyll-search';

interface ValidatorEntry {
  name: string;
  fn: () => Promise<{ passed: boolean; failures: string[] }>;
}

const VALIDATORS: ValidatorEntry[] = [
  { name: 'Font Awesome',         fn: validateFontAwesome },
  { name: 'KaTeX',                fn: validateKatex },
  { name: 'Mermaid',              fn: validateMermaid },
  { name: 'Masonry',              fn: validateMasonry },
  { name: 'imagesLoaded',         fn: validateImagesLoaded },
  { name: 'Simple-Jekyll-Search', fn: validateJekyllSearch },
];

async function main(): Promise<void> {
  console.log('\n==================================================');
  console.log('Vendor Dependencies Validation');
  console.log('==================================================\n');

  let passed = 0;
  let failed = 0;

  for (const { name, fn } of VALIDATORS) {
    logger.info(`▶ ${name}`);
    try {
      const result = await fn();
      if (result.passed) {
        passed++;
      } else {
        failed++;
        logger.error(`${name} failed: ${result.failures.join(', ')}`);
      }
    } catch (err) {
      failed++;
      logger.error(`${name} threw an error: ${(err as Error).message}`);
    }
    console.log('');
  }

  console.log('==================================================');
  console.log('Summary');
  console.log('==================================================');
  console.log(`Passed: ${passed}`);
  console.log(`Failed: ${failed}`);
  console.log('');

  if (failed === 0) {
    logger.success('All vendor dependencies validated successfully!');
  } else {
    logger.error(`${failed} validation(s) failed. Please re-download the failed files.`);
    process.exit(1);
  }
}

if (require.main === module) {
  main().catch(err => {
    logger.error((err as Error).message);
    process.exit(1);
  });
}

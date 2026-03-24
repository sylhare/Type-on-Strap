import { test } from '@playwright/test';
import { takeThemeScreenshots } from './support/helpers';

function defineScreenshotTests(): void {
  test.setTimeout(60_000);

  test('markdown article', async ({ page }) => {
    await page.goto('/syntax/2014/11/28/markdown-and-html');
    await page.waitForLoadState('load');
    await takeThemeScreenshots(page, 'markdown-article');
  });

  test('mermaid article', async ({ page }) => {
    await page.goto('/2016/12/03/Mermaid');
    await page.waitForSelector('.mermaid svg, .language-mermaid svg', { timeout: 15_000 });
    await takeThemeScreenshots(page, 'mermaid-article');
  });

  test('katex article', async ({ page }) => {
    await page.goto('/2016/12/04/Katex');
    await page.waitForSelector('.katex', { timeout: 10_000 });
    await takeThemeScreenshots(page, 'katex-article');
  });

  test('fontawesome article', async ({ page }) => {
    await page.goto('/2016/12/05/Font-Awesome');
    await page.waitForLoadState('load');
    await takeThemeScreenshots(page, 'fontawesome-article');
  });

  test('bootstrap article', async ({ page }) => {
    await page.goto('/2017/09/17/Use-Bootstrap');
    await page.waitForLoadState('load');
    await takeThemeScreenshots(page, 'bootstrap-article');
  });

  test('search with query', async ({ page }) => {
    await page.goto('/search');
    await page.waitForSelector('#search-input');
    await page.fill('#search-input', 'test');
    await page.waitForTimeout(500);
    await takeThemeScreenshots(page, 'search-query');
  });

  test('gallery page', async ({ page }) => {
    await page.goto('/gallery');
    await page.waitForLoadState('load');
    await takeThemeScreenshots(page, 'gallery');
  });

  test('tags page', async ({ page }) => {
    await page.goto('/tags');
    await page.waitForLoadState('load');
    await takeThemeScreenshots(page, 'tags');
  });

  test('portfolio page', async ({ page }) => {
    await page.goto('/portfolio');
    await page.waitForLoadState('load');
    await takeThemeScreenshots(page, 'portfolio');
  });
}

test.describe('Page Screenshots @visual', () => {
  defineScreenshotTests();
});

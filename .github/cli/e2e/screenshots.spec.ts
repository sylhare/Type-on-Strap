import { expect, test, type Page } from '@playwright/test';

/**
 * Force a theme by writing to localStorage and setting the data-theme attribute directly.
 */
async function setTheme(page: Page, theme: 'light' | 'dark'): Promise<void> {
  await page.evaluate((t) => {
    localStorage.setItem('theme', t);
    document.documentElement.setAttribute('data-theme', t);
  }, theme);
  await page.waitForTimeout(300);
}

/**
 * Capture light-mode and dark-mode screenshots for the current page state.
 * Screenshots are stored in a `screenshots.spec.ts-snapshots/` directory and
 * compared on subsequent runs (visual regression). Run with `--update-snapshots`
 * to create or refresh the baseline.
 */
async function takeThemeScreenshots(page: Page, name: string): Promise<void> {
  await setTheme(page, 'light');
  await expect(page).toHaveScreenshot(`${name}-light.png`, { fullPage: true });

  await setTheme(page, 'dark');
  await expect(page).toHaveScreenshot(`${name}-dark.png`, { fullPage: true });
}

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

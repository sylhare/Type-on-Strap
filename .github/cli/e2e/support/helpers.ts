import { expect, test, type Page } from '@playwright/test';

/**
 * Helper to open mobile menu - forces menu visibility for mobile tests
 */
export async function openMobileMenu(page: Page): Promise<void> {
  await page.waitForLoadState('domcontentloaded');
  await page.waitForTimeout(300);

  const menu = page.locator('nav ul');
  await menu.evaluate(el => {
    el.classList.remove('hide');
    (el as HTMLElement).style.opacity = '1';
    (el as HTMLElement).style.fontSize = '';
  });

  await page.waitForTimeout(200);
}

/**
 * Helper function to check if theme toggle is available on the page.
 */
export async function hasThemeToggle(page: Page): Promise<boolean> {
  const themeToggle = page.locator('#theme-toggle');
  return await themeToggle.count() > 0;
}

/**
 * Force a theme by writing to localStorage and setting the data-theme attribute directly.
 */
export async function setTheme(page: Page, theme: 'light' | 'dark'): Promise<void> {
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
export async function takeThemeScreenshots(page: Page, name: string): Promise<void> {
  await setTheme(page, 'light');
  await expect(page).toHaveScreenshot(`${name}-light.png`, { fullPage: true });

  await setTheme(page, 'dark');
  await expect(page).toHaveScreenshot(`${name}-dark.png`, { fullPage: true });
}

/**
 * Helper function to check if gallery has images and skip test if not
 */
export async function checkGalleryHasImages(page: Page) {
  const images = page.locator('.gallery img, .gallery-item img, img[class*="gallery"]');
  const count = await images.count();
  if (count === 0) {
    test.skip(true, 'Gallery has no images configured');
  }
  return { images, count };
}

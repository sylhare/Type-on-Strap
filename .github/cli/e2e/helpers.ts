import { test } from '@playwright/test';
import type { Page } from '@playwright/test';

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

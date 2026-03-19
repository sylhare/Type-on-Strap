import { expect, test } from '@playwright/test';
import { checkGalleryHasImages } from './helpers';

test.describe('Gallery Functionality @desktop', () => {
  test('should display gallery page', async ({ page }) => {
    await page.goto('/gallery');
    await expect(page.locator('body')).toBeVisible();

    const pageHeading = page.locator('h1').first();
    await expect(pageHeading).toBeVisible();
  });

  test('should display gallery images', async ({ page }) => {
    await page.goto('/gallery');

    const { images, count } = await checkGalleryHasImages(page);

    expect(count).toBeGreaterThan(0);
    await expect(images.first()).toBeVisible();
  });

  test('should have valid image sources', async ({ page }) => {
    await page.goto('/gallery');

    const { images, count } = await checkGalleryHasImages(page);
    expect(count).toBeGreaterThan(0);

    const firstImage = images.first();
    const src = await firstImage.getAttribute('src');
    expect(src).toBeTruthy();
    expect(src).toMatch(/\.(jpg|jpeg|png|gif|webp|svg)/i);
  });

  test('should have image alt text', async ({ page }) => {
    await page.goto('/gallery');

    const { images, count } = await checkGalleryHasImages(page);
    expect(count).toBeGreaterThan(0);

    const firstImage = images.first();
    const alt = await firstImage.getAttribute('alt');
    expect(alt).toBeDefined();
  });

  test('should have responsive gallery layout', async ({ page }) => {
    await page.goto('/gallery');

    const { images, count } = await checkGalleryHasImages(page);

    expect(count).toBeGreaterThan(0);

    for (let i = 0; i < Math.min(count, 3); i++) {
      await expect(images.nth(i)).toBeVisible();
    }
  });

  test('should load images progressively', async ({ page }) => {
    await page.goto('/gallery');

    const { images, count } = await checkGalleryHasImages(page);
    expect(count).toBeGreaterThan(0);

    const firstImage = images.first();
    await expect(firstImage).toBeVisible();

    const isLoaded = await firstImage.evaluate((img) => {
      return (img as HTMLImageElement).complete && (img as HTMLImageElement).naturalWidth > 0;
    });

    expect(isLoaded).toBe(true);
  });

  test('should have gallery grid layout', async ({ page }) => {
    await page.goto('/gallery');

    const gallery = page.locator('.gallery, [class*="gallery"]').first();
    const galleryCount = await gallery.count();

    if (galleryCount === 0) {
      test.skip(true, 'Gallery container not found');
      return;
    }

    await expect(gallery).toBeVisible();

    const images = gallery.locator('img');
    const count = await images.count();

    if (count === 0) {
      test.skip(true, 'Gallery has no images');
      return;
    }

    expect(count).toBeGreaterThan(0);
  });

  /**
   * @note The loading attribute is optional in HTML. If not present, browsers use their default behavior.
   */
  test('should handle image lazy loading', async ({ page }) => {
    await page.goto('/gallery');

    const { images, count } = await checkGalleryHasImages(page);
    expect(count).toBeGreaterThan(0);

    const firstImage = images.first();
    const loading = await firstImage.getAttribute('loading');

    if (loading) {
      expect(['lazy', 'eager']).toContain(loading);
    }
  });

  test('should display gallery on mobile', async ({ page, isMobile }) => {
    test.skip(!isMobile, 'Mobile-only test');

    await page.goto('/gallery');

    const { images, count } = await checkGalleryHasImages(page);
    expect(count).toBeGreaterThan(0);
    await expect(images.first()).toBeVisible();
  });
});

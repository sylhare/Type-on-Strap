import { test, expect } from '@playwright/test';

test.describe('Portfolio Functionality @desktop', () => {
  test('should display portfolio page', async ({ page }) => {
    await page.goto('/portfolio');
    await expect(page.locator('body')).toBeVisible();

    const pageHeading = page.locator('h1').first();
    await expect(pageHeading).toBeVisible();
  });

  test('should display portfolio items', async ({ page }) => {
    await page.goto('/portfolio');

    const items = page.locator('[data-testid="portfolio-item"]');
    const count = await items.count();
    expect(count).toBeGreaterThan(0);
  });

  test('should open individual portfolio item', async ({ page }) => {
    await page.goto('/portfolio');

    const firstItem = page.locator('[data-testid="portfolio-item-link"]').first();

    await expect(firstItem).toBeVisible();
    await firstItem.click();

    await expect(page.locator('body')).toBeVisible();
    expect(page.url()).toContain('portfolio');
  });

  test('should have portfolio item images', async ({ page }) => {
    await page.goto('/portfolio');

    const images = page.locator('[data-testid="portfolio-item-image"]');
    const count = await images.count();

    expect(count).toBeGreaterThan(0);

    const firstImage = images.first();

    const src = await firstImage.getAttribute('src');
    expect(src).toBeTruthy();
    expect(src).toMatch(/\.(jpg|jpeg|png|gif|webp|svg)/i);

    await firstImage.evaluate((img) => {
      return (img as HTMLImageElement).complete || new Promise((resolve) => {
        (img as HTMLImageElement).onload = resolve as () => void;
        (img as HTMLImageElement).onerror = resolve as () => void;
      });
    });
  });

  test('should have portfolio item titles', async ({ page }) => {
    await page.goto('/portfolio');

    const items = page.locator('[data-testid="portfolio-item"]');
    expect(await items.count()).toBeGreaterThan(0);

    const firstItem = items.first();
    await expect(firstItem).toBeVisible();

    const caption = firstItem.locator('[data-testid="portfolio-item-caption"]');
    const titleAttr = await caption.getAttribute('title');
    expect(titleAttr).toBeTruthy();
    expect(titleAttr!.length).toBeGreaterThan(0);
  });

  test('should have portfolio item content', async ({ page }) => {
    await page.goto('/portfolio');

    const items = page.locator('[data-testid="portfolio-item"]');
    const firstItem = items.first();

    await expect(firstItem).toBeVisible();

    await expect(firstItem.locator('[data-testid="portfolio-item-link"]')).toBeVisible();
    await expect(firstItem.locator('[data-testid="portfolio-item-image"]')).toBeVisible();
    await expect(firstItem.locator('[data-testid="portfolio-item-caption"]')).toBeVisible();

    const caption = firstItem.locator('[data-testid="portfolio-item-caption"]');
    const title = await caption.getAttribute('title');
    expect(title).toBeTruthy();
  });

  test('should have responsive portfolio grid', async ({ page }) => {
    await page.goto('/portfolio');

    const items = page.locator('[data-testid="portfolio-item"]');
    const count = await items.count();
    expect(count).toBeGreaterThan(0);

    for (let i = 0; i < Math.min(count, 3); i++) {
      await expect(items.nth(i)).toBeVisible();
    }
  });

  test('should navigate back to portfolio from detail page', async ({ page }) => {
    await page.goto('/portfolio');

    const firstItem = page.locator('[data-testid="portfolio-item-link"]').first();

    await expect(firstItem).toBeVisible();
    await firstItem.click();

    await page.goBack();

    await expect(page).toHaveURL(/portfolio/);
  });

  test('should have portfolio metadata', async ({ page }) => {
    await page.goto('/portfolio');

    const firstItem = page.locator('[data-testid="portfolio-item-link"]').first();

    await expect(firstItem).toBeVisible();
    await firstItem.click();

    await expect(page.locator('body')).toBeVisible();

    const content = page.locator('article, .post-content, main').first();
    await expect(content).toBeVisible();
  });
});

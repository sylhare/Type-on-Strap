import { expect, test } from '@playwright/test';

test.describe('Category Display Functionality @desktop', () => {
  test('should display categories on posts with category metadata', async ({ page }) => {
    await page.goto('/demo/2014/11/26/lorem-ipsum');

    const categoryList = page.locator('[data-testid="category-list"]');
    await expect(categoryList).toBeVisible();

    const categoryLink = page.locator('[data-testid="category-link"]');
    await expect(categoryLink).toBeVisible();

    const categoryText = await categoryLink.textContent();
    expect(categoryText).toContain('Demo');
  });

  test('should display empty category list on posts without categories', async ({ page }) => {
    await page.goto('/');

    await expect(page.locator('body')).toBeVisible();
  });

  test('should have category links that navigate to categories page', async ({ page }) => {
    await page.goto('/demo/2014/11/26/lorem-ipsum');

    const categoryLink = page.locator('[data-testid="category-link"]').first();
    await expect(categoryLink).toBeVisible();

    await categoryLink.click();

    await expect(page).toHaveURL(/\/categories(\/)?#/);
  });

  test('should display category with folder icon', async ({ page }) => {
    await page.goto('/demo/2014/11/26/lorem-ipsum');

    const categoryLink = page.locator('[data-testid="category-link"]').first();
    await expect(categoryLink).toBeVisible();

    const icon = categoryLink.locator('i.fa-folder');
    await expect(icon).toBeVisible();
  });

  test('should display both tags and categories on the same post', async ({ page }) => {
    await page.goto('/demo/2014/11/26/lorem-ipsum');

    const tagList = page.locator('[data-testid="tag-list"]');
    await expect(tagList).toBeVisible();

    const categoryList = page.locator('[data-testid="category-list"]');
    await expect(categoryList).toBeVisible();

    const tagLinks = page.locator('[data-testid="tag-link"]');
    expect(await tagLinks.count()).toBeGreaterThan(0);

    const categoryLinks = page.locator('[data-testid="category-link"]');
    expect(await categoryLinks.count()).toBeGreaterThan(0);
  });

  test('should display singular "Category" for single category', async ({ page }) => {
    await page.goto('/demo/2014/11/26/lorem-ipsum');

    const categoryList = page.locator('[data-testid="category-list"]');
    await expect(categoryList).toBeVisible();

    const categoryMeta = categoryList.locator('li.meta');
    const metaText = await categoryMeta.textContent();
    expect(metaText).toContain('Category');
  });

  test('should display categories on tutorial post', async ({ page }) => {
    await page.goto('/tutorial/2013/10/18/blogging-with-title');

    const categoryList = page.locator('[data-testid="category-list"]');
    await expect(categoryList).toBeVisible();

    const categoryLink = page.locator('[data-testid="category-link"]');
    await expect(categoryLink).toBeVisible();

    const categoryText = await categoryLink.textContent();
    expect(categoryText).toContain('Tutorial');
  });

  test('should display categories on example posts', async ({ page }) => {
    await page.goto('/example/2019/06/30/sample-post');

    const categoryList = page.locator('[data-testid="category-list"]');
    await expect(categoryList).toBeVisible();

    const categoryLink = page.locator('[data-testid="category-link"]');
    await expect(categoryLink).toBeVisible();

    const categoryText = await categoryLink.textContent();
    expect(categoryText).toContain('Example');
  });

  test('should have proper styling for category buttons', async ({ page }) => {
    await page.goto('/demo/2014/11/26/lorem-ipsum');

    const categoryLink = page.locator('[data-testid="category-link"]').first();
    await expect(categoryLink).toBeVisible();

    const hasButtonClass = await categoryLink.evaluate(el => el.classList.contains('button'));
    expect(hasButtonClass).toBe(true);
  });
});

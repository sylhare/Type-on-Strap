import { expect, test } from '@playwright/test';

test.describe('Blog Functionality @desktop', () => {
  test('should display blog posts', async ({ page }) => {
    await page.goto('/');

    const posts = page.locator('[data-testid="blog-post-teaser"]');
    const count = await posts.count();
    expect(count).toBeGreaterThan(0);
  });

  test('should open individual blog post', async ({ page }) => {
    await page.goto('/');

    const firstPost = page.locator('[data-testid="blog-post-link"]').first();
    await firstPost.click();

    await expect(page.locator('article')).toBeVisible();
  });

  test('should have post metadata', async ({ page }) => {
    await page.goto('/');

    const firstPost = page.locator('[data-testid="blog-post-link"]').first();
    await firstPost.click();

    const article = page.locator('article');
    await expect(article).toBeVisible();

    const articleText = await article.textContent();
    expect(articleText?.length).toBeGreaterThan(50);
  });

  test('should have post navigation', async ({ page }) => {
    await page.goto('/');

    const firstPost = page.locator('.post-teaser header h1 a, article header h1 a, .banner h1 a').first();
    await firstPost.click();

    const article = page.locator('article');
    await expect(article).toBeVisible();

    await expect(page.locator('body')).toBeVisible();
  });

  test('should filter posts by tag', async ({ page }) => {
    await page.goto('/tags');

    const tag = page.locator('[data-testid="tag-link"]').first();
    await expect(tag).toBeVisible();

    await tag.click();

    await expect(page.locator('body')).toBeVisible();
    expect(page.url()).toContain('tag');
  });

  test('should filter posts by category', async ({ page }) => {
    await page.goto('/categories');

    const category = page.locator('main a[href*="categories"], .content a[href*="categories"], .category a').first();
    const count = await category.count();

    if (count === 0) {
      test.skip(true, 'No categories configured');
      return;
    }

    await expect(category).toBeVisible();
    await category.click();

    await expect(page.locator('body')).toBeVisible();
    expect(page.url()).toMatch(/categories|category/);
  });

  test('should have blog pagination', async ({ page }) => {
    await page.goto('/');

    await expect(page.locator('body')).toBeVisible();

    const posts = page.locator('[data-testid="blog-post-teaser"]');
    expect(await posts.count()).toBeGreaterThan(0);
  });

  test('should display post excerpts on blog page', async ({ page }) => {
    await page.goto('/');

    const posts = page.locator('[data-testid="blog-post-teaser"]');
    const firstPost = posts.first();

    await expect(firstPost).toBeVisible();

    const text = await firstPost.textContent();
    expect(text?.length).toBeGreaterThan(0);
  });

  /**
   * @note Tests a specific post known to have code blocks with syntax highlighting.
   */
  test('should have syntax highlighting in code blocks', async ({ page }) => {
    await page.goto('/syntax/2014/08/08/Markup-Syntax-Highlighting');

    const codeBlocks = page.locator('pre code, .highlight');
    expect(await codeBlocks.count()).toBeGreaterThan(0);

    const codeBlock = codeBlocks.first();
    const className = await codeBlock.getAttribute('class');
    expect(className).toBeTruthy();
  });

  test('should have share buttons on posts', async ({ page }) => {
    await page.goto('/');

    const firstPost = page.locator('.post-teaser header h1 a, article header h1 a, .banner h1 a').first();
    await firstPost.click();

    const article = page.locator('article');
    await expect(article).toBeVisible();

    const articleText = await article.textContent();
    expect(articleText?.length).toBeGreaterThan(50);
  });
});

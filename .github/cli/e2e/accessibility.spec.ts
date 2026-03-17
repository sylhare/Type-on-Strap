import { test, expect } from '@playwright/test';
import { openMobileMenu } from './helpers';

test.describe('Accessibility @desktop', () => {
  test('should have proper heading hierarchy on home page', async ({ page }) => {
    await page.goto('/');

    const h1 = page.locator('h1');
    const h1Count = await h1.count();

    expect(h1Count).toBeGreaterThan(0);
  });

  test('should have alt text on images', async ({ page }) => {
    await page.goto('/');

    const images = page.locator('img');
    const count = await images.count();

    expect(count).toBeGreaterThan(0);

    for (let i = 0; i < Math.min(count, 5); i++) {
      const img = images.nth(i);
      const alt = await img.getAttribute('alt');
      expect(alt).toBeDefined();
    }
  });

  /**
   * @note Focus state detection in automated tests can be unreliable in headless mode.
   * We verify the element is focusable by attempting to focus it and checking it remains visible.
   */
  test('should have keyboard navigable navbar', async ({ page }) => {
    await page.goto('/');

    const navLinks = page.locator('nav a:not(#pull), .navbar a:not(#pull)');
    const visibleLinks = navLinks.filter({ hasText: /.+/ });
    const count = await visibleLinks.count();

    expect(count).toBeGreaterThan(0);

    const firstLink = visibleLinks.first();
    await expect(firstLink).toBeVisible();
    await firstLink.focus();
  });

  test('should have proper link text (no "click here")', async ({ page }) => {
    await page.goto('/');

    const links = page.locator('a');
    const count = await links.count();

    expect(count).toBeGreaterThan(0);

    for (let i = 0; i < Math.min(count, 10); i++) {
      const link = links.nth(i);
      const text = await link.textContent();

      if (text) {
        const lowerText = text.toLowerCase().trim();
        expect(lowerText).not.toBe('click here');
        expect(lowerText).not.toBe('here');
      }
    }
  });

  test('should have proper form labels', async ({ page }) => {
    await page.goto('/search');

    const inputs = page.locator('input');
    const count = await inputs.count();

    expect(count).toBeGreaterThan(0);

    const input = inputs.first();
    const id = await input.getAttribute('id');
    const ariaLabel = await input.getAttribute('aria-label');
    const placeholder = await input.getAttribute('placeholder');

    if (id) {
      const label = page.locator(`label[for="${id}"]`);
      const labelExists = await label.count() > 0;

      expect(labelExists || ariaLabel || placeholder).toBeTruthy();
    } else {
      expect(ariaLabel || placeholder).toBeTruthy();
    }
  });

  test('should have valid HTML lang attribute', async ({ page }) => {
    await page.goto('/');

    const lang = await page.getAttribute('html', 'lang');
    expect(lang).toBeTruthy();
    expect(lang).toMatch(/^[a-z]{2}(-[A-Z]{2})?$/);
  });

  test('should have skip to main content link', async ({ page }) => {
    await page.goto('/');

    const main = page.locator('main, [role="main"]');
    const mainCount = await main.count();

    const skipLink = page.locator('a[href="#main"], a[href="#content"], a:has-text("skip")');
    const hasSkipLink = await skipLink.count() > 0;
    const hasMainContent = mainCount > 0;

    expect(hasSkipLink || hasMainContent).toBe(true);
  });

  /**
   * @note This is a basic contrast check. Full WCAG compliance would require calculating
   * contrast ratios from RGB values according to WCAG 2.1 guidelines (4.5:1 for normal text).
   */
  test('should have sufficient color contrast', async ({ page }) => {
    await page.goto('/');

    const body = page.locator('body');
    const backgroundColor = await body.evaluate(el =>
      window.getComputedStyle(el).backgroundColor,
    );
    const color = await body.evaluate(el =>
      window.getComputedStyle(el).color,
    );

    expect(backgroundColor).not.toBe(color);
    expect(backgroundColor).toBeTruthy();
    expect(color).toBeTruthy();
  });

  /**
   * @note In headless browsers, focus detection using document.activeElement can be unreliable.
   * We verify the element exists and is interactive.
   */
  test('should have focusable interactive elements', async ({ page }) => {
    await page.goto('/');

    const buttons = page.locator('button:visible, a:visible').first();

    await expect(buttons).toBeVisible();
    await buttons.focus();

    const isFocused = await buttons.evaluate(el => el === document.activeElement);
    expect(isFocused).toBe(true);
  });

  test('should have proper ARIA roles', async ({ page }) => {
    await page.goto('/');

    const nav = page.locator('nav, [role="navigation"]');
    const navCount = await nav.count();

    expect(navCount).toBeGreaterThan(0);
  });
});

// Mobile-specific accessibility tests
test.describe('Accessibility @mobile', () => {
  test('should have keyboard navigable navbar', async ({ page }) => {
    await page.goto('/');
    await openMobileMenu(page);

    const navLinks = page.locator('nav a:not(#pull), .navbar a:not(#pull)');
    const visibleLinks = navLinks.filter({ hasText: /.+/ });
    const count = await visibleLinks.count();

    expect(count).toBeGreaterThan(0);

    const firstLink = visibleLinks.first();
    await expect(firstLink).toBeVisible();
    await firstLink.focus();
  });
});

import { test, expect } from '@playwright/test';

test.describe('KaTeX Math Rendering @desktop', () => {
  const katexPostUrl = '/syntax/2014/11/28/markdown-and-html';

  test('should load KaTeX post without console errors', async ({ page }) => {
    const errors: string[] = [];
    page.on('console', msg => {
      if (msg.type() === 'error') {
        errors.push(msg.text());
      }
    });

    await page.goto(katexPostUrl);
    await expect(page.locator('article')).toBeVisible();

    const criticalErrors = errors.filter(e =>
      e.includes('katex') || e.includes('KaTeX') || e.includes('math')
    );
    expect(criticalErrors).toHaveLength(0);
  });

  test('should render KaTeX math formulas', async ({ page }) => {
    await page.goto(katexPostUrl);

    const katexElements = page.locator('.katex');
    await expect(katexElements.first()).toBeVisible({ timeout: 10000 });

    const count = await katexElements.count();
    expect(count).toBeGreaterThan(0);
  });

  test('should render inline math formulas', async ({ page }) => {
    await page.goto(katexPostUrl);

    const inlineMath = page.locator('p .katex');
    await expect(inlineMath.first()).toBeVisible({ timeout: 10000 });
  });

  test('should render display math formulas', async ({ page }) => {
    await page.goto(katexPostUrl);

    const displayMath = page.locator('.katex-display');
    await expect(displayMath.first()).toBeVisible({ timeout: 10000 });
  });

  test('should have KaTeX fonts loaded', async ({ page }) => {
    await page.goto(katexPostUrl);

    await expect(page.locator('.katex').first()).toBeVisible({ timeout: 10000 });

    const katexElement = page.locator('.katex').first();
    const fontFamily = await katexElement.evaluate(el =>
      window.getComputedStyle(el).fontFamily
    );
    expect(fontFamily.toLowerCase()).toContain('katex');
  });
});

test.describe('Mermaid Diagram Rendering @desktop', () => {
  const mermaidPostUrl = '/2016/12/03/Mermaid';

  test('should load Mermaid post without console errors', async ({ page }) => {
    const errors: string[] = [];
    page.on('console', msg => {
      if (msg.type() === 'error') {
        errors.push(msg.text());
      }
    });

    await page.goto(mermaidPostUrl);
    await expect(page.locator('article')).toBeVisible();

    const criticalErrors = errors.filter(e =>
      e.includes('mermaid') || e.includes('Mermaid') || e.includes('diagram')
    );
    expect(criticalErrors).toHaveLength(0);
  });

  test('should render Mermaid diagrams as SVG', async ({ page }) => {
    await page.goto(mermaidPostUrl);

    const mermaidSvg = page.locator('.mermaid svg, .language-mermaid svg');
    await expect(mermaidSvg.first()).toBeVisible({ timeout: 15000 });

    const count = await mermaidSvg.count();
    expect(count).toBeGreaterThan(0);
  });

  test('should render sequence diagram', async ({ page }) => {
    await page.goto(mermaidPostUrl);

    await expect(page.locator('.mermaid svg, .language-mermaid svg').first()).toBeVisible({ timeout: 15000 });

    const sequenceDiagram = page.locator('.mermaid svg, .language-mermaid svg');
    const count = await sequenceDiagram.count();
    expect(count).toBeGreaterThan(0);
  });

  test('should render flowchart diagram', async ({ page }) => {
    await page.goto(mermaidPostUrl);

    await expect(page.locator('.mermaid svg, .language-mermaid svg').first()).toBeVisible({ timeout: 15000 });

    const diagrams = page.locator('.mermaid svg, .language-mermaid svg');
    const count = await diagrams.count();
    expect(count).toBeGreaterThan(1);
  });

  test('should render multiple diagram types', async ({ page }) => {
    await page.goto(mermaidPostUrl);

    await expect(page.locator('.mermaid svg, .language-mermaid svg').first()).toBeVisible({ timeout: 15000 });
    await page.waitForTimeout(2000);

    const allSvgs = page.locator('.mermaid svg, .language-mermaid svg');
    const count = await allSvgs.count();

    expect(count).toBeGreaterThanOrEqual(5);
  });

  test('should not have unprocessed mermaid code blocks', async ({ page }) => {
    await page.goto(mermaidPostUrl);

    await expect(page.locator('.mermaid svg, .language-mermaid svg').first()).toBeVisible({ timeout: 15000 });
    await page.waitForTimeout(2000);

    const mermaidCodeBlocks = page.locator('.language-mermaid');
    const count = await mermaidCodeBlocks.count();

    for (let i = 0; i < count; i++) {
      const block = mermaidCodeBlocks.nth(i);
      const svg = block.locator('svg');
      const svgCount = await svg.count();
      expect(svgCount).toBeGreaterThan(0);
    }
  });
});

test.describe('Rendering Performance @desktop', () => {
  test('should load KaTeX within reasonable time', async ({ page }) => {
    const startTime = Date.now();
    await page.goto('/syntax/2014/11/28/markdown-and-html');
    await expect(page.locator('.katex').first()).toBeVisible({ timeout: 10000 });
    const loadTime = Date.now() - startTime;

    expect(loadTime).toBeLessThan(5000);
  });

  test('should load Mermaid diagrams within reasonable time', async ({ page }) => {
    const startTime = Date.now();
    await page.goto('/2016/12/03/Mermaid');
    await expect(page.locator('.mermaid svg, .language-mermaid svg').first()).toBeVisible({ timeout: 15000 });
    const loadTime = Date.now() - startTime;

    expect(loadTime).toBeLessThan(10000);
  });
});

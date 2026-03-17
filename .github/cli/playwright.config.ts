import { defineConfig, devices } from '@playwright/test';

/**
 * Playwright configuration for Type-on-Strap theme e2e tests
 * @see https://playwright.dev/docs/test-configuration
 */
export default defineConfig({
  testDir: './e2e',
  timeout: 30 * 1000,  // Maximum time one test can run
  fullyParallel: true,
  forbidOnly: !!process.env['CI'],
  retries: process.env['CI'] ? 2 : 0,
  workers: process.env['CI'] ? 1 : undefined,
  reporter: [
    ['html', { outputFolder: 'playwright-report' }],
    ['list'],
  ],
  use: {
    // Base URL for tests
    // Note: Site is served from _site/ root, but links contain /Type-on-Strap/ prefix
    baseURL: process.env['BASE_URL'] ?? 'http://localhost:4000',
    trace: 'on-first-retry',
    screenshot: 'only-on-failure',
    video: 'retain-on-failure',
  },
  projects: [
    {
      name: 'chromium',
      use: { ...devices['Desktop Chrome'] },
      grep: /@desktop/,
      grepInvert: /@mobile/,
    },
    {
      name: 'firefox',
      use: { ...devices['Desktop Firefox'] },
      grep: /@desktop/,
      grepInvert: /@mobile/,
    },
    {
      name: 'webkit',
      use: { ...devices['Desktop Safari'] },
      grep: /@desktop/,
      grepInvert: /@mobile/,
    },
    // Mobile tests
    {
      name: 'mobile-chrome',
      use: { ...devices['Pixel 5'] },
      grep: /@mobile/,
      grepInvert: /@desktop/,
    },
    {
      name: 'mobile-safari',
      use: { ...devices['iPhone 12'] },
      grep: /@mobile/,
      grepInvert: /@desktop/,
    },
  ],
  webServer: {
    command: 'cd ../.. && bundle exec jekyll build --baseurl "" --quiet && cd .github/cli && npm run server',
    url: 'http://localhost:4000',
    reuseExistingServer: !process.env['CI'],
    timeout: 120 * 1000,
  },
});

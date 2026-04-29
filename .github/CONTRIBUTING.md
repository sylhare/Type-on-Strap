# Contributing

## Quick reference

**Stack:** 
- Blog: Jekyll (Ruby)
- CLI scripts to help manage the theme: Node.js (TypeScript + Playwright e2e tests)

### Known gotchas

- **Sass `@import` deprecation**: The theme intentionally uses `@import` (not `@use`). This is blocked by GitHub Pages /
  jekyll-sass-converter 1.5.2 using Ruby Sass 3.7.4. Do not migrate to `@use` yet. `quiet_deps: true` in `_config.yml`
  silences Font Awesome's own warnings.
- **Jest virtual mocks**: External packages not in `.github/cli/node_modules` must be mocked with
  `jest.mock('pkg', factory, { virtual: true })`.
- **CSS variables for theming**: Light/dark mode and skin colors are all CSS custom properties defined in
  `_sass/base/_variables.scss`. The `data-theme` attribute is set by JavaScript in `head.liquid` on page load.
- **Remote theme vs gem**: `_config.yml` has both `remote_theme` and `theme` (commented). For local dev the gem is used;
  for GitHub Pages the remote_theme is used.
- **Jekyll exclude list**: `package.json`, `package-lock.json`, and `.github/` are excluded from Jekyll build. Any new
  tooling files at root should be added to the `exclude:` list in `_config.yml`.
- **Gem includes**: The gemspec in `type-on-strap.gemspec` only includes `assets/(js|css|fonts|data)/`,
  `_(includes|layouts|sass)/`, and `_data/(icons_builder|language).yml`. Content files (`_posts/`, `_portfolio/`,
  `pages/`) are NOT shipped in the gem.

## Prerequisites

- **Node.js** >= 18 and **npm** — for the CLI tooling and tests
- **Ruby** and **Bundler** — for running Jekyll locally

Install Node dependencies from the project root:

```bash
npm install
```

Install Ruby dependencies:

```bash
bundle install
```

## Blog commands

Using bundler instead of Jekyll directly.

```bash
bundle exec jekyll serve        # Serve locally at http://localhost:4000
bundle exec jekyll build        # Build to _site/
```

## CLI scripts

All scripts live in [`.github/cli/src/`](.github/cli/src/) and are run from the project root.

### Build

Minify JS:

```bash
npm run build        # minify JS
npm run build:js     # minify JS
```

### Images

```bash
npm run compress         # compress images in place
npm run thumbnails       # create thumbnails for feature-img/ only
npm run thumbnails-all   # create thumbnails for all images
npm run webp             # convert images to WebP
```

### Create a post

```bash
npm run post 'title of the post'
```

Creates `_posts/YYYY-MM-DD-title-of-the-post.md` with default frontmatter. Does nothing if the file already exists.

### Vendor dependencies

Validate that local vendor files match their upstream sources:

```bash
npm run validate            # all vendors
npm run validate:katex
npm run validate:mermaid
npm run validate:fa
npm run validate:masonry
npm run validate:search
```

Update a vendor to its latest version:

```bash
npm run update:katex [version]
npm run update:mermaid [version]
```

## Testing

### Unit tests

```bash
npm test
```

### Typecheck

```bash
npm run typecheck
```

### End-to-end tests

Requires a built site (`bundle exec jekyll build`) and Playwright browsers (`npm run playwright:install`):

```bash
npm run test:e2e
```

#### Visual regression tests

Visual tests are tagged `@visual` and run against 5 browser/device projects (desktop Chrome, Firefox, Safari, mobile
Chrome, mobile Safari).

Run only the visual tests:

```bash
npm run test:e2e -- --grep @visual
```

Or target a single browser:

```bash
npm run test:e2e -- --project=visual-chromium
```

Update the baseline screenshots when visual changes are intentional:

```bash
npm run test:e2e -- --grep @visual --update-snapshots
```

## Git hooks

A pre-commit hook is provided that checks for non-staged assets before committing.

Enable it:

```bash
ln .github/hooks/pre-commit .git/hooks/pre-commit
```

Bypass it when needed:

```bash
git commit -n
```

## Pull requests

- PRs must remain compatible with [GitHub Pages](https://github.com/github/pages-gem)
- Include a screenshot if the change affects the layout or visual appearance
- Run `npm test` and `npm run typecheck` before submitting

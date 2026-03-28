# Contributing

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

## Running the site locally

```bash
bundle exec jekyll serve
```

The site will be available at `http://localhost:4000`.

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

Visual tests are tagged `@visual` and run against 5 browser/device projects (desktop Chrome, Firefox, Safari, mobile Chrome, mobile Safari).

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

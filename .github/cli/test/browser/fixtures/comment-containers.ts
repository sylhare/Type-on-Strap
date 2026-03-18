interface Attributes {
  [key: string]: string | null | undefined;
}

function buildAttributes(attrs: Attributes): string {
  return Object.entries(attrs)
    .filter(([_, value]) => value !== undefined && value !== null)
    .map(([key, value]) => `${key}="${value}"`)
    .join(' ');
}

interface CusdisOptions {
  lazyLoad?: string | null;
  lang?: string | null;
}

export function cusdisContainer({ lazyLoad = 'true', lang }: CusdisOptions = {}): string {
  const attrs = buildAttributes({ id: 'cusdis_thread', 'data-lazy-load': lazyLoad, 'data-lang': lang });
  return `<div ${attrs}></div>`;
}

interface DisqusOptions {
  lazyLoad?: string | null;
  shortname?: string | null;
}

export function disqusContainer({ lazyLoad = 'true', shortname }: DisqusOptions = {}): string {
  const attrs = buildAttributes({ id: 'disqus_thread', 'data-lazy-load': lazyLoad, 'data-shortname': shortname });
  return `<div ${attrs}></div>`;
}

interface GiscusOptions {
  lazyLoad?: string | null;
  repo?: string | null;
  repoId?: string | null;
  category?: string | null;
  mapping?: string | null;
  theme?: string | null;
}

export function giscusContainer({ lazyLoad = 'true', repo, repoId, category, mapping, theme }: GiscusOptions = {}): string {
  const attrs = buildAttributes({
    id: 'giscus_thread',
    'data-lazy-load': lazyLoad,
    'data-repo': repo,
    'data-repo-id': repoId,
    'data-category': category,
    'data-mapping': mapping,
    'data-theme': theme,
  });
  return `<div ${attrs}></div>`;
}

interface UtterancesOptions {
  lazyLoad?: string | null;
  repo?: string | null;
  issueTerm?: string | null;
  theme?: string | null;
  label?: string | null;
}

export function utterancesContainer({ lazyLoad = 'true', repo, issueTerm, theme, label }: UtterancesOptions = {}): string {
  const attrs = buildAttributes({
    id: 'utterances_thread',
    'data-lazy-load': lazyLoad,
    'data-repo': repo,
    'data-issue-term': issueTerm,
    'data-theme': theme,
    'data-label': label,
  });
  return `<div ${attrs}></div>`;
}

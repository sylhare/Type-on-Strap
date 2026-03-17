/**
 * Helper to build HTML attributes from an object
 */
function buildAttributes(attrs) {
    return Object.entries(attrs)
        .filter(([_, value]) => value !== undefined && value !== null)
        .map(([key, value]) => `${key}="${value}"`)
        .join(' ');
}

function cusdisContainer({ lazyLoad = 'true', lang } = {}) {
    const attrs = buildAttributes({
        id: 'cusdis_thread',
        'data-lazy-load': lazyLoad,
        'data-lang': lang,
    });
    return `<div ${attrs}></div>`;
}

function disqusContainer({ lazyLoad = 'true', shortname } = {}) {
    const attrs = buildAttributes({
        id: 'disqus_thread',
        'data-lazy-load': lazyLoad,
        'data-shortname': shortname,
    });
    return `<div ${attrs}></div>`;
}

function giscusContainer({ lazyLoad = 'true', repo, repoId, category, mapping, theme } = {}) {
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

function utterancesContainer({ lazyLoad = 'true', repo, issueTerm, theme, label } = {}) {
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

module.exports = {
    cusdisContainer,
    disqusContainer,
    giscusContainer,
    utterancesContainer,
};


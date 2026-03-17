/**
 * @fileoverview Unit tests for comments-lazy-load.js functionality
 * @module comments-lazy-load.test
 */

const {
    cusdisContainer,
    disqusContainer,
    giscusContainer,
    utterancesContainer
} = require('./fixtures/comment-containers');

describe('Comments Lazy Loading', () => {
    let mockIntersectionObserver;
    let observerCallback;

    beforeEach(() => {
        document.head.innerHTML = '';
        document.body.innerHTML = '';
        delete window.disqus_shortname;

        observerCallback = null;
        mockIntersectionObserver = jest.fn((callback) => {
            observerCallback = callback;
            this.observe = jest.fn();
            this.unobserve = jest.fn();
            this.disconnect = jest.fn();
            return this;
        });

        global.IntersectionObserver = mockIntersectionObserver;
    });

    afterEach(() => {
        jest.clearAllMocks();
    });

    function loadCommentsScript() {
        const fs = require('fs');
        const path = require('path');
        const scriptPath = path.join(__dirname, '../../../../assets/js/comments-lazy-load.js');
        const scriptContent = fs.readFileSync(scriptPath, 'utf8');
        eval(scriptContent);
    }

    describe('Initialization and Container Detection', () => {
        test.each([
            ['container does not exist', ''],
            ['data-lazy-load is not "true"', cusdisContainer({ lazyLoad: 'false' })],
            ['data-lazy-load attribute is missing', cusdisContainer({ lazyLoad: null })]
        ])('should not initialize when %s', (_, html) => {
            document.body.innerHTML = html;
            loadCommentsScript();
            expect(mockIntersectionObserver).not.toHaveBeenCalled();
        });

        test('should initialize when container exists with data-lazy-load="true"', () => {
            document.body.innerHTML = cusdisContainer();
            loadCommentsScript();
            expect(mockIntersectionObserver).toHaveBeenCalled();
        });

        test('should only initialize the comment system that exists on the page', () => {
            document.body.innerHTML = disqusContainer({ shortname: 'test' });
            loadCommentsScript();
            expect(mockIntersectionObserver).toHaveBeenCalledTimes(1);
        });

        test('should handle multiple comment systems but only initialize those present', () => {
            document.body.innerHTML = cusdisContainer() + giscusContainer({ repo: 'test/repo' });
            loadCommentsScript();
            expect(mockIntersectionObserver).toHaveBeenCalledTimes(2);
        });
    });

    describe('IntersectionObserver Configuration', () => {
        test('should configure observer with 400px rootMargin', () => {
            document.body.innerHTML = cusdisContainer();
            loadCommentsScript();

            expect(mockIntersectionObserver).toHaveBeenCalledWith(
                expect.any(Function),
                expect.objectContaining({
                    rootMargin: '400px'
                })
            );
        });

        test('should observe the container element', () => {
            document.body.innerHTML = cusdisContainer();
            const container = document.getElementById('cusdis_thread');
            loadCommentsScript();

            const observerInstance = mockIntersectionObserver.mock.results[0].value;
            expect(observerInstance.observe).toHaveBeenCalledWith(container);
        });
    });

    describe('Lazy Loading Trigger', () => {
        test.each([
            ['load comments when element intersects', true, 'cusdis.es.js', true],
            ['not load comments when element does not intersect', false, '', false]
        ])('should %s', (_, isIntersecting, scriptMatch, shouldLoad) => {
            document.body.innerHTML = cusdisContainer();
            loadCommentsScript();
            observerCallback([{ isIntersecting }]);

            const scripts = Array.from(document.querySelectorAll('script'));
            if (shouldLoad) {
                expect(scripts.some(s => s.src.includes(scriptMatch))).toBe(true);
            } else {
                expect(scripts.length).toBe(0);
            }
        });

        test('should unobserve after loading once', () => {
            document.body.innerHTML = cusdisContainer();
            const container = document.getElementById('cusdis_thread');
            loadCommentsScript();

            const observerInstance = mockIntersectionObserver.mock.results[0].value;
            observerCallback([{ isIntersecting: true }]);

            expect(observerInstance.unobserve).toHaveBeenCalledWith(container);
        });

        test('should only load once even with multiple intersection events', () => {
            document.body.innerHTML = cusdisContainer();
            loadCommentsScript();

            observerCallback([{ isIntersecting: true }]);
            observerCallback([{ isIntersecting: true }]);
            observerCallback([{ isIntersecting: true }]);

            const scripts = Array.from(document.querySelectorAll('script'));
            const cusdisScripts = scripts.filter(s => s.src.includes('cusdis.es.js'));
            expect(cusdisScripts.length).toBe(1);
        });
    });

    describe('Fallback without IntersectionObserver', () => {
        beforeEach(() => {
            delete window.IntersectionObserver;
        });

        test('should load immediately on complete document', () => {
            Object.defineProperty(document, 'readyState', {
                value: 'complete',
                writable: true,
                configurable: true
            });

            document.body.innerHTML = cusdisContainer();
            loadCommentsScript();

            const scripts = Array.from(document.querySelectorAll('script'));
            expect(scripts.some(s => s.src.includes('cusdis.es.js'))).toBe(true);
        });

        test('should load on window load event if document not complete', () => {
            Object.defineProperty(document, 'readyState', {
                value: 'loading',
                writable: true,
                configurable: true
            });

            document.body.innerHTML = cusdisContainer();
            loadCommentsScript();

            let scripts = Array.from(document.querySelectorAll('script'));
            expect(scripts.length).toBe(0);

            window.dispatchEvent(new Event('load'));

            scripts = Array.from(document.querySelectorAll('script'));
            expect(scripts.some(s => s.src.includes('cusdis.es.js'))).toBe(true);
        });
    });

    describe('Cusdis Loader', () => {
        test('should load Cusdis main script', () => {
            document.body.innerHTML = cusdisContainer();
            loadCommentsScript();
            observerCallback([{ isIntersecting: true }]);

            const scripts = Array.from(document.querySelectorAll('script'));
            const mainScript = scripts.find(s => s.src === 'https://cusdis.com/js/cusdis.es.js');
            expect(mainScript).toBeTruthy();
            expect(mainScript.async).toBe(true);
        });

        test.each([
            ['load Cusdis language script when data-lang is present', 'fr', 'https://cusdis.com/js/widget/lang/fr.js', true],
            ['not load language script when data-lang is absent', null, 'widget/lang/', false]
        ])('should %s', (_, lang, scriptMatch, shouldLoad) => {
            document.body.innerHTML = cusdisContainer({ lang });
            loadCommentsScript();
            observerCallback([{ isIntersecting: true }]);

            const scripts = Array.from(document.querySelectorAll('script'));
            if (shouldLoad) {
                const langScript = scripts.find(s => s.src === scriptMatch);
                expect(langScript).toBeTruthy();
                expect(langScript.async).toBe(true);
            } else {
                const langScripts = scripts.filter(s => s.src.includes(scriptMatch));
                expect(langScripts.length).toBe(0);
            }
        });
    });

    describe('Disqus Loader', () => {
        test('should load Disqus script with correct shortname', () => {
            document.body.innerHTML = disqusContainer({ shortname: 'mysite' });
            loadCommentsScript();
            observerCallback([{ isIntersecting: true }]);

            const scripts = Array.from(document.head.querySelectorAll('script')).concat(
                Array.from(document.body.querySelectorAll('script'))
            );
            const disqusScript = scripts.find(s => s.src.includes('mysite.disqus.com/embed.js'));
            expect(disqusScript).toBeTruthy();
            expect(disqusScript.async).toBe(true);
        });

        test('should set window.disqus_shortname', () => {
            document.body.innerHTML = disqusContainer({ shortname: 'mysite' });
            loadCommentsScript();
            observerCallback([{ isIntersecting: true }]);

            expect(window.disqus_shortname).toBe('mysite');
        });

        test('should log error when shortname is missing', () => {
            console.error = jest.fn();
            document.body.innerHTML = disqusContainer();
            loadCommentsScript();
            observerCallback([{ isIntersecting: true }]);

            expect(console.error).toHaveBeenCalledWith('Disqus shortname not provided');
        });
    });

    describe('Giscus Loader', () => {
        test('should load Giscus script and transfer data attributes', () => {
            document.body.innerHTML = giscusContainer({
                repo: 'user/repo',
                repoId: '123',
                category: 'General',
                mapping: 'pathname',
                theme: 'light'
            });
            loadCommentsScript();
            observerCallback([{ isIntersecting: true }]);

            const giscusScript = document.querySelector('#giscus_thread script[src="https://giscus.app/client.js"]');
            expect(giscusScript).toBeTruthy();
            expect(giscusScript.async).toBe(true);
            expect(giscusScript.crossOrigin).toBe('anonymous');
            expect(giscusScript.getAttribute('data-repo')).toBe('user/repo');
            expect(giscusScript.getAttribute('data-repo-id')).toBe('123');
            expect(giscusScript.getAttribute('data-category')).toBe('General');
            expect(giscusScript.getAttribute('data-theme')).toBe('light');
        });

        test('should not transfer data-lazy-load attribute to script', () => {
            document.body.innerHTML = giscusContainer({ repo: 'user/repo' });
            loadCommentsScript();
            observerCallback([{ isIntersecting: true }]);

            const giscusScript = document.querySelector('#giscus_thread script');
            expect(giscusScript.getAttribute('data-lazy-load')).toBeNull();
        });
    });

    describe('Utterances Loader', () => {
        test('should load Utterances script with attributes', () => {
            document.body.innerHTML = utterancesContainer({
                repo: 'user/repo',
                issueTerm: 'pathname',
                theme: 'github-light',
                label: 'comments'
            });
            loadCommentsScript();
            observerCallback([{ isIntersecting: true }]);

            const utterancesScript = document.querySelector('#utterances_thread script[src="https://utteranc.es/client.js"]');
            expect(utterancesScript).toBeTruthy();
            expect(utterancesScript.async).toBe(true);
            expect(utterancesScript.crossOrigin).toBe('anonymous');
            expect(utterancesScript.getAttribute('repo')).toBe('user/repo');
            expect(utterancesScript.getAttribute('issue-term')).toBe('pathname');
            expect(utterancesScript.getAttribute('theme')).toBe('github-light');
            expect(utterancesScript.getAttribute('label')).toBe('comments');
        });

        test('should handle optional label attribute', () => {
            document.body.innerHTML = utterancesContainer({
                repo: 'user/repo',
                issueTerm: 'pathname',
                theme: 'github-light'
            });
            loadCommentsScript();
            observerCallback([{ isIntersecting: true }]);

            const utterancesScript = document.querySelector('#utterances_thread script');
            expect(utterancesScript.getAttribute('label')).toBeNull();
        });
    });

    describe('Single Comment System Guarantee', () => {
        function getAllScripts() {
            return Array.from(document.head.querySelectorAll('script')).concat(
                Array.from(document.body.querySelectorAll('script'))
            );
        }

        const commentSystems = [
            {
                name: 'Cusdis',
                html: cusdisContainer(),
                scriptUrl: 'https://cusdis.com/js/cusdis.es.js',
                selector: null,
                scriptMatch: 'cusdis',
                otherSystems: ['disqus', 'giscus', 'utteranc']
            },
            {
                name: 'Disqus',
                html: disqusContainer({ shortname: 'test' }),
                scriptUrl: null,
                selector: null,
                scriptMatch: 'disqus',
                otherSystems: ['cusdis', 'giscus', 'utteranc']
            },
            {
                name: 'Giscus',
                html: giscusContainer({ repo: 'test/repo' }),
                scriptUrl: 'https://giscus.app/client.js',
                selector: '#giscus_thread script',
                scriptMatch: 'giscus',
                otherSystems: ['cusdis', 'disqus', 'utteranc']
            },
            {
                name: 'Utterances',
                html: utterancesContainer({ repo: 'test/repo', issueTerm: 'pathname', theme: 'github-light' }),
                scriptUrl: 'https://utteranc.es/client.js',
                selector: '#utterances_thread script',
                scriptMatch: 'utteranc',
                otherSystems: ['cusdis', 'disqus', 'giscus']
            }
        ];

        test.each(commentSystems)(
            'should only load $name when it is the only system present',
            ({ html, scriptUrl, selector, scriptMatch, otherSystems }) => {
                document.head.innerHTML = '';
                document.body.innerHTML = html;

                loadCommentsScript();
                observerCallback([{ isIntersecting: true }]);

                if (selector && scriptUrl) {
                    const systemScript = document.querySelector(selector);
                    expect(systemScript).toBeTruthy();
                    expect(systemScript.src).toBe(scriptUrl);
                }

                const scripts = getAllScripts();
                expect(scripts.some(s => s.src.includes(scriptMatch))).toBe(true);

                otherSystems.forEach(otherSystem => {
                    expect(scripts.some(s => s.src.includes(otherSystem))).toBe(false);
                });
            }
        );
    });
});


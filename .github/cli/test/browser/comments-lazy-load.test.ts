import { cusdisContainer, disqusContainer, giscusContainer, utterancesContainer } from './fixtures/comment-containers';
import fs from 'node:fs';
import * as path from 'node:path';

describe('Comments Lazy Loading', () => {
  let mockIntersectionObserver: jest.Mock;
  let observerCallback: ((entries: Array<{ isIntersecting: boolean }>) => void) | null;
  let mockObserverInstance: { observe: jest.Mock; unobserve: jest.Mock; disconnect: jest.Mock };

  beforeEach(() => {
    document.head.innerHTML = '';
    document.body.innerHTML = '';
    delete (window as any).disqus_shortname;

    observerCallback = null;
    mockObserverInstance = {
      observe: jest.fn(),
      unobserve: jest.fn(),
      disconnect: jest.fn(),
    };

    mockIntersectionObserver = jest.fn((callback: (entries: Array<{ isIntersecting: boolean }>) => void) => {
      observerCallback = callback;
      return mockObserverInstance as unknown as IntersectionObserver;
    });

    (global as any).IntersectionObserver = mockIntersectionObserver;
  });

  afterEach(() => {
    jest.clearAllMocks();
  });

  function loadCommentsScript() {
    const scriptPath = path.join(__dirname, '../../../../assets/js/comments-lazy-load.js');
    const scriptContent = fs.readFileSync(scriptPath, 'utf8');
    eval(scriptContent);
  }

  describe('Initialization and Container Detection', () => {
    test.each([
      ['container does not exist', ''],
      ['data-lazy-load is not "true"', cusdisContainer({ lazyLoad: 'false' })],
      ['data-lazy-load attribute is missing', cusdisContainer({ lazyLoad: null })],
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
        expect.objectContaining({ rootMargin: '400px' }),
      );
    });

    test('should observe the container element', () => {
      document.body.innerHTML = cusdisContainer();
      const container = document.getElementById('cusdis_thread');
      loadCommentsScript();

      expect(mockObserverInstance.observe).toHaveBeenCalledWith(container);
    });
  });

  describe('Lazy Loading Trigger', () => {
    test.each([
      ['load comments when element intersects', true, 'cusdis.es.js', true],
      ['not load comments when element does not intersect', false, '', false],
    ])('should %s', (_, isIntersecting, scriptMatch, shouldLoad) => {
      document.body.innerHTML = cusdisContainer();
      loadCommentsScript();
      observerCallback!([{ isIntersecting: isIntersecting as boolean }]);

      const scripts = Array.from(document.querySelectorAll('script'));
      if (shouldLoad) {
        expect(scripts.some(s => (s as HTMLScriptElement).src.includes(scriptMatch as string))).toEqual(true);
      } else {
        expect(scripts.length).toEqual(0);
      }
    });

    test('should unobserve after loading once', () => {
      document.body.innerHTML = cusdisContainer();
      const container = document.getElementById('cusdis_thread');
      loadCommentsScript();

      observerCallback!([{ isIntersecting: true }]);

      expect(mockObserverInstance.unobserve).toHaveBeenCalledWith(container);
    });

    test('should only load once even with multiple intersection events', () => {
      document.body.innerHTML = cusdisContainer();
      loadCommentsScript();

      observerCallback!([{ isIntersecting: true }]);
      observerCallback!([{ isIntersecting: true }]);
      observerCallback!([{ isIntersecting: true }]);

      const scripts = Array.from(document.querySelectorAll('script'));
      const cusdisScripts = scripts.filter(s => (s as HTMLScriptElement).src.includes('cusdis.es.js'));
      expect(cusdisScripts.length).toEqual(1);
    });
  });

  describe('Fallback without IntersectionObserver', () => {
    beforeEach(() => {
      delete (window as any).IntersectionObserver;
    });

    test('should load immediately on complete document', () => {
      Object.defineProperty(document, 'readyState', { value: 'complete', writable: true, configurable: true });

      document.body.innerHTML = cusdisContainer();
      loadCommentsScript();

      const scripts = Array.from(document.querySelectorAll('script'));
      expect(scripts.some(s => (s as HTMLScriptElement).src.includes('cusdis.es.js'))).toEqual(true);
    });

    test('should load on window load event if document not complete', () => {
      Object.defineProperty(document, 'readyState', { value: 'loading', writable: true, configurable: true });

      document.body.innerHTML = cusdisContainer();
      loadCommentsScript();

      let scripts = Array.from(document.querySelectorAll('script'));
      expect(scripts.length).toEqual(0);

      window.dispatchEvent(new Event('load'));

      scripts = Array.from(document.querySelectorAll('script'));
      expect(scripts.some(s => (s as HTMLScriptElement).src.includes('cusdis.es.js'))).toEqual(true);
    });
  });

  describe('Cusdis Loader', () => {
    test('should load Cusdis main script', () => {
      document.body.innerHTML = cusdisContainer();
      loadCommentsScript();
      observerCallback!([{ isIntersecting: true }]);

      const scripts = Array.from(document.querySelectorAll('script'));
      const mainScript = scripts.find(s => (s as HTMLScriptElement).src === 'https://cusdis.com/js/cusdis.es.js') as HTMLScriptElement | undefined;
      expect(mainScript).toBeTruthy();
      expect(mainScript!.async).toEqual(true);
    });

    test.each([
      ['load Cusdis language script when data-lang is present', 'fr', 'https://cusdis.com/js/widget/lang/fr.js', true],
      ['not load language script when data-lang is absent', null, 'widget/lang/', false],
    ])('should %s', (_, lang, scriptMatch, shouldLoad) => {
      document.body.innerHTML = cusdisContainer({ lang: lang as string | null });
      loadCommentsScript();
      observerCallback!([{ isIntersecting: true }]);

      const scripts = Array.from(document.querySelectorAll('script'));
      if (shouldLoad) {
        const langScript = scripts.find(s => (s as HTMLScriptElement).src === scriptMatch) as HTMLScriptElement | undefined;
        expect(langScript).toBeTruthy();
        expect(langScript!.async).toEqual(true);
      } else {
        const langScripts = scripts.filter(s => (s as HTMLScriptElement).src.includes(scriptMatch as string));
        expect(langScripts.length).toEqual(0);
      }
    });
  });

  describe('Disqus Loader', () => {
    test('should load Disqus script with correct shortname', () => {
      document.body.innerHTML = disqusContainer({ shortname: 'mysite' });
      loadCommentsScript();
      observerCallback!([{ isIntersecting: true }]);

      const scripts = [
        ...Array.from(document.head.querySelectorAll('script')),
        ...Array.from(document.body.querySelectorAll('script')),
      ];
      const disqusScript = scripts.find(s => (s as HTMLScriptElement).src.includes('mysite.disqus.com/embed.js')) as HTMLScriptElement | undefined;
      expect(disqusScript).toBeTruthy();
      expect(disqusScript!.async).toEqual(true);
    });

    test('should set window.disqus_shortname', () => {
      document.body.innerHTML = disqusContainer({ shortname: 'mysite' });
      loadCommentsScript();
      observerCallback!([{ isIntersecting: true }]);

      expect((window as any).disqus_shortname).toEqual('mysite');
    });

    test('should log error when shortname is missing', () => {
      console.error = jest.fn();
      document.body.innerHTML = disqusContainer();
      loadCommentsScript();
      observerCallback!([{ isIntersecting: true }]);

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
        theme: 'light',
      });
      loadCommentsScript();
      observerCallback!([{ isIntersecting: true }]);

      const giscusScript = document.querySelector('#giscus_thread script[src="https://giscus.app/client.js"]') as HTMLScriptElement | null;
      expect(giscusScript).toBeTruthy();
      expect(giscusScript!.async).toEqual(true);
      expect(giscusScript!.crossOrigin).toEqual('anonymous');
      expect(giscusScript!.getAttribute('data-repo')).toEqual('user/repo');
      expect(giscusScript!.getAttribute('data-repo-id')).toEqual('123');
      expect(giscusScript!.getAttribute('data-category')).toEqual('General');
      expect(giscusScript!.getAttribute('data-theme')).toEqual('light');
    });

    test('should not transfer data-lazy-load attribute to script', () => {
      document.body.innerHTML = giscusContainer({ repo: 'user/repo' });
      loadCommentsScript();
      observerCallback!([{ isIntersecting: true }]);

      const giscusScript = document.querySelector('#giscus_thread script') as HTMLScriptElement | null;
      expect(giscusScript!.getAttribute('data-lazy-load')).toBeNull();
    });
  });

  describe('Utterances Loader', () => {
    test('should load Utterances script with attributes', () => {
      document.body.innerHTML = utterancesContainer({
        repo: 'user/repo',
        issueTerm: 'pathname',
        theme: 'github-light',
        label: 'comments',
      });
      loadCommentsScript();
      observerCallback!([{ isIntersecting: true }]);

      const utterancesScript = document.querySelector('#utterances_thread script[src="https://utteranc.es/client.js"]') as HTMLScriptElement | null;
      expect(utterancesScript).toBeTruthy();
      expect(utterancesScript!.async).toEqual(true);
      expect(utterancesScript!.crossOrigin).toEqual('anonymous');
      expect(utterancesScript!.getAttribute('repo')).toEqual('user/repo');
      expect(utterancesScript!.getAttribute('issue-term')).toEqual('pathname');
      expect(utterancesScript!.getAttribute('theme')).toEqual('github-light');
      expect(utterancesScript!.getAttribute('label')).toEqual('comments');
    });

    test('should handle optional label attribute', () => {
      document.body.innerHTML = utterancesContainer({
        repo: 'user/repo',
        issueTerm: 'pathname',
        theme: 'github-light'
      });
      loadCommentsScript();
      observerCallback!([{ isIntersecting: true }]);

      const utterancesScript = document.querySelector('#utterances_thread script') as HTMLScriptElement | null;
      expect(utterancesScript!.getAttribute('label')).toBeNull();
    });
  });

  describe('Single Comment System Guarantee', () => {
    function getAllScripts(): HTMLScriptElement[] {
      return [
        ...Array.from(document.head.querySelectorAll('script')),
        ...Array.from(document.body.querySelectorAll('script')),
      ] as HTMLScriptElement[];
    }

    const commentSystems = [
      {
        name: 'Cusdis',
        html: cusdisContainer(),
        scriptUrl: 'https://cusdis.com/js/cusdis.es.js',
        selector: null,
        scriptMatch: 'cusdis',
        otherSystems: ['disqus', 'giscus', 'utteranc'],
      },
      {
        name: 'Disqus',
        html: disqusContainer({ shortname: 'test' }),
        scriptUrl: null,
        selector: null,
        scriptMatch: 'disqus',
        otherSystems: ['cusdis', 'giscus', 'utteranc'],
      },
      {
        name: 'Giscus',
        html: giscusContainer({ repo: 'test/repo' }),
        scriptUrl: 'https://giscus.app/client.js',
        selector: '#giscus_thread script',
        scriptMatch: 'giscus',
        otherSystems: ['cusdis', 'disqus', 'utteranc'],
      },
      {
        name: 'Utterances',
        html: utterancesContainer({ repo: 'test/repo', issueTerm: 'pathname', theme: 'github-light' }),
        scriptUrl: 'https://utteranc.es/client.js',
        selector: '#utterances_thread script',
        scriptMatch: 'utteranc',
        otherSystems: ['cusdis', 'disqus', 'giscus'],
      },
    ];

    test.each(commentSystems)(
      'should only load $name when it is the only system present',
      ({ html, scriptUrl, selector, scriptMatch, otherSystems }) => {
        document.head.innerHTML = '';
        document.body.innerHTML = html;

        loadCommentsScript();
        observerCallback!([{ isIntersecting: true }]);

        if (selector && scriptUrl) {
          const systemScript = document.querySelector(selector) as HTMLScriptElement | null;
          expect(systemScript).toBeTruthy();
          expect(systemScript!.src).toEqual(scriptUrl);
        }

        const scripts = getAllScripts();
        expect(scripts.some(s => s.src.includes(scriptMatch))).toEqual(true);

        otherSystems.forEach(otherSystem => {
          expect(scripts.some(s => s.src.includes(otherSystem))).toEqual(false);
        });
      },
    );
  });
});

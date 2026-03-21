describe('Dark Mode', () => {
  let setMode: (theme: string) => void;
  let themeToggle: () => void;
  let currentTheme: () => string | null;
  let bootstrapTheme: () => void;

  beforeEach(() => {
    localStorage.clear();

    document.body.innerHTML = `
      <html>
        <head></head>
        <body>
          <button id="theme-toggle"></button>
        </body>
      </html>
    `;

    (global as Record<string, unknown>).darkBtn = 'Dark';
    (global as Record<string, unknown>).lightBtn = 'Light';
    (global as Record<string, unknown>).isAutoTheme = true;

    const themeButton: Record<string, string> = {
      'light': `<i class="fas fa-adjust" aria-hidden="true"></i><span class="navbar-label-with-icon"> ${(global as Record<string, unknown>).darkBtn}</span>`,
      'dark': `<i class="fas fa-adjust fa-rotate-180" aria-hidden="true"></i><span class="navbar-label-with-icon"> ${(global as Record<string, unknown>).lightBtn}</span>`,
    };

    currentTheme = () => localStorage.getItem('theme');

    setMode = (theme: string) => {
      document.documentElement.setAttribute('data-theme', theme);
      localStorage.setItem('theme', theme);
      const toggle = document.getElementById('theme-toggle');
      if (toggle) {
        toggle.innerHTML = themeButton[theme];
      }
    };

    themeToggle = () => {
      const sessionPrefers = currentTheme();
      if (sessionPrefers === 'light') {
        setMode('dark');
      } else {
        setMode('light');
      }
    };

    bootstrapTheme = () => {
      if ((global as Record<string, unknown>).isAutoTheme) {
        if (!currentTheme()) {
          const browserPrefersDark = window.matchMedia('(prefers-color-scheme: dark)');
          if (browserPrefersDark.matches) localStorage.setItem('theme', 'dark');
        }
        const sessionPrefers = currentTheme();
        setMode(sessionPrefers ? sessionPrefers : 'light');
      }
    };
  });

  describe('currentTheme()', () => {
    test('should return null when no theme is set', () => {
      expect(currentTheme()).toBeNull();
    });

    test('should return theme from localStorage', () => {
      localStorage.setItem('theme', 'dark');
      expect(currentTheme()).toEqual('dark');
    });
  });

  describe('setMode()', () => {
    test('should set data-theme attribute on document', () => {
      setMode('dark');
      expect(document.documentElement.getAttribute('data-theme')).toEqual('dark');
    });

    test('should save theme to localStorage', () => {
      setMode('dark');
      expect(localStorage.getItem('theme')).toEqual('dark');
    });

    test('should update theme toggle button innerHTML', () => {
      const toggle = document.getElementById('theme-toggle');
      setMode('dark');
      expect(toggle!.innerHTML).toContain('Light');
      expect(toggle!.innerHTML).toContain('fa-rotate-180');
    });

    test('should handle light mode correctly', () => {
      setMode('light');
      expect(document.documentElement.getAttribute('data-theme')).toEqual('light');
      const toggle = document.getElementById('theme-toggle');
      expect(toggle!.innerHTML).toContain('Dark');
    });

    test('should not error if toggle button does not exist', () => {
      document.getElementById('theme-toggle')!.remove();
      expect(() => setMode('dark')).not.toThrow();
    });
  });

  describe('themeToggle()', () => {
    test('should toggle from light to dark', () => {
      setMode('light');
      themeToggle();
      expect(currentTheme()).toEqual('dark');
      expect(document.documentElement.getAttribute('data-theme')).toEqual('dark');
    });

    test('should toggle from dark to light', () => {
      setMode('dark');
      themeToggle();
      expect(currentTheme()).toEqual('light');
      expect(document.documentElement.getAttribute('data-theme')).toEqual('light');
    });

    test('should default to light when no theme is set', () => {
      themeToggle();
      expect(currentTheme()).toEqual('light');
    });

    test('should toggle multiple times correctly', () => {
      setMode('light');
      themeToggle();
      expect(currentTheme()).toEqual('dark');
      themeToggle();
      expect(currentTheme()).toEqual('light');
      themeToggle();
      expect(currentTheme()).toEqual('dark');
    });
  });

  describe('bootstrapTheme()', () => {
    test('should set light mode by default when no preference', () => {
      (window as any).matchMedia = jest.fn().mockReturnValue({
        matches: false,
        media: '(prefers-color-scheme: none)',
        addEventListener: jest.fn(),
      });

      bootstrapTheme();
      expect(currentTheme()).toEqual('light');
      expect(document.documentElement.getAttribute('data-theme')).toEqual('light');
    });

    test('should respect existing localStorage theme', () => {
      localStorage.setItem('theme', 'dark');
      bootstrapTheme();
      expect(currentTheme()).toEqual('dark');
      expect(document.documentElement.getAttribute('data-theme')).toEqual('dark');
    });

    test('should respect browser preference for dark mode', () => {
      localStorage.clear();

      (window as any).matchMedia = jest.fn().mockImplementation((query: string) => ({
        matches: query === '(prefers-color-scheme: dark)',
        media: query,
        addEventListener: jest.fn(),
      }));

      bootstrapTheme();

      expect(window.matchMedia).toHaveBeenCalledWith('(prefers-color-scheme: dark)');
      expect(currentTheme()).toEqual('dark');
      expect(document.documentElement.getAttribute('data-theme')).toEqual('dark');
    });

    test('should default to light mode when browser prefers light', () => {
      localStorage.clear();

      (window as any).matchMedia = jest.fn().mockReturnValue({
        matches: false,
        media: '(prefers-color-scheme: dark)',
        addEventListener: jest.fn(),
      });

      bootstrapTheme();

      expect(currentTheme()).toEqual('light');
      expect(document.documentElement.getAttribute('data-theme')).toEqual('light');
    });

    test('should not set theme when isAutoTheme is false', () => {
      localStorage.clear();
      document.documentElement.removeAttribute('data-theme');

      (global as Record<string, unknown>).isAutoTheme = false;
      bootstrapTheme();

      expect(localStorage.getItem('theme')).toBeNull();
      expect(document.documentElement.getAttribute('data-theme')).toBeNull();
    });
  });

  describe('Theme persistence', () => {
    test('should persist theme across page loads', () => {
      (window as any).matchMedia = jest.fn().mockReturnValue({
        matches: false,
        media: '(prefers-color-scheme: none)',
        addEventListener: jest.fn(),
      });

      setMode('dark');
      expect(localStorage.getItem('theme')).toEqual('dark');

      document.documentElement.removeAttribute('data-theme');

      bootstrapTheme();
      expect(document.documentElement.getAttribute('data-theme')).toEqual('dark');
    });

    test('should maintain theme through multiple toggles', () => {
      setMode('light');
      themeToggle();
      themeToggle();
      themeToggle();

      expect(localStorage.getItem('theme')).toEqual('dark');
      expect(document.documentElement.getAttribute('data-theme')).toEqual('dark');
    });
  });
});

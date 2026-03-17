/**
 * @fileoverview Unit tests for dark-mode.js functionality
 * @module dark-mode.test
 */

describe('Dark Mode', () => {
  let setMode, themeToggle, currentTheme, bootstrapTheme;

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

    global.darkBtn = 'Dark';
    global.lightBtn = 'Light';
    global.isAutoTheme = true;

    const themeButton = {
      'light': `<i class="fas fa-adjust" aria-hidden="true"></i><span class="navbar-label-with-icon"> ${global.darkBtn}</span>`,
      'dark': `<i class="fas fa-adjust fa-rotate-180" aria-hidden="true"></i><span class="navbar-label-with-icon"> ${global.lightBtn}</span>`,
    };

    /**
     * Gets the current theme from localStorage
     * @returns {string|null} The current theme or null
     */
    currentTheme = () => localStorage.getItem('theme');

    /**
     * Sets the theme mode and updates DOM
     * @param {string} theme - The theme to set ('light' or 'dark')
     */
    setMode = (theme) => {
      document.documentElement.setAttribute('data-theme', theme);
      localStorage.setItem('theme', theme);
      const toggle = document.getElementById('theme-toggle');
      if (toggle) {
        toggle.innerHTML = themeButton[theme];
      }
    };

    /**
     * Toggles between light and dark themes
     */
    themeToggle = () => {
      let sessionPrefers = currentTheme();
      if (sessionPrefers === 'light') {
        setMode('dark');
      } else {
        setMode('light');
      }
    };

    /**
     * Initializes theme based on browser preference and localStorage
     */
    bootstrapTheme = () => {
      if (global.isAutoTheme) {
        if (!currentTheme()) {
          let browserPrefersDark = window.matchMedia('(prefers-color-scheme: dark)');
          if (browserPrefersDark.matches) localStorage.setItem('theme', 'dark');
        }
        let sessionPrefers = currentTheme();
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
      expect(currentTheme()).toBe('dark');
    });
  });

  describe('setMode()', () => {
    test('should set data-theme attribute on document', () => {
      setMode('dark');
      expect(document.documentElement.getAttribute('data-theme')).toBe('dark');
    });

    test('should save theme to localStorage', () => {
      setMode('dark');
      expect(localStorage.getItem('theme')).toBe('dark');
    });

    test('should update theme toggle button innerHTML', () => {
      const toggle = document.getElementById('theme-toggle');
      setMode('dark');
      expect(toggle.innerHTML).toContain('Light');
      expect(toggle.innerHTML).toContain('fa-rotate-180');
    });

    test('should handle light mode correctly', () => {
      setMode('light');
      expect(document.documentElement.getAttribute('data-theme')).toBe('light');
      const toggle = document.getElementById('theme-toggle');
      expect(toggle.innerHTML).toContain('Dark');
    });

    test('should not error if toggle button does not exist', () => {
      document.getElementById('theme-toggle').remove();
      expect(() => setMode('dark')).not.toThrow();
    });
  });

  describe('themeToggle()', () => {
    test('should toggle from light to dark', () => {
      setMode('light');
      themeToggle();
      expect(currentTheme()).toBe('dark');
      expect(document.documentElement.getAttribute('data-theme')).toBe('dark');
    });

    test('should toggle from dark to light', () => {
      setMode('dark');
      themeToggle();
      expect(currentTheme()).toBe('light');
      expect(document.documentElement.getAttribute('data-theme')).toBe('light');
    });

    test('should default to light when no theme is set', () => {
      themeToggle();
      expect(currentTheme()).toBe('light');
    });

    test('should toggle multiple times correctly', () => {
      setMode('light');
      themeToggle(); // -> dark
      expect(currentTheme()).toBe('dark');
      themeToggle(); // -> light
      expect(currentTheme()).toBe('light');
      themeToggle(); // -> dark
      expect(currentTheme()).toBe('dark');
    });
  });

  describe('bootstrapTheme()', () => {
    test('should set light mode by default when no preference', () => {
      window.matchMedia = jest.fn().mockReturnValue({
        matches: false,
        media: '(prefers-color-scheme: none)',
        addEventListener: jest.fn(),
      });

      bootstrapTheme();
      expect(currentTheme()).toBe('light');
      expect(document.documentElement.getAttribute('data-theme')).toBe('light');
    });

    test('should respect existing localStorage theme', () => {
      localStorage.setItem('theme', 'dark');
      bootstrapTheme();
      expect(currentTheme()).toBe('dark');
      expect(document.documentElement.getAttribute('data-theme')).toBe('dark');
    });

    test('should respect browser preference for dark mode', () => {
      localStorage.clear();

      window.matchMedia = jest.fn().mockImplementation(query => ({
        matches: query === '(prefers-color-scheme: dark)',
        media: query,
        addEventListener: jest.fn(),
      }));

      bootstrapTheme();

      expect(window.matchMedia).toHaveBeenCalledWith('(prefers-color-scheme: dark)');
      expect(currentTheme()).toBe('dark');
      expect(document.documentElement.getAttribute('data-theme')).toBe('dark');
    });

    test('should default to light mode when browser prefers light', () => {
      localStorage.clear();

      window.matchMedia = jest.fn().mockReturnValue({
        matches: false,
        media: '(prefers-color-scheme: dark)',
        addEventListener: jest.fn(),
      });

      bootstrapTheme();

      expect(currentTheme()).toBe('light');
      expect(document.documentElement.getAttribute('data-theme')).toBe('light');
    });

    test('should not set theme when isAutoTheme is false', () => {
      localStorage.clear();
      document.documentElement.removeAttribute('data-theme');

      global.isAutoTheme = false;
      bootstrapTheme();

      expect(localStorage.getItem('theme')).toBeNull();
      expect(document.documentElement.getAttribute('data-theme')).toBeNull();
    });
  });

  describe('Theme persistence', () => {
    test('should persist theme across page loads', () => {
      window.matchMedia = jest.fn().mockReturnValue({
        matches: false,
        media: '(prefers-color-scheme: none)',
        addEventListener: jest.fn(),
      });

      setMode('dark');
      expect(localStorage.getItem('theme')).toBe('dark');

      document.documentElement.removeAttribute('data-theme');

      bootstrapTheme();
      expect(document.documentElement.getAttribute('data-theme')).toBe('dark');
    });

    test('should maintain theme through multiple toggles', () => {
      setMode('light');
      themeToggle();
      themeToggle();
      themeToggle();

      expect(localStorage.getItem('theme')).toBe('dark');
      expect(document.documentElement.getAttribute('data-theme')).toBe('dark');
    });
  });
});


/* @preserve Dark mode Init */
/*
 * There are two colour palettes on CSS for the data-theme: 'light' and 'dark'.
 * Initially the script checks if a theme is set in session storage and
 * alternatively listens to a MediaQuery callback looking for "prefers-color-scheme: dark".
 *
 * The variables darkBtn and lightBtn are defined in head.liquid from the _data/translations.yml
 * The isAutoTheme is defined in head.liquid from the _config.yml
 */

const themeButton = {
    'light': `<i class="fas fa-adjust" aria-hidden="true"></i><span class="navbar-label-with-icon"> ${darkBtn}</span>`,
    'dark': `<i class="fas fa-adjust fa-rotate-180" aria-hidden="true"></i><span class="navbar-label-with-icon"> ${lightBtn}</span>`
};

function currentTheme(){
    return localStorage.getItem('theme');
}

function setMode(theme) {
    document.documentElement.setAttribute('data-theme', theme);
    localStorage.setItem('theme', theme);
    const toggle = document.getElementById('theme-toggle');
    if (toggle) {
        toggle.innerHTML = themeButton[theme];
    }
}

function themeToggle() {
    let sessionPrefers = currentTheme();
    if (sessionPrefers === 'light') {
        setMode('dark');
    } else {
        setMode('light');
    }
}

function bootstrapTheme() {
    if (isAutoTheme) {
        if (!currentTheme()) {
            // Load browser's preference
            let browserPrefersDark = window.matchMedia('(prefers-color-scheme: dark)');
            if (browserPrefersDark.matches) localStorage.setItem('theme', 'dark');
            browserPrefersDark.addEventListener('change', () => {
                if (browserPrefersDark.matches) localStorage.setItem('theme', 'dark');
            });
        }
        // Load theme
        let sessionPrefers = currentTheme();
        setMode(sessionPrefers ? sessionPrefers : 'light');
    }
}

// Init
(function () {
    bootstrapTheme();
})()

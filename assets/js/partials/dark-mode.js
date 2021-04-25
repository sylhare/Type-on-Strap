// There are two color palettes on CSS, one for 'data-theme === light' and one for 'data-theme === dark'.
// Initially the script check if a theme is set in session storage and alternatively listens to a MediaQuery callback looking for "prefers-color-scheme: dark".

// HTML of the dark theme button
const themeButton = {
    'light': '<i class="fas fa-adjust" aria-hidden="true"></i>',
    'dark': '<i class="fas fa-adjust fa-rotate-180" aria-hidden="true"></i>'
}

// change theme
function setMode(theme) {
    document.documentElement.setAttribute('data-theme', theme)
    sessionStorage.setItem('theme', theme)
    document.getElementById('theme-toggle').innerHTML = themeButton[theme]
}

// theme switcher button callback
function themeToggle() {
    let sessionPrefers = sessionStorage.getItem('theme')
    if (sessionPrefers === 'light') {
        setMode('dark')
    } else {
        setMode('light')
    }
}

window.onload = function bootstrapTheme() {
    if (isAutoTheme) {
        // Load browser's preference
        let browserPrefersDark = window.matchMedia('(prefers-color-scheme: dark)');
        browserPrefersDark.addEventListener('change', () => {
            sessionStorage.setItem('theme', browserPrefersDark.matches ? 'dark' : 'light')
        });

        // Load theme
        let sessionPrefers = sessionStorage.getItem('theme')
        setMode(sessionPrefers ? sessionPrefers : 'light')
    }
}

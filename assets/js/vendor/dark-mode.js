// There are two color palettes on CSS, one for 'data-theme === light' and one for 'data-theme === dark'.
// Initially the script check if a theme is set in session storage and alternatively listens to a MediaQuery callback looking for "prefers-color-scheme: dark".
// If the script is not loaded the toggle button will not be show and hopefully the browser will be smart enough to choose the right palette.
// If the script is loaded, users can use the toggle to change theme.

// TODO: Test it in other environments to see if it's working properly.

/* LIB
   --- */

// HTML of the dark theme button
const buttonLightHTML = '<i class="fas fa-adjust" aria-hidden="true"></i>' 
const buttonDarkHTML = '<i class="fas fa-adjust fa-rotate-180" aria-hidden="true"></i>' 

// update the rotation of the button
function setButton(theme) {
    if (theme === 'light') {
        document.getElementById('theme-toggle').innerHTML = buttonLightHTML
    } else if (theme === 'dark') {
        document.getElementById('theme-toggle').innerHTML = buttonDarkHTML 
    }
}

// change theme
function setMode(theme) {
    // change the theme 
	document.documentElement.setAttribute('data-theme', theme)
    // store the change
    sessionStorage.setItem('theme', theme)
    // update button
    setButton(theme)
}

// browser preference callback
function browserPrefers(browserPrefersDark) {
    if (browserPrefersDark.matches) {
        sessionPrefers = 'dark'
    } else {
        sessionPrefers = 'light'
    }
}

// theme switcher button callback
function themeToggle() {
    // recheck the session preference
	let sessionPrefers = sessionStorage.getItem('theme')
	if (sessionPrefers === 'dark') {
		setMode('light')
	} else if (sessionPrefers === 'light') {
		setMode('dark')
    }
}

/* MAIN
   ---- */
window.onload = function bootstrapTheme() {
    // check for session stored preference (string)
    let sessionPrefers = sessionStorage.getItem('theme')

    // check if browser prefers dark theme (boolean)
    // FIXME: method is deprecated but still working
    // https://developer.mozilla.org/en-US/docs/Web/API/MediaQueryList/addListener
    let browserPrefersDark = window.matchMedia('(prefers-color-scheme: dark)'); 

    // if there's no session preference check if browser prefers dark
    if (!sessionPrefers) {
        browserPrefersDark.addListener(browserPrefers)
    }

    // if still no preference is found, set it to light (default)
    if (!sessionPrefers) {
        sessionPrefers = 'light'
    }

    // check what the session prefers
    if (sessionPrefers === 'dark') {
        setMode('dark')
    } else if (sessionPrefers === 'light') {
        setMode('light')
    }
}

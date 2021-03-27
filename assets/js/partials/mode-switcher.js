//@preserve Mode Switcher by derekkedziora
//@preserve License Creative Commons Zero

let systemInitiatedDark = window.matchMedia("(prefers-color-scheme: dark)"); 
let theme = sessionStorage.getItem('theme');

if (systemInitiatedDark.matches) {
	document.getElementById("theme-toggle").innerHTML = '<i class="fas fa-adjust fa-rotate-180" aria-hidden="true"></i>';
} else {
	document.getElementById("theme-toggle").innerHTML = '<i class="fas fa-adjust" aria-hidden="true"></i>';
}

function prefersColorTest(systemInitiatedDark) {
  if (systemInitiatedDark.matches) {
  	document.documentElement.setAttribute('data-theme', 'dark');		
   	document.getElementById("theme-toggle").innerHTML = '<i class="fas fa-adjust fa-rotate-180" aria-hidden="true"></i>';
   	sessionStorage.setItem('theme', '');
  } else {
  	document.documentElement.setAttribute('data-theme', 'light');
    document.getElementById("theme-toggle").innerHTML = '<i class="fas fa-adjust" aria-hidden="true"></i>';
    sessionStorage.setItem('theme', '');
  }
}
systemInitiatedDark.addListener(prefersColorTest);


function modeSwitcher() {
	let theme = sessionStorage.getItem('theme');
	if (theme === "dark") {
		document.documentElement.setAttribute('data-theme', 'light');
		sessionStorage.setItem('theme', 'light');
		document.getElementById("theme-toggle").innerHTML = '<i class="fas fa-adjust" aria-hidden="true"></i>';
	}	else if (theme === "light") {
		document.documentElement.setAttribute('data-theme', 'dark');
		sessionStorage.setItem('theme', 'dark');
		document.getElementById("theme-toggle").innerHTML = '<i class="fas fa-adjust fa-rotate-180" aria-hidden="true"></i>';
	} else if (systemInitiatedDark.matches) {	
		document.documentElement.setAttribute('data-theme', 'light');
		sessionStorage.setItem('theme', 'light');
		document.getElementById("theme-toggle").innerHTML = '<i class="fas fa-adjust" aria-hidden="true"></i>';
	} else {
		document.documentElement.setAttribute('data-theme', 'dark');
		sessionStorage.setItem('theme', 'dark');
		document.getElementById("theme-toggle").innerHTML = '<i class="fas fa-adjust fa-rotate-180" aria-hidden="true"></i>';
	}
}

if (theme === "dark") {
	document.documentElement.setAttribute('data-theme', 'dark');
	sessionStorage.setItem('theme', 'dark');
	document.getElementById("theme-toggle").innerHTML = '<i class="fas fa-adjust fa-rotate-180" aria-hidden="true"></i>';
} else if (theme === "light") {
	document.documentElement.setAttribute('data-theme', 'light');
	sessionStorage.setItem('theme', 'light');
	document.getElementById("theme-toggle").innerHTML = '<i class="fas fa-adjust" aria-hidden="true"></i>';
}

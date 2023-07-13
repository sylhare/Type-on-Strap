/* @preserve Cookie Consent Init */
function createCookie(name, value, days) {
  var expires = "";
  if (days) {
    var date = new Date();
    date.setTime(date.getTime() + (days * 24 * 60 * 60 * 1000));
    expires = "; expires=" + date.toUTCString();
  }
  document.cookie = `${name}=${value}${expires}; path=/`;
}

function readCookie(name) {
  var nameEQ = name + "=";
  var ca = document.cookie.split(';');
  for (var i = 0; i < ca.length; i++) {
    var c = ca[i];
    while (c.charAt(0) === ' ') c = c.substring(1, c.length);
    if (c.indexOf(nameEQ) === 0) return c.substring(nameEQ.length, c.length);
  }
  return null;
}

function addCookieConsentListener() {
  document.getElementById('cookie-notice-accept').addEventListener("click", function () {
    createCookie(cookieName, 'true', 31);
    document.getElementById('cookie-notice').style.display = 'none';
    location.reload();
  });
}

function googleAnalytics() {
  if (analyticsName.toLowerCase() !== '') {
    // Google tag manager
    window.dataLayer = window.dataLayer || [];
    function gtag() { dataLayer.push(arguments); }
    gtag('js', new Date());
    gtag('config', analyticsName, { 'anonymize_ip': true });
    if (analyticsNameGA4) {
      gtag('config', analyticsNameGA4, { 'anonymize_ip': true });
    }
  }
}

if (isCookieConsent.toLowerCase() === 'true') {
  addCookieConsentListener();
  if (readCookie(cookieName) === 'true') {
    googleAnalytics();
  } else {
    document.getElementById('cookie-notice').style.display = 'block';
  }
} else {
  googleAnalytics();
}



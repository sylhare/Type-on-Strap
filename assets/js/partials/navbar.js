document.addEventListener("DOMContentLoaded", function (event) {

  /*
   * Display the menu items on smaller screens
   */
  var pull = document.getElementById('pull');
  var menu = document.querySelector('nav ul');

  pull.addEventListener('click', function (e) {
    menu.classList.toggle('hide');
  });

  /*
   * Make the header images move on scroll
   */
  window.addEventListener('scroll', function () {
    var x = window.pageYOffset || document.body.scrollTop;
    var main = document.getElementById("main");
    var mainStyle = main.style;

    mainStyle.backgroundPosition = '100% ' + parseInt(-x / 3) + 'px' + ', 0%, center top';
  });
});
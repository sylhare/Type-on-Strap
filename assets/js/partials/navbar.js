document.addEventListener("DOMContentLoaded", function(event) { 

    /*
    * Display the menu items on smaller screens
    */
    var pull = document.getElementById('pull');
    var menu = document.querySelector('nav ul');

    pull.addEventListener('click', function(e) {
        menu.classList.toggle('hide');
    });

    /*
    * Display the navbar back to normal after resize
    */
    window.addEventListener('resize', function() {
        var w = window.innerWidth;
        if (w > 320) {
            menu.removeAttribute('style');
        }
    });

    /*
    * Make the header images move on scroll
    */
    window.addEventListener('scroll', function() {
        var x = window.pageYOffset | document.body.scrollTop;
        var m = document.getElementById("main"), c = m.style;
        
        c.backgroundPosition = '100% ' + parseInt(-x/3) + 'px' + ', 0%, center top';
    });
});
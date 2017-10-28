/*
 * Display the menu items on smaller screens
 */
$(function () {
    var pull = $('#pull');
    menu = $('nav ul');
    menuHeight = menu.height();

    $(pull).on('click', function (e) {
        e.preventDefault();
        menu.slideToggle();
    });
});

/*
 * Display the navbar back to normal after resize
 */
$(window).resize(function () {
    var w = $(window).width();
    if (w > 320 && menu.is(':hidden')) {
        menu.removeAttr('style');
    }
});


/*
$('.call-out').scroll(function () {
    var x = $(this).scrollTop();
    $(this).css('background-position', '0% ' + parseInt(-x / 10) + 'px');
});
*/
$(window).scroll(function() {
    var x = $(this).scrollTop();
    $('#main').css('background-position', '100% ' + parseInt(-x / 1) + 'px' + ', 0% ' + parseInt(-x / 2) + 'px, center top');
});
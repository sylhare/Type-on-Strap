/*
 * Display the menu items on smaller screens
 */
$(function() {
  var pull    = $('#pull');
    menu    = $('nav ul');
    menuHeight  = menu.height();
 
  $(pull).on('click', function(e) {
    e.preventDefault();
    menu.slideToggle();
  });
});

/*
 * Display the navbar back to normal after resize
 */
$(window).resize(function(){
  var w = $(window).width();
  if(w > 320 && menu.is(':hidden')) {
    menu.removeAttr('style');
  }
}); 

$(window).scroll(function(){
    var h = 20;    
    var $scrollTop = $(window).scrollTop();
    
    if($scrollTop < 60){
        $(".navbar").css({
          'paddingTop': h - ($scrollTop/3) + "px",
        });
    }
});

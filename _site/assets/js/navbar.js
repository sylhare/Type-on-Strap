$(window).scroll(function(){
    var h = 20;    
    var $scrollTop = $(window).scrollTop();
    
    if($scrollTop < 60){
        $(".scrollbar").css({
          'paddingTop': h - ($scrollTop/3) + "px",
        });
    }
});

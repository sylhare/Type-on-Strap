function toggle() {
    var x = document.getElementById("navbar");
    if (x.className === "navbar") {
        x.className += " responsive";
    } else {
        x.className = "navbar";
    }
}

/* Toggle the active class when clicked on the bar icon */
jQuery(document).ready(function() {
    jQuery('.toggle-nav').click(function(e) {
        jQuery(this).toggleClass('active');
        jQuery('.menu ul').toggleClass('active');
 
        e.preventDefault();
    });
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

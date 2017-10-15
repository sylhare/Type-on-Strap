function toggle() {
    var x = document.getElementById("navbar");
    if (x.className === "navbar") {
        x.className += " responsive";
    } else {
        x.className = "navbar";
    }
}

$(window).scroll(function(){
    var h = 20;    
    var $scrollTop = $(window).scrollTop();
    
    if($scrollTop < 60){
        $(".navbar").css({
          'paddingTop': h - ($scrollTop/3) + "px",
        });
    }
});

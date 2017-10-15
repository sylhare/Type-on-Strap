$(window).scroll(function(){
    var h = 20;    
    var $scrollTop = $(window).scrollTop();
    
    if($scrollTop < 60){
        $(".scrollbar").css({
          'paddingTop': h - ($scrollTop/3) + "px",
        });
    }
});

function myFunction() {
    var x = document.getElementById("myTopnav");
    if (x.className === "topnav") {
        x.className += " responsive";
    } else {
        x.className = "topnav";
    }
}
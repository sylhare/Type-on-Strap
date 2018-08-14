// var wrap = $("header-bar");

// wrap.on("scroll", function(e) {
    
//   if (this.scrollTop > 147) {
//     wrap.addClass("fixed");
//     wrap.removeClass("content");
//   } else {
//     wrap.removeClass("fixed");
//     wrap.addClass("content");
//   }
  
// });

$( document ).ready(function() {
// 	$(window).scroll(function(){
// 	  var fixedHeader = $('.header-bar-fixed'),
// 	  	  fixedSearch = $('.search'),
// 	      scroll = $(window).scrollTop();

// 	  if (scroll >= 100) fixedHeader.removeClass('hide');
// 	  else fixedHeader.addClass('hide');

	  // if (scroll >= 650) fixedSearch.addClass('search-fixed');
	  // else fixedSearch.removeClass('search-fixed');
// 	});


// Sticky Header
	$(window).scroll(function() {

	    if ($(window).scrollTop() > 100) {
	        $('.main_h').addClass('sticky');
	    } else {
	        $('.main_h').removeClass('sticky');
	    }
	});

	// Mobile Navigation
	$('.mobile-toggle').click(function() {
	    if ($('.main_h').hasClass('open-nav')) {
	        $('.main_h').removeClass('open-nav');
	    } else {
	        $('.main_h').addClass('open-nav');
	    }
	});

	$('.main_h li a').click(function() {
	    if ($('.main_h').hasClass('open-nav')) {
	        $('.navigation').removeClass('open-nav');
	        $('.main_h').removeClass('open-nav');
	    }
	});

	// Navigation Scroll
	$('nav a').click(function(event) {
	    var id = $(this).attr("href");
	    var offset = 70;
	    var target = $(id).offset().top - offset;
	    $('html, body').animate({
	        scrollTop: target
	    }, 500);
	    event.preventDefault();
	});

});
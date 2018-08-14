var gData;
/*
    Draw the toolbox with an array of spreadsheet data
    and cache the data in a global variable.
*/
function initialize(data) {
    gData = data
    drawToolBox(data)
}

function drawToolBox(data) {
  var tools = ich.tools({
    'rows': data
  });
    var out = ''
    $.each(data, function(k, v) {
	//console.log(k,v);
	if (v.category2 != '') {
	    var category = v.category +', '+ v.category2
	}
	else {
	    var category = v.category
	}
	out += '<div id="'+ v.id +'" class="tool-box cat-'+ v.categoryid +'">'
        out += '<span class="quote-heading">'+ v.description +'</span>'
    	out += '<div class="tool-box-bottom '+ v.id +'"><ul>'
    	out += '<li class="quote">'+ v.quote +'</li>'
    	out += '<li><strong>'+ v.interviewee +', '+ v.title +'</strong></li>'
    	out += '<li class="category"><small>'+ category +'</small></li>'
    	out += '<li><a href="https://t4sj.co/quotes#'+ v.id +'"><small>LINK</small></a></li>'
    	out += '</ul></div></div>'
    });

    $('#tools').html(out);
    // open a quote if it is in the URL
    if (window.location.hash != '') {
	var rowNumber = window.location.hash.substring(1)
	$('.tool-box-bottom' + '.' + rowNumber).css('display', 'inherit')
	$('div#'+ rowNumber).addClass('selected-tool')
	var position = ($('div#'+ rowNumber).offset().top - 80)
	$(window).scrollTop(position)
    }    
}

$(document).on( 'click', '#showAvailable', toggleAvailable)

$(document).on( 'click', '.clear', function(e) {
  clearSearch(e)
  toggleAvailable(false)
})

$(document).on('keyup search', '#toolSearch', function(e) {
  var text = $(e.target).val().trim().toLowerCase()
  if (text === '') return clearSearch(e)
  filterTools(text)
})

// $(document).on('click submit', '#catSearch', function(e) {
//   var text = $(e.target).val().trim().toLowerCase()
//   console.log(text)
//   if (text === '') return clearSearch(e)
//   filterCategory(text)
// })

$(document).on('keyup', 'form[name=catSearchForm]', function(e) {
    // prevent the form from submitting
    var text = $(e.target).val().trim().toLowerCase()
    console.log(text)
    if (text === '') return clearSearch(e)
    filterCategory(text)
});


$(document).on( 'click', '.quote-heading', function(e) {
  var rowNumber = $(this).closest("div").attr("id")
  if ($(this).closest('div').hasClass('selected-tool')) {
    $('.tool-box-bottom' + '.' + rowNumber).css('display', 'none')
    $(this).closest('div').removeClass('selected-tool')
  }
  else {
    $('.tool-box-bottom' + '.' + rowNumber).css('display', 'inherit')
    $(this).closest('div').addClass('selected-tool')
  }
})

/*
    Toggle whether unavailable items are displayed or not

    [state] - optionally specify the desired state, either true or false. If
    not defined, it will just toggle.
*/
function toggleAvailable(state) {
    // Since #showAvailable and #tools both start out without the class, they
    // will remain in sync
    let button = $('#showAvailable').toggleClass('button-pressed', state)
    // CSS in site.css hides tools that are '.show-available .not-available'
    $('#tools').toggleClass('show-available', state)

    // Update the button text
    if (button.hasClass('button-pressed')) {
        button.html("Show All")
    } else {
        button.html("Show Available")
    }
}

function clearSearch(e) {
  console.log('clear')
  $('#toolSearch').val('')
  $('form :input').val('');
  drawToolBox(gData)
}

function filterTools(text) {
  $('.quote').each(function() {
  var tool = $(this).html().toLowerCase()
  if (tool.match(text)) {
    $(this).parent().parent().parent().removeClass('filtered')
} else $(this).parent().parent().parent().addClass('filtered')
  })
}

function filterCategory(text) {
  $('.category').each(function() {
  var tool = $(this).html().toLowerCase()
  if (tool.match(text)) {
    $(this).parent().parent().parent().removeClass('filtered')
} else $(this).parent().parent().parent().addClass('filtered')
  })
}



/* @preserve Masonry Init */
try {
  var elem = document.querySelector('.grid');
  var msnry = new Masonry(elem, {
    itemSelector: '.grid-item',
    columnWidth: '.grid-sizer',
    gutter: '.gutter-sizer',
    percentPosition: true
  });

  // layout Masonry after each image loads
  var imgLoad = imagesLoaded(elem);
  imgLoad.on('progress', function (instance, image) {
    msnry.layout();
  });
} catch (err) {
  if (err instanceof ReferenceError) {
    // Do nothing, Masonry is defined only in the gallery page
  } else {
    throw err;
  }
}

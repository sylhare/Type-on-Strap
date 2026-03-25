;(function () {
  var grid = document.querySelector('.grid');
  if (!grid) return;

  var columns = [];
  var columnHeights = [];
  var items = [];
  var resizeTimeout;

  function getColumns() {
    return Array.from(grid.querySelectorAll('.grid-col')).filter(function (col) {
      return getComputedStyle(col).display !== 'none';
    });
  }

  function layout() {
    columns = getColumns();
    columnHeights = columns.map(function () { return 0; });
    items = Array.from(grid.querySelectorAll('.grid-item'));
    items.forEach(function (item) {
      var minHeight = Math.min.apply(Math, columnHeights);
      var index = columnHeights.indexOf(minHeight);
      columns[index].appendChild(item);
      columnHeights[index] += item.offsetHeight || 1;
    });
  }

  function measureColumnHeight(elem) {
    var rect = grid.getBoundingClientRect();
    columns.forEach(function (col, i) {
      if (!elem || col.contains(elem)) {
        var last = col.lastElementChild;
        if (last) columnHeights[i] = last.getBoundingClientRect().bottom - rect.top;
      }
    });
  }

  function onResize() {
    clearTimeout(resizeTimeout);
    resizeTimeout = setTimeout(function () {
      var newColumns = getColumns();
      if (newColumns.length !== columns.length ||
          newColumns.some(function (c, i) { return c !== columns[i]; })) {
        columns = newColumns;
        layout();
      }
    }, 100);
  }

  layout();
  window.addEventListener('resize', onResize);
  grid.addEventListener('load', function (e) { measureColumnHeight(e.target); }, true);
})();

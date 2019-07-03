(function(i, s, o, g, r, a, m) {
  i['GoogleAnalyticsObject'] = r;
  i[r] = i[r] || function() {
    (i[r].q = i[r].q || []).push(arguments)
  }, i[r].l = 1 * new Date();
  a = s.createElement(o),
    m = s.getElementsByTagName(o)[0];
  a.async = 1;
  a.src = g;
  m.parentNode.insertBefore(a, m)
})(window, document, 'script', 'https://www.google-analytics.com/analytics.js', 'ga');

ga('create', 'UA-39373211-1', 'auto');
ga('send', 'pageview');

var graphToResize = document.getElementById('count_graph');

var resized = false;

window.onload = function(){
  setTimeout(function() {
    var graphToResize = document.getElementById('count_graph');
    graphToResize.on('plotly_afterplot', triggerResize, {once: true})
  }, 1000);
}


var triggerResize = function() {
    if (!resized) {
      setTimeout(function() {
        window.dispatchEvent(new Event('resize'));
        console.log("fired resize");
      }, 1000);
      resized = true;
  }
}
//
// window.onload = function(){
//   triggerResize();
// }

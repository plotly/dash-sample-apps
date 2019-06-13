var graphToResize = document.getElementById("count_graph");

var resized = false;

window.onload = function() {
  setTimeout(function() {
    var graphToResize = document.getElementById("count_graph");
    graphToResize.on("plotly_afterplot", triggerResize, { once: true });
  }, 1000);
  a;
};

var triggerResize = function() {
  if (!resized) {
    setTimeout(function() {
      window.dispatchEvent(new Event("resize"));
      console.log("fired resize");
    }, 1000);
    resized = true;
  }
};

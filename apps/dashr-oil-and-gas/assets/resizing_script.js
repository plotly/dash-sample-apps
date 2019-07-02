if (!window.dash_clientside) {
  window.dash_clientside = {};
}
window.dash_clientside.clientside = {
  resize: function(value) {
    console.log("resizing..."); // for testing
    setTimeout(function() {
      window.dispatchEvent(new Event("resize"));
      console.log("fired resize");
    }, 1000);
    return null;
  }
};

// var graphToResize = document.getElementById("count_graph");
// var resized = false;
// window.onload = function() {
//   setTimeout(function() {
//     var graphToResize = document.getElementById("count_graph");
//     graphToResize.on("plotly_afterplot", triggerResize, { once: true });
//   }, 1000);
// };
// var triggerResize = function() {
//   if (!resized) {
//     setTimeout(function() {
//       window.dispatchEvent(new Event("resize"));
//       console.log("fired resize");
//     }, 1000);
//     resized = true;
//   }
// };

// if (!window.dash_clientside) {
//   window.dash_clientside = {};
// }
// window.dash_clientside.clientside = {
//   resize: function(value) {
//     console.log("resizing..."); // for testing
//     window.dispatchEvent(new Event("resize"));
//     return null;
//   }
// };

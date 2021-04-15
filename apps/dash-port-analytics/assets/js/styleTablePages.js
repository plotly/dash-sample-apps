// work around with removing bottom border at table current-page(through css styling it is impossible)
document.addEventListener('DOMContentLoaded', function () {
  let timeInterval = setInterval(changeCurrentPageStyle, 400);

  function changeCurrentPageStyle() {
    const element = document.querySelector('.current-page');
    if (element) {
      element.setAttribute('style', 'border:none !important');
      //clearInterval(timeInterval);
    }
  }
});
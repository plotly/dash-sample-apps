function updateScroll(){
    var element = document.getElementById("worker-output");
    element.scrollTop = element.scrollHeight;
}


setInterval(updateScroll,5000);
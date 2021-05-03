window.dash_clientside = Object.assign({}, window.dash_clientside, {
    clientside: {
        get_window_width: function() {
            return window.innerWidth;
        },

        get_window_height: function() {
            return window.innerHeight;
        },

        move_hover: function() {
            document.onmousemove = handleMouseMove;
            function handleMouseMove(event) {
                var eventDoc, doc, body;

                event = event || window.event; // IE-ism

                if (event.pageX == null && event.clientX != null) {
                    eventDoc = (event.target && event.target.ownerDocument) || document;
                    doc = eventDoc.documentElement;
                    body = eventDoc.body;

                    event.pageX = event.clientX +
                      (doc && doc.scrollLeft || body && body.scrollLeft || 0) -
                      (doc && doc.clientLeft || body && body.clientLeft || 0);
                    event.pageY = event.clientY +
                      (doc && doc.scrollTop  || body && body.scrollTop  || 0) -
                      (doc && doc.clientTop  || body && body.clientTop  || 0 );
                }
                //Use event.PageX ...
                document.getElementById('hovered_location').style.top = event.pageY+5+"px";
                document.getElementById('hovered_location').style.left = event.pageX+5+"px";

            }
            return "";
        }
    }
});

window.dash_clientside = Object.assign({}, window.dash_clientside, {
    clientside: {
        canvas_toggle: (n, c) => c === '#509188' ? '#FF7070' : '#509188',
        datasets_params_store: (...rest) => rest[0] > 0 && Object.values(rest).splice(1)
        // I actually want it to have a cache right after initialization.
    }
});

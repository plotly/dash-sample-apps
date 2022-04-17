window.dash_clientside = Object.assign({}, window.dash_clientside, {
    clientside: {
        canvas_toggle: (n, c) => c === '#509188' ? '#FF7070' : '#509188',
        datasets_params_store: (...rest) => rest[0] > 0 ? Object.values(rest).splice(1) : window.dash_clientside.no_update,
        open_offcanvas: (n, sn, is_open) => n || sn ? !is_open : is_open,
        btn_disabled: tabs_table => tabs_table.some(Boolean),
        disable_param_degree: kernel => kernel != 'poly',
        reset_threshold: (n_clicks, fig) => {
            if (n_clicks) {
                let Z = fig['data'][0]['z'].flat(Infinity);
                return -Math.min(...Z) / (Math.max(...Z) - Math.min(...Z));
            } else {
                return 0.5;
            }
        },
        kernel_formula: kernel => ({
            'rbf': '$K(x, z) = exp(-\\gamma||x-z||^2)$',
            'linear': '$K(x, z) = x \\bullet z$',
            'poly': '$K(x,z) = (\\gamma x \\bullet z+r)^d$',
            'sigmoid': '$K(x,z) = tanh(\\gamma x \\bullet z+r)$'
        }[kernel]),
        disable_param_gamma: kernel => { let _ = ['rbf', 'poly', 'sigmoid'].includes(kernel); return [_, _]; },
        scale_param: power => {
            let labels = {};
            for (i of Array.from(Array(5).keys(), n => 2 * n + 1)) {
                labels[i] = power < 0 ? (i / 10 ** -power).toString() : (i * 10 ** power).toString()
            };
            return labels;
        }
    }
});

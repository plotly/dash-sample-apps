import dash_core_components as dcc
import dash_html_components as html

### Operating things
operating_slider_components = [
    html.Div(id="alpha_slider_output"),
    dcc.Slider(id="alpha_slider_input", min=-20, max=20, step=1e-6, value=8.755,),
    html.Div(id="height_slider_output"),
    dcc.Slider(id="height_slider_input", min=0, max=1, step=1e-6, value=0),
    dcc.Checklist(
        id="operating_checklist",
        options=[{"label": " Ground Effect", "value": "ground_effect"}],
        value=[],
    ),
    html.Div(id="streamline_density_slider_output"),
    dcc.Slider(id="streamline_density_slider_input", min=0.1, max=5, step=0.1, value=1),
]

### Shape things
n_kulfan_inputs_per_side = 3
sides = ["Upper", "Lower"]

kulfan_slider_components = []
for side in sides:
    kulfan_slider_components.append(dcc.Markdown(f"##### {side} Surface"))
    for i in range(n_kulfan_inputs_per_side):

        slider_id = f"kulfan_{side.lower()}_{i}"
        slider_min = None
        slider_max = None
        slider_init = None

        if side == "Upper":
            slider_min = -0.3
            slider_max = 0.7
            slider_init = 0.2
            if i == 0:
                slider_min = 0.06
            if n_kulfan_inputs_per_side == 3:
                slider_init = [0.254094, 0.47475, 0.023816][i]
        else:
            slider_min = -0.7
            slider_max = 0.3
            slider_init = -0.2
            if i == 0:
                slider_max = -0.06
            if n_kulfan_inputs_per_side == 3:
                slider_init = [-0.107972, 0.061474, -0.055904][i]

        kulfan_slider_components.extend(
            [
                html.Div(id=slider_id + "_output"),
                dcc.Slider(
                    id=slider_id + "_input",
                    min=slider_min,
                    max=slider_max,
                    step=1e-6,  # continuous
                    value=slider_init,
                ),
            ]
        )

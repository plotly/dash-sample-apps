from dash.testing.application_runners import import_app


# Basic test for the component rendering.
# The dash_duo pytest fixture is installed with dash (v1.0+)
def test_render_component(dash_duo):
    # Start a dash app contained as the variable `app` in `usage.py`
    app = import_app('usage')
    dash_duo.start_server(app)

    # Get the generated component input with selenium
    # The html input will be a children of the #input dash component
    my_component = dash_duo.find_element('#input > input')

    assert 'my-value' == my_component.get_attribute('value')

    # Clear the input
    dash_duo.clear_input(my_component)

    # Send keys to the custom input.
    my_component.send_keys('Hello dash')

    # Wait for the text to equal, if after the timeout (default 10 seconds)
    # the text is not equal it will fail the test.
    dash_duo.wait_for_text_to_equal('#output', 'You have entered Hello dash')

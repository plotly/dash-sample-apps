import dash
import os


def new_dash_app(f, assets_in_cwd=True, verbose=False, **dash_args):
    """
    f is the __file__ variable from the calling script
    returns a new dash app with name = '__main__'
    if assets_in_cwd True, the assets folder is made relative to where the
    command was invoked from.
    """
    assets_folder = "assets"
    if assets_in_cwd:
        script_dir = os.path.dirname(os.path.realpath(f))
        invoke_dir_rel = os.path.relpath(os.getcwd(), script_dir)
        assets_folder = os.path.join(invoke_dir_rel, "assets")
        if verbose:
            print("script_dir", script_dir)
            print("invoke_dir_rel", invoke_dir_rel)
            print("assets should be:", assets_folder)
    ret = dash.Dash(assets_folder=assets_folder, **dash_args)
    return ret

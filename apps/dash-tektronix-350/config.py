import os

# Replace with the name of your Dash app
# This will end up being part of the URL of your deployed app,
# so it can't contain any spaces, capitalizations, or special characters
#
# This name MUST match the name that you specified in the
# Dash App Manager
DASH_APP_NAME = 'dash-daq-tektronix350'
# Set to 'private' if you want to add a login screen to your app
# You can choose who can view the app in your list of files
# at <your-plotly-server>/organize.
# Set to 'public' if you want your app to be accessible to
# anyone who has access to your Plotly server on your network without
# a login screen.
# Set to 'secret' if you want to add a login screen, but allow it
# to be bypassed by using a secret "share_key" parameter.
DASH_APP_PRIVACY = 'public'
# Dash On-Premise is configured with either "Path based routing"
# or "Domain based routing"
# Ask your server administrator which version was set up.
# If a separate subdomain was created,
# then set this to `False`. If it was not, set this to 'True'.
# Path based routing is the default option and most On-Premise
# users use this option.
PATH_BASED_ROUTING = True
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This section only needs to be filled out                                #
# if DASH_APP_PRIVACY is set to 'private' or 'secret'                     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Fill in with your Plotly On-Premise username
os.environ['PLOTLY_USERNAME'] = 'daq-developers'
# Fill in with your Plotly On-Premise API key
# See <your-plotly-server>/settings/api to generate a key
# If you have already created a key and saved it on your own machine
# (from the Plotly-Python library instructions at https://plot.ly/python/getting-started)
# then you can view that key in your ~/.plotly/.config file
# or inside a Python session with these commands:
# import plotly
# print(plotly.tools.get_config_file())
os.environ['PLOTLY_API_KEY'] = 'your-plotly-api-key'
# Fill in with your Plotly On-Premise domain
os.environ['PLOTLY_DOMAIN'] = 'https://your-plotly-domain.com'
os.environ['PLOTLY_API_DOMAIN'] = os.environ['PLOTLY_DOMAIN']
# Fill in with the domain of your Dash subdomain.
# This matches the domain of the Dash App Manager
PLOTLY_DASH_DOMAIN = 'https://your-dash-manager-plotly-domain.com'
# Keep as True if your SSL certificates are valid.
# If you are just trialing Plotly On-Premise with self signed certificates,
# then you can set this to False. Note that self-signed certificates are not
# safe for production.
os.environ['PLOTLY_SSL_VERIFICATION'] = 'True'

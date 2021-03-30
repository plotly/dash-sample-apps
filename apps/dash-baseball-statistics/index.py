# import dash-core, dash-html, dash io, bootstrap
import os

import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

# Dash Bootstrap components
import dash_bootstrap_components as dbc

# Navbar, layouts, custom callbacks
from navbar import Navbar
from layouts import (
    appMenu,
    menuSlider,
    playerMenu,
    teamLayout,
    battingLayout,
    fieldingLayout,
)
import callbacks

# Import app
from app import app

# Import server for deployment
from app import srv as server


app_name = os.getenv("DASH_APP_PATH", "/dash-baseball-statistics")

# Layout variables, navbar, header, content, and container
nav = Navbar()

header = dbc.Row(
    dbc.Col(
        html.Div(
            [
                html.H2(children="Major League Baseball History"),
                html.H3(children="A Visualization of Historical Data"),
            ]
        )
    ),
    className="banner",
)

content = html.Div([dcc.Location(id="url"), html.Div(id="page-content")])

container = dbc.Container([header, content])


# Menu callback, set and return
# Declair function  that connects other pages with content to container
@app.callback(Output("page-content", "children"), [Input("url", "pathname")])
def display_page(pathname):
    if pathname in [app_name, app_name + "/"]:
        return html.Div(
            [
                dcc.Markdown(
                    """
            ### The Applicaiton
            This application is a portfolio project built by [Matt Parra](https://devparra.github.io/) using Plotly's Dash,
            faculty.ai's Dash Bootstrap Components, and Pandas. Using historical MLB (Major League Baseball) data,
            this application provides visualizations for team and player statistics dating from 1903 to 2020. Selecting
            from a dropdown menu, the era will update the list of available teams and players in the range set on the years
            slider. The slider allows the user to adjust the range of years with which the data is presented.

            ### The Analysis
            The applicaiton breaks down each baseballs teams win/loss performance within a range of the teams history.
            Additionally, the application will break down the batting performance with the team batting average, BABIP, and strikeout
            rate. The application also brakes down the piching perfomance using the teams ERA and strikeout to walk ratio. Finally the feilding
            performance of each team is illustrated with total errors and double plays. The applicaiton will also breakdown
            each of teams players statistics within the given era.

            ### The Data
            The data used in this application was retrieved from [Seanlahman.com](http://www.seanlahman.com/baseball-archive/statistics/).
            Provided by [Chadwick Baseball Bureau's GitHub](https://github.com/chadwickbureau/baseballdatabank/) .
            This database is copyright 1996-2021 by Sean Lahman. This data is licensed under a Creative Commons Attribution-ShareAlike
            3.0 Unported License. For details see: [CreativeCommons](http://creativecommons.org/licenses/by-sa/3.0/)
        """
                )
            ],
            className="home",
        )
    elif pathname.endswith("/team"):
        return appMenu, menuSlider, teamLayout
    elif pathname.endswith("/player"):
        return appMenu, menuSlider, playerMenu, battingLayout
    elif pathname.endswith("/field"):
        return appMenu, menuSlider, playerMenu, fieldingLayout
    else:
        return "ERROR 404: Page not found!"


# Main index function that will call and return all layout variables
def index():
    layout = html.Div([nav, container])
    return layout


# Set layout to index function
app.layout = index()

# Call app server
if __name__ == "__main__":
    # set debug to false when deploying app
    app.run_server(debug=True)

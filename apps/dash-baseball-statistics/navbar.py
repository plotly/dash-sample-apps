# Import Bootstrap from Dash
import os

import dash_bootstrap_components as dbc


app_name = os.getenv("DASH_APP_PATH", "/dash-baseball-statistics")

# Navigation Bar fucntion
def Navbar():
    navbar = dbc.NavbarSimple(
        children=[
            dbc.NavItem(dbc.NavLink("Team Analysis", href=f"{app_name}/team")),
            dbc.NavItem(dbc.NavLink("Batting Analysis", href=f"{app_name}/player")),
            dbc.NavItem(
                dbc.NavLink("Pitching/Fielding Analysis", href=f"{app_name}/field")
            ),
        ],
        brand="Home",
        brand_href=f"{app_name}",
        sticky="top",
        color="light",
        dark=False,
        expand="lg",
    )
    return navbar

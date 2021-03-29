# Import Bootstrap from Dash
import dash_bootstrap_components as dbc


# Navigation Bar fucntion
def Navbar():
    navbar = dbc.NavbarSimple(
        children=[
            dbc.NavItem(dbc.NavLink("Team Analysis", href="/team")),
            dbc.NavItem(dbc.NavLink("Batting Analysis", href="/player")),
            dbc.NavItem(dbc.NavLink("Pitching/Fielding Analysis", href="/field")),
        ],
        brand="Home",
        brand_href="/",
        sticky="top",
        color="light",
        dark=False,
        expand="lg",
    )
    return navbar

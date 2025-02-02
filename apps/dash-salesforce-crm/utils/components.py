from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc
from constants import states
from datetime import date


def Header(app):
    name = [
        html.Span("CRM App"),
        html.Span(" using Salesforce API", style={"font-size": "1.8rem", "margin-top": "15px"}),
    ]
    title = html.H2(name, style={"margin-top": 5})
    logo = html.Img(src=app.get_asset_url("images/plotly-logo.png"), style={"float": "right", "height": 60})
    link = html.A(logo, href="https://plotly.com/dash/", target="_blank") 
    demo_link = html.A("ENTERPRISE DEMO", href="https://plotly.com/get-demo/", target="_blank", className="demo-button")
    return dbc.Row([dbc.Col(title, md=8), dbc.Col([demo_link, link], md=4, className="header-logos")], className="header")

def dbc_indicator(text, id_value, width):
    return dbc.Col(
        dbc.Card([
            html.H1(id=id_value, className="card-title"),
            html.P(text),
        ], className="align-items-center"), width = width
    )

def dbc_card(header, child_id, width, table=None):
    if table is None:
        child = dcc.Graph(id=child_id, config=dict(displayModeBar=False))
    else:
        child = dash_table.DataTable(id=child_id, style_table={'overflowX': 'auto'},)
    return dbc.Col(
        dbc.Card(
            dbc.CardBody([
                html.H4(header, className="card-title"),
                child
            ])
        ), width = width
    )

### Leads page components ### ### ### ###


leads_controls = [
        dbc.Col([
            dcc.Dropdown(
                id="converted_leads_dropdown",
                options=[
                    {"label": "By day", "value": "D"},
                    {"label": "By week", "value": "W-MON"},
                    {"label": "By month", "value": "M"},
                ],
                value="M",
                clearable=False,
                searchable=False,
            )
        ], width=2),
        dbc.Col([
            dcc.Dropdown(
                id="lead_source_dropdown",
                options=[
                    {"label": "All status", "value": "all"},
                    {"label": "Open leads", "value": "open"},
                    {"label": "Converted leads", "value": "converted"},
                    {"label": "Lost leads", "value": "lost"},
                ],
                value="all",
                clearable=False,
                searchable=False,
            )
        ], width=2),
        dbc.Col(width=6),
        dbc.Col([
            dbc.Button("Add New Lead", color="primary", size="lg", id="new_lead"),
        ], width=2)
    ]
leads_data_cards = [
    dbc_indicator("Converted Leads", "left_leads_indicator", 4),
    dbc_indicator("Open Leads", "middle_leads_indicator", 4),
    dbc_indicator("Conversion Rates", "right_leads_indicator", 4),
]

leads_graphs = [
    dbc_card("Leads Count per State", "leads_map", width=4),
    dbc_card("Leads by Source", "lead_source", width=4),
    dbc_card("Converted Leads Count", "converted_leads",  width=4),
    dbc_card("Table of Leads", "leads_table", width=12, table=True),
]


def leads_modal():
    return dbc.Modal([
        dbc.ModalHeader(dbc.ModalTitle("New Lead")),
        dbc.ModalBody(dbc.Row([
            dbc.Col([
                dbc.Label("Company Name"),
                dcc.Input(
                    id="new_lead_company",
                    placeholder="Enter company name",
                    type="text",
                    value="",
                    style={"width": "100%"},
                ),
            ], width=6),
            dbc.Col([
                dbc.Label("Company State"),
                dcc.Dropdown(
                    id="new_lead_state",
                    options=states,
                    value="NY",
                    searchable=False,
                ),
            ], width=6),
            dbc.Col([
                dbc.Label("Status"),
                dcc.Dropdown(
                    id="new_lead_status",
                    options=["Open - Not Contacted", "Working - Contacted", "Closed - Converted", "Closed - Not Converted"],
                    value="Open - Not Contacted",
                    searchable=False,
                ),
            ], width=6),
            dbc.Col([
                dbc.Label("Source"),
                dcc.Dropdown(
                    id="new_lead_source",
                    options=["Web", "Phone Inquiry", "Partner Referral", "Purchased List", "Other"],
                    value="Web",
                    searchable=False,
                ),
            ], width=6),
        ])),
        dbc.ModalFooter(dbc.Button("Submit", id="submit_new_lead", className="ms-auto")),
    ],
    id="leads_modal",
    is_open=False,
)



### Cases page components ### ### ### ###

cases_controls = [
        dbc.Col([
            dcc.Dropdown(
                id="cases_period_dropdown",
                options=[
                    {"label": "By day", "value": "D"},
                    {"label": "By week", "value": "W-MON"},
                    {"label": "By month", "value": "M"},
                ],
                value="M",
                clearable=False,
                searchable=False,
            ),
        ], width=2),
        dbc.Col([
            dcc.Dropdown(
                id="priority_dropdown",
                options=[
                    {"label": "All priority", "value": "all_p"},
                    {"label": "High priority", "value": "High"},
                    {"label": "Medium priority", "value": "Medium"},
                    {"label": "Low priority", "value": "Low"},
                ],
                value="all_p",
                clearable=False,
                searchable=False,
            ),
        ], width=2),
        dbc.Col([
            dcc.Dropdown(
                id="origin_dropdown",
                options=[
                    {"label": "All origins", "value": "all"},
                    {"label": "Phone", "value": "Phone"},
                    {"label": "Web", "value": "Web"},
                    {"label": "Email", "value": "Email"},
                ],
                value="all",
                clearable=False,
                searchable=False,
            ),
        ], width=2),
        dbc.Col(width=4),
        dbc.Col([
            dbc.Button("Add New Case", color="primary", size="lg", id="new_case"),
        ], width=2)
    ]
cases_data_cards = [
    dbc_indicator("Low priority cases", "left_cases_indicator", 4),
    dbc_indicator("Medium priority cases", "middle_cases_indicator", 4),
    dbc_indicator("High priority cases", "right_cases_indicator", 4),
]

cases_graphs = [
    dbc_card("Cases Type", "cases_types", width=6),
    dbc_card("Cases Reasons", "cases_reasons", width=6),
    dbc_card("Cases over Time", "cases_by_period",  width=6),
    dbc_card("Cases by Company", "cases_by_account", width=6),
]

def cases_modal(salesforce_manager):
    accounts = salesforce_manager.get_accounts()
    contacts = salesforce_manager.get_contacts()

    contacts["Name"] = (
        contacts["Salutation"]
        + " "
        + contacts["FirstName"]
        + " "
        + contacts["LastName"]
    )
    return dbc.Modal([
        dbc.ModalHeader(dbc.ModalTitle("New Case")),
        dbc.ModalBody(dbc.Row([
            dbc.Col([
                dbc.Label("Account Name"),
                dcc.Dropdown(
                    id="new_case_account",
                    options=[
                        {
                            "label": row["Name"],
                            "value": row["Id"],
                        }
                        for index, row in accounts.iterrows()
                    ],
                    clearable=False,
                    searchable=False,
                    value=accounts.iloc[0].Id,
                )
            ], width=6),
            dbc.Col([
                dbc.Label("Contact Name"),
                dcc.Dropdown(
                    id="new_case_contact",
                    options=[
                        {
                            "label": row["Name"],
                            "value": row["Id"],
                        }
                        for index, row in contacts.iterrows()
                    ],
                    clearable=False,
                    searchable=False,
                    value=contacts.iloc[0].Id,
                )
            ], width=6),
            dbc.Col([
                dbc.Label("Priority"),
                dcc.Dropdown(
                    id="new_case_priority",
                    options=[
                        {"label": "High", "value": "High"},
                        {"label": "Medium", "value": "Medium"},
                        {"label": "Low", "value": "Low"},
                    ],
                    value="Medium",
                    clearable=False,
                    searchable=False,
                ),
            ], width=6),
            dbc.Col([
                dbc.Label("Type"),
                dcc.Dropdown(
                    id="new_case_type",
                    options=[
                        {
                            "label": "Electrical",
                            "value": "Electrical",
                        },
                        {
                            "label": "Mechanical",
                            "value": "Mechanical",
                        },
                        {
                            "label": "Electronic",
                            "value": "Electronic",
                        },
                        {
                            "label": "Structural",
                            "value": "Structural",
                        },
                        {"label": "Other", "value": "Other"},
                    ],
                    value="Electrical",
                    searchable=False,
                ),
            ], width=6),
            dbc.Col([
                dbc.Label("Origin"),
                dcc.Dropdown(
                    id="new_case_origin",
                    options=[
                        {"label": "Phone", "value": "Phone"},
                        {"label": "Web", "value": "Web"},
                        {"label": "Email", "value": "Email"},
                    ],
                    value="Phone",
                    clearable=False,
                    searchable=False,
                ),
            ], width=6),
            dbc.Col([
                dbc.Label("Status"),
                dcc.Dropdown(
                    id="new_case_status",
                    options=[
                        {"label": "New", "value": "New"},
                        {
                            "label": "Working",
                            "value": "Working",
                        },
                        {
                            "label": "Escalated",
                            "value": "Escalated",
                        },
                        {"label": "Closed", "value": "Closed"},
                    ],
                    value="New",
                    searchable=False,
                ),
            ], width=6),
            dbc.Col([
                dbc.Label("Reason"),
                dcc.Dropdown(
                    id="new_case_reason",
                    options=[
                        {
                            "label": "Installation",
                            "value": "Installation",
                        },
                        {
                            "label": "Equipment Complexity",
                            "value": "Equipment Complexity",
                        },
                        {
                            "label": "Performance",
                            "value": "Performance",
                        },
                        {
                            "label": "Breakdown",
                            "value": "Breakdown",
                        },
                        {
                            "label": "Equipment Design",
                            "value": "Equipment Design",
                        },
                        {
                            "label": "Feedback",
                            "value": "Feedback",
                        },
                        {"label": "Other", "value": "Other"},
                    ],
                    value="Installation",
                    clearable=False,
                    searchable=False,
                ),
            ], width=6),
            dbc.Col([
                dbc.Label("Supplied Email"),
                dcc.Input(
                    id="new_case_email",
                    placeholder="email",
                    type="email",
                    value="",
                    style={"width": "100%"},
                ),
            ], width=6),
            dbc.Col([
                dbc.Label("Subject"),
                dcc.Input(
                    id="new_case_subject",
                    placeholder="The Subject of the case",
                    type="text",
                    value="",
                    style={"width": "100%"},
                ),
            ], width=6),
            dbc.Col([
                dbc.Label("Description"),
                dcc.Textarea(
                    id="new_case_description",
                    placeholder="Description of the case",
                    value="",
                    style={"width": "100%"},
                ),
            ], width=6),
        ])),
        dbc.ModalFooter(dbc.Button("Submit", id="submit_new_case", className="ms-auto")),
    ],
    id="cases_modal",
    is_open=False,
)


### Opportunities page components ### ### ### ###

opportunities_controls = [
        dbc.Col([
            dcc.Dropdown(
                id="converted_opportunities_dropdown",
                options=[
                    {"label": "By day", "value": "D"},
                    {"label": "By week", "value": "W-MON"},
                    {"label": "By month", "value": "M"},
                ],
                value="M",
                clearable=False,
                searchable=False,
            ),
        ], width=2),
        dbc.Col([
            dcc.Dropdown(
                id="heatmap_dropdown",
                options=[
                    {"label": "All stages", "value": "all_s"},
                    {"label": "Cold stages", "value": "cold"},
                    {"label": "Warm stages", "value": "warm"},
                    {"label": "Hot stages", "value": "hot"},
                ],
                value="all_s",
                clearable=False,
                searchable=False,
            ),
        ], width=2),
        dbc.Col([
            dcc.Dropdown(
                id="source_dropdown",
                options=[
                    {"label": "All sources", "value": "all_s"},
                    {"label": "Web", "value": "Web"},
                    {"label": "Word of Mouth", "value": "Word of mouth"},
                    {"label": "Phone Inquiry", "value": "Phone Inquiry"},
                    {"label": "Partner Referral", "value": "Partner Referral"},
                    {"label": "Purchased List", "value": "Purchased List"},
                    {"label": "Other", "value": "Other"},
                ],
                value="all_s",
                clearable=False,
                searchable=False,
            ),
        ], width=2),
        dbc.Col(width=4),
        dbc.Col([
            dbc.Button("Add New Opportunity", color="primary", size="lg", id="new_opportunity"),
        ], width=2)
    ]
opportunities_data_cards = [
    dbc_indicator("Won opportunities", "left_opportunities_indicator", 4),
    dbc_indicator("Open opportunities", "middle_opportunities_indicator", 4),
    dbc_indicator("Lost opportunities", "right_opportunities_indicator", 4),
]

opportunities_graphs = [
    dbc_card("Converted Opportunities count", "converted_count", width=4),
    dbc_card("Probabilty heatmap per Stage and Type", "opportunities_heatmap", width=8),
    dbc_card("Top Open opportunities", "top_open_opportunities",  width=6, table=True),
    dbc_card("Top Lost opportunities", "top_lost_opportunities", width=6, table=True),
]


def opportunities_modal():
    return dbc.Modal([
        dbc.ModalHeader(dbc.ModalTitle("New Opportunity")),
        dbc.ModalBody(dbc.Row([
            dbc.Col([
                dbc.Label("Name"),
                dcc.Input(
                    id="new_opportunity_name",
                    placeholder="Name of the opportunity",
                    type="text",
                    value="",
                    style={"width": "100%"},
                ),
            ], width=6),
            dbc.Col([
                dbc.Label("StageName"),
                dcc.Dropdown(
                    id="new_opportunity_stage",
                    options=[
                        {
                            "label": "Prospecting",
                            "value": "Prospecting",
                        },
                        {
                            "label": "Qualification",
                            "value": "Qualification",
                        },
                        {
                            "label": "Needs Analysis",
                            "value": "Needs Analysis",
                        },
                        {
                            "label": "Value Proposition",
                            "value": "Value Proposition",
                        },
                        {
                            "label": "Id. Decision Makers",
                            "value": "Closed",
                        },
                        {
                            "label": "Perception Analysis",
                            "value": "Perception Analysis",
                        },
                        {
                            "label": "Proposal/Price Quote",
                            "value": "Proposal/Price Quote",
                        },
                        {
                            "label": "Negotiation/Review",
                            "value": "Negotiation/Review",
                        },
                        {
                            "label": "Closed/Won",
                            "value": "Closed Won",
                        },
                        {
                            "label": "Closed/Lost",
                            "value": "Closed Lost",
                        },
                    ],
                    clearable=False,
                    value="Prospecting",
                    searchable=False,
                )
            ], width=6),
            dbc.Col([
                dbc.Label("Source"),
                dcc.Dropdown(
                    id="new_opportunity_source",
                    options=["Web", "Phone Inquiry", "Partner Referral", "Purchased List", "Other"],
                    value="Web",
                    searchable=False,
                ),
            ], width=6),
            dbc.Col([
                dbc.Label("Close Date"),
                html.Div(
                    dcc.DatePickerSingle(
                        id="new_opportunity_date",
                        min_date_allowed=date.today(),
                        # max_date_allowed=dt(2017, 9, 19),
                        initial_visible_month=date.today(),
                        date=date.today(),
                    ),
                    style={"textAlign": "left"},
                ),
            ], width=6),
            dbc.Col([
                dbc.Label("Type"),
                dcc.Dropdown(
                    id="new_opportunity_type",
                    options=["Existing Customer - Replacement", "New Customer", "Existing Customer - Upgrade", "Existing Customer - Downgrade"],
                    value="New Customer",
                    searchable=False,
                ),
            ], width=6),
            dbc.Col([
                dbc.Label("Amount"),
                dcc.Input(
                    id="new_opportunity_amount",
                    placeholder="0",
                    type="number",
                    value="",
                    style={"width": "100%"},
                ),
            ], width=6),
            dbc.Col([
                dbc.Label("Probability"),
                dcc.Input(
                    id="new_opportunity_probability",
                    placeholder="0",
                    type="number",
                    max=100,
                    step=1,
                    value="",
                    style={"width": "100%"},
                ),
            ], width=6),
        ])),
        dbc.ModalFooter(dbc.Button("Submit", id="submit_new_opportunity", className="ms-auto")),
    ],
    id="opportunities_modal",
    is_open=False,
)
    
    
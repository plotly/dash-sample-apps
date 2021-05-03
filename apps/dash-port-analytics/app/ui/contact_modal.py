import dash_bootstrap_components as dbc
import dash_html_components as html
from config import strings


def make_contact_modal() -> html.Div:
    """
    Makes a modal that opens on the "Request a port" button in the sidebar.

    :return: HTML div
    """
    return html.Div(
        [
            dbc.Modal(
                id="modal-contact",
                children=[
                    dbc.ModalHeader(strings.MODAL_HEADER),
                    dbc.ModalBody(
                        children=[
                            html.H3(strings.MODAL_BODY_HEADING),
                            html.P(strings.MODAL_BODY_TEXT),
                            html.Div(
                                className="modal-contact-contact-row",
                                children=[
                                    html.I(className="fas fa-phone"),
                                    html.P(strings.MODAL_PHONE),
                                ],
                            ),
                            html.Div(
                                className="modal-contact-contact-row",
                                children=[
                                    html.I(className="fas fa-envelope"),
                                    html.P(strings.MODAL_EMAIL),
                                ],
                            ),
                        ]
                    ),
                    dbc.ModalFooter(
                        children=[
                            html.A(
                                className="a-btn",
                                href="mailto:hello@appsilon.com",
                                children=[strings.MODAL_BTN_MAIL],
                            ),
                            dbc.Button(
                                strings.MODAL_BTN_CANCEL,
                                id="btn-modal-close",
                                className="ml-auto",
                            ),
                        ]
                    ),
                ],
            )
        ]
    )

import os
from textwrap import dedent


def StaticUrlPath(resource):
    if not os.path.exists(os.path.join("static", resource)):
        raise Exception(
            dedent(
                """
            The file "{}" does not exist in the "static" folder.
        """.format(
                    resource, resource
                )
            )
        )
    if "DYNO" in os.environ:
        return "/{}/static/{}".format(os.environ["DASH_APP_NAME"], resource)

    return "/static/{}".format(resource)

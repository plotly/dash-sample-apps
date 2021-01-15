from os import environ


def get_env(name, default=None, conv=None, check_if_none=False):
    try:
        ret = environ[name]
    except KeyError:
        ret = default
        if check_if_none and ret is None:
            raise Exception("Specify " + name + ".")
        return ret
    if conv is not None:
        return conv(ret)
    return ret

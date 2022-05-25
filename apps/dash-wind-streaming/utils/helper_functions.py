import datetime as dt

def get_current_time():
    """Helper function to get the current time in seconds."""

    now = dt.datetime.now()
    total_time = (now.hour * 3600) + (now.minute * 60) + (now.second)
    return total_time
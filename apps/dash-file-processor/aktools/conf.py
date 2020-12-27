import random
import time
import glob
import os
import sys
import collections
from collections import Counter

# Default list of hosts
hosts = ["Hannibal", "Samatchi", "Hanny", "Steeve", "Mustafa"]

# Default host
def_host = "Hannibal"

# Default backward processing time (s)
dbp_time = 10

# Default number of genrerated lines per log file
nb_lines = 6

# log out of order time. Useful for including log files created slightly before init_datetime but still have entry lines in the interval.
log_ofo_time = 9

# dashboard variables
connected_hosts, received_hosts, active_hosts = Counter(), Counter(), Counter()

# dashboard variables
dashboard_refresh=2000 #2s

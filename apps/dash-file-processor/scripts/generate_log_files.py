from aktools.conf import *
from aktools.loggen import *


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Please call with nb of lines parameter")
        sys.exit()
    elif int(sys.argv[1]) > 10000:
        print("Please enter a number of lines less than 10000")
        sys.exit()
    else:
        generate_log_files(int(sys.argv[1]))

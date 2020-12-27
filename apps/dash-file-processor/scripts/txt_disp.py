from aktools.conf import *
from aktools.tools import *


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Please call with two parameter")
        sys.exit()

    elif sys.argv[1] not in hosts:
        print("Please enter an existing host")
        sys.exit()

    else:
        process_log_files(sys.argv[1], int(sys.argv[2]))

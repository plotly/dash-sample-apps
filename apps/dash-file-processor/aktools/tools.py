from aktools.conf import *


# LIST OF HOSTNAMES CONNECTED TO A HOST FOR A GIVEN PERIOD


def connected_hostnames(logpath, init_datetime, end_datetime, Hostname):
    connected_hosts = []

    input_log = open(logpath)
    for line in input_log:
        # Check if within the interval
        if init_datetime <= int(line.split()[0]) <= end_datetime:
            # check if the host initialized the connection and append the receiver if true.
            if line.split()[1] == Hostname:
                connected_hosts.append(line.split()[2])
            # check if the host received the connection and append the initializer if true.
            elif line.split()[2] == Hostname:
                connected_hosts.append(line.split()[1])
        # Exit early: finish the process if the interval is exceeded
        elif int(line.split()[0]) > end_datetime:
            break
    input_log.close()
    return collections.Counter(connected_hosts)


# HOSTNAMES THAT INITIALIZED CONNECTION TO A GIVEN HOST FOR A GIVEN PERIOD


def connected_to(logpath, init_datetime, end_datetime, Hostname):
    hostnames = []
    input_log = open(logpath)
    for line in reversed(list(input_log)):
        # print(''.join(['parsed line: ',line]))     #For debug
        # Check if within the interval and correct hotsname
        if (
            int(line.split()[0]) >= init_datetime
            and int(line.split()[0]) <= end_datetime
            and line.split()[2] == Hostname
        ):
            # print(''.join(['----> considered line: ',line]))  #For debug
            hostnames.append(line.split()[1])
        # Exit early: finish the process if the interval is exceeded
        if int(line.split()[0]) < init_datetime:
            break

    # print('------------------ \n\n')
    input_log.close()
    return collections.Counter(hostnames)


# HOSTNAMES THAT RECEIVED CONNECTION FROM A GIVEN HOST FOR A GIVEN PERIOD


def received_from(logpath, init_datetime, end_datetime, Hostname):
    hostnames = []
    input_log = open(logpath)
    for line in reversed(list(input_log)):
        # print(''.join(['parsed line: ',line]))  #For debug
        # Check if within the interval and correct hotsname
        if (
            int(line.split()[0]) >= init_datetime
            and int(line.split()[0]) <= end_datetime
            and line.split()[1] == Hostname
        ):
            # print(''.join(['----> considered line: ',line]))   #For debug
            hostnames.append(line.split()[2])

        if int(line.split()[0]) < init_datetime:
            break
    # print('------------------ \n\n')
    input_log.close()
    return collections.Counter(hostnames)


#  HOSTNAMES THAT GENERATED THE HIGHEST NUMBER OF CONNECTIONS DURING A PERIOD


def generated_conn(logpath, init_datetime, end_datetime):
    hostnames = []
    input_log = open(logpath)
    for line in reversed(list(input_log)):
        # print(''.join(['parsed line: ',line]))    #For debug
        # Check if within the interval
        if (
            int(line.split()[0]) >= init_datetime
            and int(line.split()[0]) <= end_datetime
        ):
            # print(''.join(['----> considered line: ',line]))   #For debug
            hostnames.append(line.split()[1])

        if int(line.split()[0]) < init_datetime:
            break

    # print('------------------ \n\n')
    input_log.close()
    return collections.Counter(hostnames)


# LOG FILES PROCESSING TOOL


def process_log_files(Hostname, past_time):

    # can achieve the same effect slightly faster by using while 1. This is a single jump operation, as it is a numerical comparison.
    while 1:
        connected_hosts, received_hosts, active_hosts = Counter(), Counter(), Counter()

        # Format unix timestamp to the given one
        init_datetime = int((time.time() - past_time) * 1000)
        end_datetime = int(time.time() * 1000)

        # Filter logfiles by the interval and sort  them by datatime
        past_files = sorted(
            [
                filename
                for filename in glob.glob("output/*.txt")
                if os.path.getmtime(filename) >= init_datetime / 1000 - log_ofo_time
            ],
            key=os.path.getmtime,
        )[::-1]

        # For each file call necessary tools and append the results
        for filename in past_files:
            connected_hosts += connected_to(
                filename, init_datetime, end_datetime, Hostname
            )
            received_hosts += received_from(
                filename, init_datetime, end_datetime, Hostname
            )
            active_hosts += generated_conn(filename, init_datetime, end_datetime)

        # take the hostname with the highest number of occurences and considerif others have the same amount
        active_hosts = [
            h
            for h in active_hosts.most_common()
            if h[1] == active_hosts.most_common(1)[0][1]
        ]

        # Display
        print(
            " ".join(
                [
                    "Hosts that connected to ",
                    Hostname,
                    "in the last",
                    str(past_time),
                    "s are: ",
                    str(dict(connected_hosts)),
                    "\n",
                ]
            )
        )
        print(
            " ".join(
                [
                    "Hosts that received connection from",
                    Hostname,
                    "in the last",
                    str(past_time),
                    "s are: ",
                    str(dict(received_hosts)),
                    "\n",
                ]
            )
        )
        print(
            " ".join(
                [
                    "The hostname that generated most connections in the last",
                    str(past_time),
                    "s is: ",
                    str(dict(active_hosts)),
                    "\n",
                ]
            )
        )

        print("--------------------------------\n\n")

        print(
            "".join(
                [
                    "It is :  ",
                    time.strftime("%X %x"),
                    ".  the next output is in ",
                    str(past_time),
                    " s. \n",
                ]
            )
        )
        # Sleep time
        time.sleep(past_time)

from aktools.conf import *


def generate_log_files(nb_log_lines):
    log_count = 1
    while 1:
        file = open("output/" + "log_" + str(log_count) + ".txt", "a")

        # Loop over and create log entries
        for i in range(0, nb_log_lines):
            timeout = time.time() + 60 * 5  # 5 min timeout
            time.sleep(1)  # 1s sleep
            source, destination = random.sample(hosts, 2)
            file.write(
                str(int(time.time() * 1000)) + " " + source + " " + destination + "\n"
            )
            file.flush()

            if (
                time.time() > timeout
            ):  # break if this file creation takes more than 5 min
                break

        print("log file " + str(log_count) + " finished")
        log_count += 1

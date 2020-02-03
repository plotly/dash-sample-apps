from random import randint

#######################################################################################################################
# Setup
#######################################################################################################################
# Satellite H45-K1 data
file_h_0 = open("non_gps_data_h_0.csv", "w")
file_m_0 = open("non_gps_data_m_0.csv", "w")

# Satellite L12-5 data
file_h_1 = open("non_gps_data_h_1.csv", "w")
file_m_1 = open("non_gps_data_m_1.csv", "w")

# Initialize the first column
file_h_0.write("elevation,temperature,speed,fuel,battery\n")
file_m_0.write("elevation,temperature,speed,fuel,battery\n")

file_h_1.write("elevation,temperature,speed,fuel,battery\n")
file_m_1.write("elevation,temperature,speed,fuel,battery\n")


# Data randomizer function
def randomize(start, variance):
    max_diff = int(start * variance)
    sign = 1 - 2 * randint(0, 1)
    return start + sign * randint(0, max_diff)


# Add data to file
def add_data(data, file):
    file.write(str(data) + ",")
    return data


######################################################################################################################
# Data Generation
#######################################################################################################################

# Data updated by the minute (H45-K1 data)
for i in range(60):
    add_data(randomize(650, 0.2), file_h_0)
    add_data(randomize(290, 0.2), file_h_0)
    add_data((randomize(280, 0.2) / 10), file_h_0)
    add_data(randomize(80, 0.05), file_h_0)
    file_h_0.write(str(randomize(80, 0.05)) + "\n")

# Data updated by the second (H45-K1 data)
for i in range(60):
    add_data(randomize(650, 0.025), file_m_0)
    add_data(randomize(290, 0.025), file_m_0)
    add_data((randomize(280, 0.05) / 10), file_m_0)
    add_data(randomize(80, 0.02), file_m_0)
    file_m_0.write(str(randomize(80, 0.02)) + "\n")

# Data updated by the minute (L12-5 data)
for i in range(60):
    add_data(randomize(750, 0.2), file_h_1)
    add_data(randomize(350, 0.2), file_h_1)
    add_data((randomize(280, 0.2) / 10), file_h_1)
    add_data(randomize(80, 0.05), file_h_1)
    file_h_1.write(str(randomize(80, 0.05)) + "\n")

# Data updated by the second (L12-5 data)
for i in range(60):
    add_data(randomize(650, 0.025), file_m_1)
    add_data(randomize(290, 0.025), file_m_1)
    add_data((randomize(330, 0.05) / 10), file_m_1)
    add_data(randomize(60, 0.02), file_m_1)
    file_m_1.write(str(randomize(65, 0.02)) + "\n")

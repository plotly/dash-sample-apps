from aktools.conf import os, nb_lines, def_host, dbp_time


# Run the file generator in the background
os.system("python setup.py install")

# Run the file generator in the background
os.system("python scripts/generate_log_files.py " + str(nb_lines) + " &")


# Run the processing tools: either the text or the graphic one
# txt display
os.system("python scripts/txt_disp.py " + str(def_host) + " " + str(dbp_time) + " &")
# graphic display
os.system("python scripts/dashboard_disp.py")

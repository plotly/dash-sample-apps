import os
import sys


app_index_string = """app.index_string = \'\'\'
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
        <script>(function(w,d,s,l,i){w[l]=w[l]||[];w[l].push({'gtm.start':new Date().getTime(),event:'gtm.js'});var f=d.getElementsByTagName(s)[0],j=d.createElement(s),dl=l!='dataLayer'?'&l='+l:'';j.async=true;j.src='https://www.googletagmanager.com/gtm.js?id='+i+dl;f.parentNode.insertBefore(j,f);})(window,document,'script','dataLayer','GTM-N6T2RXG');</script>
    </head>
    <body>
    <!-- Google Tag Manager (noscript) -->
        <noscript><iframe src="https://www.googletagmanager.com/ns.html?id=GTM-N6T2RXG"
        height="0" width="0" style="display:none;visibility:hidden"></iframe></noscript>
    <!-- End Google Tag Manager (noscript) -->
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>
\'\'\'\n\n"""

app_file_name = ""

# get the name of the file used to run the app;
# this is not always 'app.py'

with open(os.path.join("/app", "Procfile"), "r") as f:
    contents = f.read().split(" ")
    for item in contents:
        if "server" in item:
            app_file_name = item.split(":")[0]

if app_file_name == "":
    print("Could not find app file")
    sys.exit(1)

full_app_path = os.path.join("/app", app_file_name + ".py")

lines = []

# find the line with the conditional used to run
# the app server; anything after the `app.run_server`
# call will not get executed because the app is running
with open(full_app_path, "r") as f:
    lines = f.readlines()
    name_main_index = 0
    for line in lines:
        if '__name__ == "__main__"' in line:
            name_main_index = lines.index(line)
            break

    # insert the index_string declaration just above
    # the run_server conditional
    lines.insert(name_main_index, app_index_string)

# if something has gone wrong, don't overwrite
# the file; otherwise, write the new lines into it
# (which include the tags)
if len(lines) > 0:
    with open(full_app_path, "w") as f:
        f.writelines(lines)

import requests
import os


deck_demos = [n for n in sorted(os.listdir("./demos")) if ".py" in n]

for demo in deck_demos:
    url = f"https://raw.githubusercontent.com/plotly/dash-deck/master/demos/{demo}"

    r = requests.get(url, allow_redirects=True)

    if r.status_code != 200:
        print("FAILED WITH", url)
        print("Stopping now")
        break

    with open(os.path.join("demos", demo), "wb") as f:
        f.write(r.content)

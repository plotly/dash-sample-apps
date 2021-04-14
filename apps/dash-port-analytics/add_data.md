## How to add Data to the app
1. Inside the `/data` folder locate the CSV file of interest - there are 3 CSV files for separate tabs + 1 CSV file - **add your data there**
2. Under `config/constants.py` add a latitude and lognitude pair for your location
3. Under `config/strings.py` add a variable name for your city - `CITY_<cityname>`
4. Under `app/helpers.py` change `get_port_dropdown_values()` and `get_lat_long_for_port()` functions accordingly
5. Launch the app - `python main.py`
6. Enjoy!
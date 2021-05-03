import pandas as pd
import pickle
import os

from definitions import *
from datamanager.processing import (
    prepare_unloc,
    prepare_ccodes,
    check_ccode_loccode,
    prepare_world_cities,
    enrich_cities,
)


rawdatainfo = pd.read_csv(joinpath(DATAMGR_PATH, "overview_rawdata.csv"), sep=";")

# country codes
if not os.path.isfile(CCODE_FILE_CLEAN_ABS):
    country_codes = prepare_ccodes(ccode_file=CCODE_FILE_ABS)
    with open(CCODE_FILE_CLEAN_ABS, "wb") as f:
        pickle.dump(country_codes, f)
else:
    with open(CCODE_FILE_CLEAN_ABS, "rb") as f:
        country_codes = pickle.load(f)
    print(country_codes.shape)


# city codes
if not os.path.isfile(UNLOC_FILE_CLEAN_ABS):
    city_codes = prepare_unloc(unloc_files=UNLOC_FILES_ABS)
    check_ccode_loccode(country_codes, city_codes)

    cities = prepare_world_cities(world_city_file=WORLD_CITIES_FILE_ABS)
    final_cities = enrich_cities(
        cities=cities, country_codes=country_codes, city_codes=city_codes
    )

    with open(UNLOC_FILE_CLEAN_ABS, "wb") as f:
        pickle.dump(final_cities, f)
else:
    with open(UNLOC_FILE_CLEAN_ABS, "rb") as f:
        final_cities = pickle.load(f)
    print(final_cities.shape)

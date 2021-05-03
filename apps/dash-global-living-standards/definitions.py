import os

joinpath = os.path.join

BASE_PATH = os.path.dirname(os.path.realpath(__file__))
print(f"set base dir to {BASE_PATH}")

DATAMGR_PATH = joinpath(BASE_PATH, "datamanager/")

DATA_BASE_PATH = joinpath(BASE_PATH, "data/")
RAWDATA_PATH = joinpath(DATA_BASE_PATH, "raw_data/")

if not os.path.exists(RAWDATA_PATH):
    raise Exception(
        f"Rawdata path not found. Please create dir and move rawdata there. {os.path.relpath(RAWDATA_PATH)}"
    )


COMPUTED_DATA_PATH = joinpath(DATA_BASE_PATH, "computed_data")
if not os.path.exists(COMPUTED_DATA_PATH):
    os.mkdir(COMPUTED_DATA_PATH)

# city codes
UNLOC_PATH = joinpath(RAWDATA_PATH, "country_city_codes/")

UNLOC_FILES_ABS = [
    joinpath(UNLOC_PATH, "CodeListPart1.csv"),
    joinpath(UNLOC_PATH, "CodeListPart2.csv"),
    joinpath(UNLOC_PATH, "CodeListPart3.csv"),
]

for file in UNLOC_FILES_ABS:
    if not os.path.isfile(file):
        raise Exception(
            f"{os.path.relpath(file)} is missing. Please download from"
            f" https://unece.org/trade/cefact/UNLOCODE-Download or https://datahub.io/core/un-locode "
            f"or define location."
        )

UNLOC_FILE_CLEAN = "city_codes.pickle"
UNLOC_FILE_CLEAN_ABS = joinpath(COMPUTED_DATA_PATH, UNLOC_FILE_CLEAN)

# ISO Country codes
CCODE_FILE_ABS = joinpath(UNLOC_PATH, "ISO3166alpha2CountryCode.csv")

if not os.path.isfile(CCODE_FILE_ABS):
    raise Exception(
        f"{os.path.relpath(CCODE_FILE_ABS)} is missing. Please download from"
        f" https://datahub.io/core/country-list or define location."
    )

CCODE_FILE_CLEAN_ABS = joinpath(COMPUTED_DATA_PATH, "country_codes.pickle")

#
DIMS_BASE_PATH = joinpath(RAWDATA_PATH, "dimensions/")


# world cities
WORLD_CITIES_FILE_ABS = joinpath(DIMS_BASE_PATH, "worldcities_clean.csv")

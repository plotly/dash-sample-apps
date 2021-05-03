import pandas as pd
import string


def filter_rows(df, condition, reason):
    """
    :param reason:
    :param df:
    :param condition: boolean, true for row to keep
    :return: filter country_city_codes df
    """
    n_dropped = (condition == False).sum()
    print(
        f"\nexcluding {n_dropped} locations ({n_dropped / df.shape[0]:.1%}) due to {reason}"
    )
    return df[condition]


def clean_unloc(unloc):
    colidx_to_keep = [1, 2, 3, 4, 7, 10]
    unloc = unloc.iloc[:, colidx_to_keep].copy()
    colnames_dict = dict(
        zip(
            colidx_to_keep,
            ["ccode", "loccode", "locname", "subdiv", "status", "geocode"],
        )
    )
    unloc.rename(columns=colnames_dict, inplace=True)
    print(f"head after import:\n{unloc.head()}")

    status_to_exclude = ["RQ", "RR", "QQ", "UR", "XX"]
    filter_status = ~unloc.status.isin(status_to_exclude)

    unloc = filter_rows(unloc, condition=filter_status, reason="invalid status")

    print(f"\nNAs:\n{unloc.isna().sum()}")
    unloc = filter_rows(
        unloc, condition=(~unloc.status.isna()), reason="missing status"
    )
    unloc = filter_rows(
        unloc, condition=(~unloc.loccode.isna()), reason="missing loccode"
    )
    unloc = filter_rows(unloc, condition=(~unloc.ccode.isna()), reason="missing ccode")
    unloc = filter_rows(
        unloc, condition=(~unloc.ccode.isin(["XZ"])), reason="undefined ccode XZ"
    )

    unloc.insert(2, "citycode", unloc.ccode + unloc.loccode)

    return unloc


def prepare_unloc(unloc_files):
    if not isinstance(unloc_files, list):
        unloc_files = [unloc_files]
    city_codes = pd.DataFrame()
    for file_path in unloc_files:
        print(file_path)
        unloc_raw = pd.read_csv(file_path, encoding="ISO-8859-1", header=None)
        unloc_tmp = clean_unloc(unloc_raw)
        city_codes = pd.concat([city_codes, unloc_tmp])

    # cleaning
    city_codes = city_codes[["ccode", "locname", "loccode"]]
    city_codes = city_codes.rename(
        columns={"ccode": "country_code", "locname": "city", "loccode": "city_code"}
    )
    city_codes = (
        city_codes.groupby(["country_code", "city"]).city_code.first().reset_index()
    )

    return city_codes


def prepare_ccodes(ccode_file):
    country = pd.read_csv(ccode_file, names=["cname", "ccode"], header=0)
    filter_cname = ~country.cname.isin(["Namibia", "Bouvet Island"])
    country = country[filter_cname]

    country.rename(columns={"cname": "country", "ccode": "country_code"}, inplace=True)

    country.loc[country.country_code == "KR", "country"] = "South Korea"
    country.loc[country.country_code == "RU", "country"] = "Russia"
    country.loc[country.country_code == "IR", "country"] = "Iran"
    country.loc[country.country_code == "VN", "country"] = "Vietnam"
    country.loc[country.country_code == "TZ", "country"] = "Tanzania"
    country.loc[country.country_code == "BO", "country"] = "Bolivia"
    country.loc[country.country_code == "TW", "country"] = "Taiwan"
    country.loc[country.country_code == "VE", "country"] = "Venezuela"
    country.loc[country.country_code == "SY", "country"] = "Syria"
    country.loc[country.country_code == "CZ", "country"] = "Czechia"
    country.loc[country.country_code == "LA", "country"] = "Laos"
    country.loc[country.country_code == "MD", "country"] = "Moldova"
    country.loc[country.country_code == "MK", "country"] = "Macedonia"
    country.loc[country.country_code == "TW", "country"] = "Taiwan"
    # country = country.append({'country': 'Namibia', 'country_code': 'NA'}, ignore_index=True)
    country.loc[country.country_code == "BN", "country"] = "Brunei"
    country.loc[country.country_code == "FM", "country"] = "Micronesia"
    country.loc[country.country_code == "RE", "country"] = "Reunion"
    country.loc[country.country_code == "VA", "country"] = "Vatican City"
    country.loc[country.country_code == "FK", "country"] = "Falkland Islands"

    obs_per_country = country.groupby("country").size()
    obs_per_country
    if (obs_per_country > 1).sum() != 0:
        raise Exception("Country name not unique!")

    country

    return country


def check_ccode_loccode(ccodes, unloc):

    unloc_check = pd.DataFrame(
        {"country_code": unloc.country_code.unique(), "country_city_codes": 1}
    )
    # print(unloc_check.isna().sum())
    ccode_loccode_check = ccodes.merge(unloc_check, how="outer", on="country_code")

    # ccode_loccode_check[ccode_loccode_check.isnull().any(axis=1)]
    # unloc_check[unloc_check.ccode.str.contains('NA')]
    # country_city_codes[country_city_codes.ccode=='XZ']
    # ccodes[ccodes.cname=='Namibia']

    if ccode_loccode_check.isna().sum().sum() != 0:
        raise Exception(
            "Missmatch between ccodes and country_city_codes. Investigate both files."
        )
    else:
        return "ok"


def filter_per_group(df, id_col, filter_col, filter_fct):
    idx = df.groupby(id_col)[filter_col].transform(filter_fct) == df[filter_col]
    return df[idx]


def prepare_world_cities(world_city_file):
    cities = pd.read_csv(world_city_file)
    (cities.city + "_" + cities.country).nunique() / cities.shape[0]

    # cities['city_country'] = cities.city +'_'+ cities.country

    # rows_per_citycountry = cities.groupby('city_country').size()
    # rows_per_citycountry[rows_per_citycountry > 1]

    cities = filter_per_group(
        df=cities, id_col=["city", "country"], filter_col="population", filter_fct=max
    )

    # rename

    # inspection
    # country[country.country.str.contains('falk', case=False)]
    # country.loc[country.ccode == 'KO']

    cities.loc[cities.country == "Korea, South", "country"] = "South Korea"
    cities.loc[cities.country.str.contains("congo", case=False), "country"] = "Congo"
    cities.loc[
        cities.country.str.contains("ivoire", case=False), "country"
    ] = "Côte d'Ivoire"
    cities.loc[
        cities.country.str.contains("bosnia", case=False), "country"
    ] = "Bosnia and Herzegovina"
    cities.loc[
        cities.country.str.contains("bahamas", case=False), "country"
    ] = "Bahamas"
    cities.loc[
        cities.country.str.contains("verde", case=False), "country"
    ] = "Cape Verde"
    cities.loc[
        cities.country.str.contains("tobago", case=False), "country"
    ] = "Trinidad and Tobago"
    cities.loc[cities.country.str.contains("Gambia", case=False), "country"] = "Gambia"
    cities.loc[
        cities.country.str.contains("tome", case=False), "country"
    ] = "Sao Tome and Principe"
    cities.loc[
        cities.country.str.contains("Grenadines", case=False), "country"
    ] = "Saint Vincent and the Grenadines"
    cities.loc[
        cities.country.str.contains("Antigua", case=False), "country"
    ] = "Antigua and Barbuda"
    cities.loc[cities.country.str.contains("Macau", case=False), "country"] = "Macao"
    cities.loc[
        cities.country.str.contains("Micronesia", case=False), "country"
    ] = "Micronesia"
    cities.loc[
        cities.country.str.contains("Reunion", case=False), "country"
    ] = "Reunion"
    cities.loc[
        cities.country.str.contains("Isle of Man", case=False), "country"
    ] = "Isle of Man"
    cities.loc[
        cities.country.str.contains("falkland", case=False), "country"
    ] = "Falkland Islands"

    return cities


def merge_country_codes_on_cities(cities, country):
    # check join country code in cities
    cities_check = pd.DataFrame(
        {"country_name": cities.country.unique(), "in_cities": 1}
    )
    cities_check

    country_check = pd.DataFrame(
        {"country_name": country.country.unique(), "in_countries": 1}
    )

    check = cities_check.merge(country_check, how="left", on="country_name")
    check.in_cities.isna().sum()
    check.in_countries.isna().sum()
    check[check.in_countries.isna()]

    cities = cities.merge(country, how="left", on="country")

    cities_no_countrycode = cities[cities.country_code.isna()].country.unique()
    print(cities.country.nunique())
    print(cities.shape)
    print(len(cities_no_countrycode))
    if len(cities_no_countrycode) == 7:
        cities = cities.loc[
            (~cities.country_code.isna()),
        ]
    else:
        raise Exception("Eliminating more countries than expected")
    print(cities.shape)
    print(cities.country.nunique())

    print(f"\nNAs in cities:\n{cities.isna().sum()}")

    # remove dubplicate cities
    cities = cities[cities.city != "Darhan"]
    cities = cities[cities.city != "Crato"]

    # check if city country unique
    (cities.country + cities.city).nunique()
    obs_per_citycountry = cities.groupby(["country", "city"]).size()
    if len(obs_per_citycountry[obs_per_citycountry > 1]) != 0:
        raise Exception("city country not unique!")

    return cities


def merge_city_codes(cities, city_codes):
    cities_check = cities[["country_code", "city"]].copy()
    cities_check["in_cities"] = 1
    cities_check = cities_check.merge(
        city_codes, how="left", on=["country_code", "city"]
    )

    obs_per_citycountry = cities_check.groupby(["country_code", "city"]).size()
    obs_per_citycountry[obs_per_citycountry > 1]
    print("cities_check", cities_check.shape)
    print("cities", cities.shape)

    # merge city codes
    if len(obs_per_citycountry[obs_per_citycountry > 1]) == 0:
        cities = cities.merge(city_codes, how="left", on=["country_code", "city"])
        print("cities", cities.shape)

    else:
        raise Exception("Can not merge city_codes in cities. Keys ambigiouse")

    cities["cc_code"] = cities.country_code + "-" + cities.city_code

    cities = cities[
        [
            "cc_code",
            "city",
            "city_code",
            "country",
            "country_code",
            "lat",
            "lng",
            "population",
        ]
    ]

    # create artificial city codes
    print("creating artificial city codes")
    cc_codes = cities[~cities.cc_code.isna()].cc_code.to_list()
    for idx, row in cities[cities.city_code.isna()].iterrows():
        city_name_tmp_clean = (
            row.city.translate(str.maketrans("", "", string.punctuation))
            .replace(" ", "")
            .upper()
        )
        city_code_tmp = city_name_tmp_clean[0:3]
        cc_code_tmp = row.country_code + "-" + city_code_tmp
        if cc_code_tmp not in cc_codes and len(city_name_tmp_clean) >= 3:
            cc_codes.append(cc_code_tmp + "-X")
            cities.loc[idx, "cc_code"] = cc_code_tmp + "-X"
            cities.loc[idx, "city_code"] = city_code_tmp + "-X"
        else:
            if len(city_name_tmp_clean) >= 4:
                city_code_tmp = city_name_tmp_clean[0:2] + city_name_tmp_clean[3]
                cc_code_tmp = row.country_code + "-" + city_code_tmp
                if cc_code_tmp not in cc_codes:
                    cc_codes.append(cc_code_tmp + "-X")
                    cities.loc[idx, "cc_code"] = cc_code_tmp + "-X"
                    cities.loc[idx, "city_code"] = city_code_tmp + "-X"

                else:
                    pass
            else:
                pass

    if cities.cc_code.isna().sum() <= 278:
        cities = cities[~cities.cc_code.isna()]
    else:
        raise Exception("More more missing cc_codes than expected!")

    print("done")

    return cities


def enrich_cities(cities, country_codes, city_codes):

    cities = merge_country_codes_on_cities(cities, country_codes)

    cities = merge_city_codes(cities, city_codes)

    return cities


def clean_citytemp():
    print("clean_citytemp")


cleaner = globals()["clean_citytemp"]
cleaner()

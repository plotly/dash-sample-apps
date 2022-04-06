import io
import base64
import json
import numpy as np
import pandas as pd


def parse_contents(contents, filename, header, usecols=None):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        if filename.endswith('csv'):
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')),
                             usecols=usecols)
        elif filename.endswith('xls') or filename.endswith('xlsx'):
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded), usecols=usecols)

    except:
        print('Somthing wrong with uploader.')

    if header:
        return df.columns

    else:
        return df


def handle_json(js):
    df = pd.DataFrame(json.loads(js)['objects'])
    df.fillna(value={'path': ''}, inplace=True)
    df['dot'] = df['path'].apply(lambda x: 1 if len(x) == 2 else 0)
    df['c'] = df['stroke'].apply(lambda x: 0 if x == '#509188' else 1)
    X = df[df['dot'] == 1]['pathOffset'].apply(pd.Series).to_numpy()
    y = df[df['dot'] == 1]['c'].to_numpy()
    X = (X - 250) / 500
    X[:, 1] = -X[:, 1]
    return X, y
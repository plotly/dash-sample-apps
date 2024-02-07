from constants import millnames
import math


# returns most significant part of a number
def millify(n):
    n = float(n)
    millidx = max(
        0,
        min(
            len(millnames) - 1, int(math.floor(0 if n == 0 else math.log10(abs(n)) / 3))
        ),
    )
    return "{:.0f}{}".format(n / 10 ** (3 * millidx), millnames[millidx])



# returns top 5 open opportunities
def top_open_opportunities(df):
    df = df.sort_values("Amount", ascending=True)
    cols = ["CreatedDate", "Name", "Amount", "StageName"]
    df = df[cols].iloc[:5]
    # only display 21 characters
    df["Name"] = df["Name"].apply(lambda x: x[:30])
    return df.to_dict('records'), [{"name": i, "id": i} for i in df.columns]


# returns top 5 lost opportunities
def top_lost_opportunities(df):
    df = df[df["StageName"] == "Closed Lost"]
    cols = ["CreatedDate", "Name", "Amount", "StageName"]
    df = df[cols].sort_values("Amount", ascending=False).iloc[:5]
    # only display 21 characters
    df["Name"] = df["Name"].apply(lambda x: x[:30])
    return df.to_dict('records'), [{"name": i, "id": i} for i in df.columns]


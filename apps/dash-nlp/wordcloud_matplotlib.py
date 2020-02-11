import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt
from PIL import Image
from io import BytesIO
import base64
from os import path, getcwd
from wordcloud import WordCloud, ImageColorGenerator, STOPWORDS


def fig_to_uri(in_fig, close_all=True, **save_args):
    # type: (plt.Figure) -> str
    """
    Save a figure as a URI
    :param in_fig:
    :return:
    """
    out_img = BytesIO()
    in_fig.savefig(out_img, format="png", **save_args)
    if close_all:
        in_fig.clf()
        plt.close("all")
    out_img.seek(0)  # rewind file
    encoded = base64.b64encode(out_img.read()).decode("ascii").replace("\n", "")
    return "data:image/png;base64,{}".format(encoded)


def create_wordcloud(df):
    complaints_text = list(df["Consumer complaint narrative"].dropna().values)

    # join all documents in corpus
    text = " ".join(list(complaints_text))
    print("Complaints received")
    print(len(complaints_text))

    d = getcwd()
    mask = np.array(Image.open(path.join(d, "thumbs-down.png")))

    STOPWORDS.add("XXXX")
    STOPWORDS.add("XX")
    STOPWORDS.add("xx")
    STOPWORDS.add("xxxx")
    # TODO exclude name of all banks here
    STOPWORDS.add("wells")
    STOPWORDS.add("fargo")

    wc = WordCloud(
        background_color="white",
        stopwords=STOPWORDS,
        max_words=1000,
        mask=mask,
        max_font_size=90,
        random_state=42,
        contour_width=1,
        contour_color="#119DFF",
    )
    wc.generate(text)

    # create wordcloud shape from image
    fig = plt.figure(figsize=[8, 8])
    ax = plt.imshow(wc.recolor(), interpolation="bilinear")
    plt.axis("off")
    out_url = fig_to_uri(fig, bbox_inches="tight")
    return out_url

import pathlib
import pandas as pd
from ldacomplaints import lda_analysis
from wordcloud import STOPWORDS
import multiprocessing
import re
import time

DATA_PATH = pathlib.Path(__file__).parent.resolve()
EXTERNAL_STYLESHEETS = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]
FILENAME = "data/customer_complaints_narrative_sample.csv"
GLOBAL_DF = pd.read_csv(DATA_PATH.joinpath(FILENAME), header=0)


def add_stopwords(selected_bank):
    """
    In order to make a more useful NLP-data based graphs, it helps to remove
    common useless words. In this case XXXX usually represents a redacted name
    We also exlude more standard words defined in STOPWORDS which is provided by
    the Wordcloud dash component.
    """
    selected_bank_words = re.findall(r"[\w']+", selected_bank)
    for word in selected_bank_words:
        STOPWORDS.add(word.lower())

    print("Added %s stopwords:" % selected_bank)
    for word in selected_bank_words:
        print("\t", word)
    return STOPWORDS


def precomputer_one_lda(bank):
    print("crunching LDA for: ", bank, "...")
    add_stopwords(bank)
    bank_df = GLOBAL_DF[GLOBAL_DF["Company"] == bank]
    # try:
    tsne_lda, lda_model, topic_num, df_dominant_topic = lda_analysis(
        bank_df, list(STOPWORDS)
    )
    if tsne_lda is None:
        print("crunching LDA for: ", bank, ". DONE (SAD)")
        return {bank: "NOT ENOUGH DOCS/DICTS TO RUN LDA"}
    topic_top3words = [
        (i, topic)
        for i, topics in lda_model.show_topics(formatted=False)
        for j, (topic, wt) in enumerate(topics)
        if j < 3
    ]

    df_top3words_stacked = pd.DataFrame(topic_top3words, columns=["topic_id", "words"])
    df_top3words = df_top3words_stacked.groupby("topic_id").agg(", \n".join)
    df_top3words.reset_index(level=0, inplace=True)

    print(len(tsne_lda))
    print(len(df_dominant_topic))
    tsne_df = pd.DataFrame(
        {
            "tsne_x": tsne_lda[:, 0],
            "tsne_y": tsne_lda[:, 1],
            "topic_num": topic_num,
            "doc_num": df_dominant_topic["Document_No"],
        }
    )

    results_bank = {
        str(bank): {
            "tsne_df": tsne_df.to_json(),
            "df_dominant_topic": df_dominant_topic.to_json(),
        }
    }
    print("crunching LDA for: ", bank, ". DONE")
    return results_bank


def precompute_all_lda():
    """ QD function for precomputing all necessary LDA results
     to allow much faster load times when the app runs. """
    print("precomupting LDA for all banks...")
    start_time = time.time()

    file = open("precomupted", "w")
    file.close()
    failed_banks = []
    counter = 0
    bank_names = GLOBAL_DF["Company"].value_counts().keys().tolist()
    results = {}
    nbr_workers = 4
    print("Making a worker pool with: %d workers" % nbr_workers)
    print("Starting the run at:", start_time)
    pool = multiprocessing.Pool(processes=nbr_workers)
    results = pool.map(precomputer_one_lda, bank_names)
    end_time = time.time()

    print("It took: %s" % end_time - start_time)

    print("Saving data to file...")
    file = open("precomupted", "a")
    for result in results:
        file.write(str(result))
    file.close()
    print("Saving data to file. DONE")
    save_time = time.time()
    print("It took: %s" % save_time - end_time)

    print("DONE")
    print("did %d banks" % counter)
    print("failed %d:" % len(failed_banks))
    for fail in failed_banks:
        print(fail)


if __name__ == "__main__":
    precompute_all_lda()

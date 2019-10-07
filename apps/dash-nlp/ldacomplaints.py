import spacy  # for our NLP processing
import nltk  # to use the stopwords library
import string  # for a list of all punctuation
from nltk.corpus import stopwords  # for a list of stopwords
import gensim
from sklearn.manifold import TSNE
import pandas as pd
import numpy as np

# Now we can load and use spacy to analyse our complaint
nlp = spacy.load("en_core_web_sm")


def format_topics_sentences(ldamodel, corpus, texts):
    sent_topics_df = pd.DataFrame()

    # Get main topic in each document
    for i, row in enumerate(ldamodel[corpus]):
        row = sorted(row, key=lambda x: (x[1]), reverse=True)
        # Get the Dominant topic, Perc Contribution and Keywords for each document
        for j, (topic_num, prop_topic) in enumerate(row):
            if j == 0:  # => dominant topic
                wp = ldamodel.show_topic(topic_num)
                topic_keywords = ", ".join([word for word, prop in wp])
                sent_topics_df = sent_topics_df.append(
                    pd.Series([int(topic_num), round(prop_topic, 4), topic_keywords]),
                    ignore_index=True,
                )
            else:
                break
    sent_topics_df.columns = ["Dominant_Topic", "Perc_Contribution", "Topic_Keywords"]

    # Add original text to the end of the output
    contents = pd.Series(texts)
    sent_topics_df = pd.concat([sent_topics_df, contents], axis=1)
    return sent_topics_df


def lda_analysis(docs, stop_words):
    # stop_words = stopwords.words('english')

    # TODO: fix a custom stop words file
    # stop_words.append('XXXX')
    # stop_words.append('XX')
    # stop_words.append('xx')
    ##stop_words.append('xxxx')
    # stop_words.append('wells')
    # stop_words.append('fargo')
    def cleanup_text(doc):
        doc = nlp(doc, disable=["parser", "ner"])
        tokens = [tok.lemma_.lower().strip() for tok in doc if tok.lemma_ != "-PRON-"]
        tokens = [
            tok for tok in tokens if tok not in stop_words and tok not in punctuations
        ]
        return tokens

    punctuations = string.punctuation

    processed_docs = list(map(cleanup_text, docs))
    print(len(processed_docs))

    dictionary = gensim.corpora.Dictionary(processed_docs)
    dictionary.filter_extremes(no_below=10, no_above=0.95, keep_n=100000)
    bow_corpus = [dictionary.doc2bow(doc) for doc in processed_docs]

    lda_model = gensim.models.LdaModel(
        bow_corpus, num_topics=5, id2word=dictionary, passes=10
    )

    df_topic_sents_keywords = format_topics_sentences(
        ldamodel=lda_model, corpus=bow_corpus, texts=docs
    )

    df_dominant_topic = df_topic_sents_keywords.reset_index()
    df_dominant_topic.columns = [
        "Document_No",
        "Dominant_Topic",
        "Topic_Perc_Contrib",
        "Keywords",
        "Text",
    ]

    topic_num, tsne_lda = tsne_analysis(lda_model, bow_corpus)

    return (tsne_lda, lda_model, topic_num, df_dominant_topic)


def tsne_analysis(ldamodel, corpus):
    topic_weights = []
    for i, row_list in enumerate(ldamodel[corpus]):
        topic_weights.append([w for i, w in row_list])

    # Array of topic weights
    df_topics = pd.DataFrame(topic_weights).fillna(0).values

    # Keep the well separated points (optional)
    # arr = arr[np.amax(arr, axis=1) > 0.35]

    # Dominant topic number in each doc
    topic_nums = np.argmax(df_topics, axis=1)

    # tSNE Dimension Reduction
    tsne_model = TSNE(n_components=2, verbose=1, random_state=0, angle=0.99, init="pca")
    tsne_lda = tsne_model.fit_transform(df_topics)

    return (topic_nums, tsne_lda)

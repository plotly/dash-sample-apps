# ========== (c) JP Hwang 2020-03-17  ==========

import logging

# ===== START LOGGER =====
logger = logging.getLogger(__name__)
root_logger = logging.getLogger()
root_logger.setLevel(logging.INFO)
sh = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
sh.setFormatter(formatter)
root_logger.addHandler(sh)

import pandas as pd
import plotly.express as px
from sklearn.manifold import TSNE

desired_width = 320
pd.set_option("display.max_columns", 20)
pd.set_option("display.width", desired_width)


def main():

    comp_df = pd.read_csv("data/customer_complaints_narrative_sample.csv", index_col=0)

    # ===============================================
    # ========== EXPLORATORY DATA ANALYSIS ==========
    # ===============================================
    # Complaints by date
    fig = px.histogram(
        comp_df, x="datetime", template="plotly_white", title="Complaint counts by date"
    )
    fig.update_xaxes(categoryorder="category descending", title="Date").update_yaxes(
        title="Number of complaints"
    )
    fig.update_layout(width=1200, height=500)
    fig.show()

    # Complaints by words
    fig = px.histogram(
        comp_df,
        x="Words_clipped",
        template="plotly_white",
        title="Complain counts by length",
    )
    fig.update_xaxes(
        categoryorder="total descending",
        title="Number of words (clipped at 1000 words)",
    ).update_yaxes(title="Number of complaints")
    fig.update_layout(width=1200, height=500)
    fig.show()

    # There are too many companies to display - let's just show a few companies
    # ========== Top n companies only ==========
    # Pre-processing data
    top_comps = (
        comp_df.groupby("Company")["Date received"]
        .count()
        .sort_values(ascending=False)[:10]
        .index
    )
    top_comps_df = comp_df[comp_df["Company"].isin(top_comps)]

    # Top companies by complaints
    fig = px.histogram(
        top_comps_df,
        x="Company",
        template="plotly_white",
        title="Complaint counts by company",
    )
    fig.update_xaxes(categoryorder="total descending").update_yaxes(
        title="Number of complaints"
    )
    fig.update_layout(width=1200, height=500)
    fig.show()

    # Complaints by company & date
    fig = px.histogram(
        top_comps_df,
        x="datetime",
        template="plotly_white",
        title="Complaint counts by date & company",
        color="Company",
        nbins=6,
        log_y=True,
        barmode="group",
    )
    fig.update_xaxes(categoryorder="category descending", title="Date").update_yaxes(
        title="Number of complaints"
    )
    fig.update_layout(width=1200, height=500)
    fig.show()

    # =====================================================
    # ========== PLOT N-GRAM RELATED DATA HERE ==========
    # =====================================================

    bigram_df = pd.read_csv("data/bigram_data.csv", index_col=0)
    fig = px.bar(
        bigram_df[:20],
        x="ngram",
        y="count",
        title="Counts of top bigrams",
        template="plotly_white",
        labels={"ngram": "Bigram", "count": "Count"},
    )
    fig.update_layout(width=1200, height=500)
    fig.show()

    # Hierarchical Treemap
    fig = px.treemap(
        top_comps_df,
        title="Treemap chart by companies and whether complaint mentions credit report",
        path=["Company", "credit_report"],
        color="Words",
        color_continuous_scale=px.colors.sequential.GnBu,
    )
    fig.update_layout(width=1200, height=600)
    fig.show()

    # Visualising proportions
    comp_grp_df = pd.read_csv("data/comp_bigram_data.csv", index_col=0)
    # fig = px.scatter(comp_grp_df, x='bigram', y='company', size='count', color='Words', template='plotly_white',
    #                  labels={'Words':'Length<BR>(words)', 'bigram': 'Bigram', 'company': 'Company'},
    #                  category_orders=top_comps, range_color=[150, 450], color_continuous_scale=px.colors.sequential.YlOrRd)
    # fig.update_traces(marker=dict(line=dict(width=1, color='Gray')))
    # fig.update_layout(width=1200, height=500)
    # fig.show()

    fig = px.bar(
        comp_grp_df,
        x="portion",
        y="company",
        template="plotly_white",
        orientation="h",
        labels={"portion": "% of Complaints", "bigram": "Bigram", "company": "Company"},
        color="bigram",
        color_discrete_sequence=px.colors.qualitative.Safe,
    )
    fig.update_layout(font=dict(size=10, color="DarkSlateGray"))
    fig.update_layout(width=1200, height=500)
    fig.show()

    # =====================================================
    # ========== PLOT N-GRAM RELATED DATA HERE ==========
    # =====================================================

    vects_df = pd.read_csv("data/bigram_vectors.csv", index_col=0)
    embed_df = pd.read_csv(
        "data/tsne_bigram_data.csv", index_col=0
    )  # Bigram embedding dataframe, with placeholder tsne values (at perplexity=3)

    # Try different t-SNE values here
    X_embedded = TSNE(n_components=2, perplexity=3).fit_transform(vects_df)
    embed_df["tsne_1"] = X_embedded[:, 0]
    embed_df["tsne_2"] = X_embedded[:, 1]

    # Plot t-SNE graph
    fig = px.scatter(
        embed_df,
        x="tsne_1",
        y="tsne_2",
        hover_name="bigram",
        text="bigram",
        size="count",
        color="words",
        size_max=45,
        template="plotly_white",
        title="Bigram similarity and frequency",
        labels={"words": "Avg. Length<BR>(words)"},
        color_continuous_scale=px.colors.sequential.Sunsetdark,
    )
    fig.update_traces(marker=dict(line=dict(width=1, color="Gray")))
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    fig.update_layout(width=1200, height=500)
    fig.show()


if __name__ == "__main__":
    main()

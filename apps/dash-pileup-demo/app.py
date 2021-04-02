import os
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output
import dash_bio
import pandas as pd
import numpy as np

from layout_helper import run_standalone_app

text_style = {"color": "#506784", "font-family": "Open Sans"}

_COMPONENT_ID = "pileup-browser"


def description():
    return "An interactive in-browser track viewer."


def header_colors():
    return {
        "bg_color": "#0F5BA7",
        "font_color": "white",
    }


def rna_differential(app):
    basal_bam = {
        "url": app.get_asset_url("data/rna/SRR1552454.fastq.gz.sampled.converted.bam"),
        "indexUrl": app.get_asset_url(
            "data/rna/SRR1552454.fastq.gz.sampled.converted.bam.bai"
        ),
    }

    luminal_bam = {
        "url": app.get_asset_url("data/rna/SRR1552448.fastq.gz.sampled.bam"),
        "indexUrl": app.get_asset_url("data/rna/SRR1552448.fastq.gz.sampled.bam.bai"),
    }

    return {
        "range": {"contig": "chr1", "start": 54986297, "stop": 54991347},
        "tracks": [
            {"viz": "scale", "label": "Scale"},
            {"viz": "location", "label": "Location"},
            {
                "viz": "genes",
                "label": "genes",
                "source": "bigBed",
                "sourceOptions": {
                    "url": app.get_asset_url("data/rna/mm10.chr1.ncbiRefSeq.sorted.bb")
                },
            },
            {
                "viz": "coverage",
                "label": "Basal Mouse cells",
                "source": "bam",
                "sourceOptions": basal_bam,
            },
            {
                "viz": "pileup",
                "vizOptions": {"viewAsPairs": True},
                "label": "Basal Mouse cells",
                "source": "bam",
                "sourceOptions": basal_bam,
            },
            {
                "viz": "coverage",
                "label": "Luminal Mouse cells",
                "source": "bam",
                "sourceOptions": luminal_bam,
            },
            {
                "viz": "pileup",
                "label": "Luminal Mouse cells",
                "source": "bam",
                "sourceOptions": luminal_bam,
            },
        ],
    }


REFERENCE = {
    "label": "mm10",
    "url": "https://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit",
}

DATAPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets/data")

# Differentially expressed genes (identified in R, see assets/data/rna/README.md)
DE_dataframe = pd.read_csv(os.path.join(DATAPATH, "rna", "DE_genes.csv"))
# add SNP column
DE_dataframe["SNP"] = "NA"


def layout(app):
    HOSTED_CASE_DICT = {
        "rna-differential": rna_differential(app),
    }

    HOSTED_USE_CASES = [
        {"value": "rna-differential", "label": "Differential RNA-seq"},
    ]

    return html.Div(
        id="pileup-body",
        className="app-body",
        children=[
            html.Div(
                id="pileup-control-tabs",
                className="control-tabs",
                children=[
                    dcc.Tabs(
                        id="pileup-tabs",
                        value="data",
                        children=[
                            dcc.Tab(
                                label="Data",
                                value="data",
                                children=html.Div(
                                    className="control-tab",
                                    children=[
                                        dcc.Graph(
                                            id="pileup-dashbio-volcanoplot",
                                            figure=dash_bio.VolcanoPlot(
                                                dataframe=DE_dataframe,
                                                effect_size="log2FoldChange",
                                                title="Differentially Expressed Genes",
                                                genomewideline_value=-np.log10(0.05),
                                                p="padj",
                                                snp="SNP",
                                                gene="Gene",
                                            ),
                                        )
                                    ],
                                ),
                            ),
                            dcc.Tab(
                                label="About",
                                value="what-is",
                                children=html.Div(
                                    className="control-tab",
                                    children=[
                                        html.H4(
                                            className="what-is",
                                            children="What is pileup.js?",
                                        ),
                                        dcc.Markdown(
                                            """
                                The Dash pileup.js component is a high-performance genomics
                                data visualization component developed originally by the Hammer Lab
                                (https://github.com/hammerlab/pileup.js). pileup.js
                                supports visualization of genomic file formats, such as vcfs,
                                bam, and bigbed files. pileup.js additionally allows flexible
                                interaction with non-standard data formats. Users can visualize
                                GA4GH JSON formatted alignments, features and variants. Users can
                                also connect with and visualize data stored in GA4GH formatted data
                                stores.
                                """
                                        ),
                                    ],
                                ),
                            ),
                        ],
                    )
                ],
            ),
            dcc.Loading(
                parent_className="dashbio-loading",
                id="pileup-output",
                children=html.Div(
                    [
                        dash_bio.Pileup(
                            id=_COMPONENT_ID,
                            range=HOSTED_CASE_DICT["rna-differential"]["range"],
                            reference=REFERENCE,
                            tracks=HOSTED_CASE_DICT["rna-differential"]["tracks"],
                        )
                    ]
                ),
            ),
        ],
    )


def callbacks(_app):
    HOSTED_CASE_DICT = {
        "rna-differential": rna_differential(_app),
    }

    HOSTED_USE_CASES = [
        {"value": "rna-differential", "label": "Differential RNA-seq"},
    ]

    @_app.callback(
        Output(_COMPONENT_ID, "range"), Input("pileup-dashbio-volcanoplot", "clickData")
    )
    def update_range(point):

        data = HOSTED_CASE_DICT["rna-differential"]

        if point is None:
            range = data["range"]
        else:

            # get genomic location of selected genes and goto
            pointText = point["points"][0]["text"]
            gene = pointText.split("GENE: ")[-1]

            row = DE_dataframe[DE_dataframe["Gene"] == gene].iloc[0]

            range = {"contig": row["chr"], "start": row["start"], "stop": row["end"]}

        return range


app = run_standalone_app(layout, callbacks, header_colors, __file__)
server = app.server

if __name__ == "__main__":
    app.run_server(debug=True, port=8050)

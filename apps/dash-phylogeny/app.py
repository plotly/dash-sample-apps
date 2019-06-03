# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
import numpy as np
import random

from dash.dependencies import Input, Output
from Bio import Phylo


app = dash.Dash(__name__)
server = app.server

virus_name = "measles"
species = ['Avian', 'Ebola', 'Lassa', 'Measles', 'Mumps', 'Zika']
tree_fig = {}
mapbox_access_token = "pk.eyJ1IjoicGxvdGx5bWFwYm94IiwiYSI6ImNqdnBvNDMyaTAxYzkzeW5ubWdpZ2VjbmMifQ.TXcBE-xg9BFdV2ocecc_7g"


def get_x_coordinates(tree):
    """Associates to  each clade an x-coord.
       returns dict {clade: x-coord}
    """
    xcoords = tree.depths()
    # tree.depth() maps tree clades to depths (by branch length).
    # returns a dict {clade: depth} where clade runs over all Clade instances of the tree, and depth
    # is the distance from root to clade

    #  If there are no branch lengths, assign unit branch lengths
    if not max(xcoords.values()):
        xcoords = tree.depths(unit_branch_lengths=True)
    return xcoords


def get_y_coordinates(tree, dist=1.3):
    """
       returns  dict {clade: y-coord}
       The y-coordinates are  (float) multiple of integers (i*dist below)
       dist depends on the number of tree leafs
    """
    maxheight = tree.count_terminals()  # Counts the number of tree leafs.
    # Rows are defined by the tips/leafs
    ycoords = dict((leaf, maxheight - i * dist) for i, leaf in enumerate(reversed(tree.get_terminals())))

    def calc_row(clade):
        for subclade in clade:
            if subclade not in ycoords:
                calc_row(subclade)
        ycoords[clade] = (ycoords[clade.clades[0]] +
                          ycoords[clade.clades[-1]]) / 2

    if tree.root.clades:
        calc_row(tree.root)
    return ycoords


def get_clade_lines(orientation='horizontal', y_curr=0, x_start=0, x_curr=0, y_bot=0, y_top=0,
                    line_color='rgb(25,25,25)', line_width=0.5):
    """define a shape of type 'line', for branch
    """
    branch_line = dict(type='line',
                       layer='below',
                       line=dict(color=line_color,
                                 width=line_width)
                       )
    if orientation == 'horizontal':
        branch_line.update(x0=x_start,
                           y0=y_curr,
                           x1=x_curr,
                           y1=y_curr)
    elif orientation == 'vertical':
        branch_line.update(x0=x_curr,
                           y0=y_bot,
                           x1=x_curr,
                           y1=y_top)
    else:
        raise ValueError("Line type can be 'horizontal' or 'vertical'")

    return branch_line


def draw_clade(clade, x_start, line_shapes, line_color='rgb(15,15,15)', line_width=1, x_coords=0, y_coords=0):
    """Recursively draw the tree branches, down from the given clade"""

    x_curr = x_coords[clade]
    y_curr = y_coords[clade]

    # Draw a horizontal line from start to here
    branch_line = get_clade_lines(orientation='horizontal', y_curr=y_curr, x_start=x_start, x_curr=x_curr,
                                  line_color=line_color, line_width=line_width)

    line_shapes.append(branch_line)

    if clade.clades:
        # Draw a vertical line connecting all children
        y_top = y_coords[clade.clades[0]]
        y_bot = y_coords[clade.clades[-1]]

        line_shapes.append(get_clade_lines(orientation='vertical', x_curr=x_curr, y_bot=y_bot, y_top=y_top,
                                           line_color=line_color, line_width=line_width))

        # Draw descendants
        for child in clade:
            draw_clade(child, x_curr, line_shapes, x_coords=x_coords, y_coords=y_coords)


def read_treefile(filename):
    tree = Phylo.read(filename, "newick")
    return tree


def read_metadata(filename):
    df = pd.read_csv(filename)
    return df


def create_title(virus, nb_genome):
    graph_title = "Phylogeny of " + virus.title() + " Virus<br>" + str(
        nb_genome) + " genomes colored according to region and country"
    return graph_title


def compute_expensive_data(chemin):
    dir = dir + chemin
    return dir


def create_map_bubble_year(virus_name, metadata_file_stat, map_choice, min_date, max_date):
    df = pd.read_csv(metadata_file_stat)
    # To select only the data between min_date and max_date
    df = df[df["Year"] >= min_date]
    df = df[df["Year"] <= max_date]
    # min_date, max_date = min_max_date(df)
    df.head()

    cases = []
    colors = ['rgb(231,229,204)', 'rgb(255,255,204)', 'rgb(255,178,102)', 'rgb(255,153,51)',
              'rgb(204,0,0)', 'rgb(189,215,231)', 'rgb(107,174,214)', 'rgb(33,113,181)',
              'rgb(255,102,255)', 'rgb(189,15,255)', 'rgb(121,74,244)', 'rgb(133,13,181)',
              'rgb(239,243,255)', 'rgb(189,215,231)', 'rgb(107,174,214)', 'rgb(33,113,181)']
    months = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun', 7: 'Jul', 8: 'Aug', 9: 'Sept', 10: 'Oct',
              11: 'Nov', 12: 'Dec'}

    if map_choice == 2:
        df_old = df
        df = pd.read_csv('data/2014_world_gdp_with_codes.csv')
        for index2, row2 in df.iterrows():
            for index, row in df_old.iterrows():
                if row[7] == row2[2]:
                    df.loc[index2, 'GDP (BILLIONS)'] += df_old.loc[index, 'Value']

        data = [dict(
            type='choropleth',
            locations=df['CODE'],
            z=df['GDP (BILLIONS)'],
            text=df['COUNTRY'],
            colorscale=[[0, "rgb(5, 10, 172)"], [0.35, "rgb(40, 60, 190)"], [0.5, "rgb(70, 100, 245)"], \
                        [0.6, "rgb(90, 120, 245)"], [0.7, "rgb(106, 137, 247)"], [1, "rgb(220, 220, 220)"]],
            autocolorscale=False,
            reversescale=True,
            marker=dict(
                line=dict(
                    color='rgb(180,180,180)',
                    width=0.5
                )),
            colorbar=dict(
                autotick=False,
                titleside="right",
                title='Number of ' + virus_name.title() + ' Virus'),
        )]

        layout = dict(
            title=virus_name.title() + ' Virus cases in West Africa <br>Between ' + str(min_date) + " and " + str(
                max_date),
            margin=go.layout.Margin(l=0, r=0, b=0, t=35),
            geo=dict(
                showframe=False,
                showcoastlines=False,
                projection=dict(
                    type='Mercator'
                )
            ),
            autosize=True
        )

        fig = dict(data=data, layout=layout)
        return fig


country_colors = ["red", "green", "blue", "orange", "purple", "pink", "yellow", "black", "gray", "silver", "violet",
                  "yellowgreen", "turquoise", "sienna", "salmon"]

def create_curve_line(df, virus_name, min_date, max_date):

    results_country_by_year = df.groupby(['Country', 'Year'])['Value'].sum()
    # Translate pandas.core.series.series object in dataframe
    df_group_by_country = results_country_by_year.to_frame()
    # Rename the first column in Value
    df_group_by_country.columns = ['Value']
    # Move the index values (i.e. Country column) in column and reset index
    df_group_by_country = df_group_by_country.reset_index()
    country_list = df_group_by_country.Country.unique()
    p_data = []
    p_name = []
    p_year = []


    for country_name in country_list:
        country_data = np.extract(df_group_by_country.Country == country_name, df_group_by_country.Value)
        country_year = np.extract(df_group_by_country.Country == country_name, df_group_by_country.Year)

        if country_year[0] < min_date:
            min_date = country_year[0]
        if country_year[len(country_year)-1] > max_date:
            max_date = country_year[len(country_year)-1]

        # Stock country name in array structure
        if len(country_name) != 0:
            p_name.append(country_name)

        # Stock country data in array structure
        if len(country_data) != 0:
            p_data.append(country_data)

        # Stock country year in array structure
        if len(country_year) != 0:
            p_year.append(country_year)

    # Add step in x axe
    step = 1
    if 5 < max_date - min_date <= 10:
        step = 2
    elif 10 < max_date - min_date <= 50:
        step = 5
    elif max_date - min_date > 50:
        step = 10
    marks_data = []
    for i in range(int(min_date), int(max_date + 1), step):
        marks_data.append(str(i))

    if i < int(max_date):
        marks_data.append(str(max_date))

    year = marks_data

    i = 0
    data = []
    for l_country in p_data:
        trace = go.Scatter(
            x=p_year[i],
            y=l_country,
            name=p_name[i],
            line=dict(
                color=country_colors[i],
                width=1
            )
        )
        i = i + 1
        data.append(trace)

    # Edit the layout
    layout = dict(
        title="Evolution of " + virus_name.title() + " virus <br>Between " + str(min_date) + " and " + str(max_date) + " by country.",
        xaxis=dict(title='Year'),
        yaxis=dict(title='Number of ' + virus_name.title() + " virus"),
        legend=dict(overflowY="scroll"),
        margin=go.layout.Margin(l=50, r=0, b=100, t=100),
    )

    fig = dict(data=data, layout=layout)
    return fig


def create_tree(virus_name, tree_file, metadata_file, ord_by):
    tree = read_treefile(tree_file)
    x_coords = get_x_coordinates(tree)
    y_coords = get_y_coordinates(tree)
    line_shapes = []
    draw_clade(tree.root, 0, line_shapes, line_color='rgb(25,25,25)', line_width=1, x_coords=x_coords,
               y_coords=y_coords)
    my_tree_clades = x_coords.keys()
    X = []
    Y = []
    text = []

    for cl in my_tree_clades:
        X.append(x_coords[cl])
        Y.append(y_coords[cl])
        text.append(cl.name)

    df = read_metadata(metadata_file)
    data_metadata_stat_csv = df.groupby(ord_by)['Strain'].count()

    # for index_val, series_val in data_metadata_stat_csv.iteritems():
    df.columns
    nb_genome = len(df)

    graph_title = create_title(virus_name, nb_genome)
    intermediate_node_color = 'rgb(100,100,100)'

    NA_color = {'Cuba': 'rgb(252, 196, 174)',  # from cm.Reds color 0.2, ... 0.8
                'Dominican Republic': 'rgb(201, 32, 32)',
                'El Salvador': 'rgb(253, 202, 181)',
                'Guadeloupe': 'rgb(253, 202, 181)',
                'Guatemala': 'rgb(252, 190, 167)',
                'Haiti': 'rgb(252, 145, 114)',
                'Honduras': 'rgb(239, 66, 49)',
                'Jamaica': 'rgb(252, 185, 161)',
                'Martinique': 'rgb(252, 190, 167)',
                'Mexico': 'rgb(247, 109, 82)',
                'Nicaragua': 'rgb(249, 121, 92)',
                'Panama': 'rgb(252, 185, 161)',
                'Puerto Rico': 'rgb(252, 174, 148)',
                'Saint Barthelemy': 'rgb(253, 202, 181)',
                'USA': 'rgb(188, 20, 26)',
                'Canada': 'rgb(188, 20, 26)',
                'USVI': 'rgb(206, 36, 34)',
                'Dominica': 'rgb(209, 39, 37)',
                'Trinidad And Tobago': 'rgb(201, 19, 51)',
                'Belize': 'rgb(241, 49, 61)',
                'Grenada': 'rgb(222, 32, 38)',
                'Costa Rica': 'rgb(211, 12, 48)',
                'Bermuda': 'rgb(231, 53, 24)'
                }

    SAmer_color = {'Brazil': 'rgb(21, 127, 59)',  # from cm.Greens colors 0.2, 0.4, 0.6, 0.8
                   'Colombia': 'rgb(153, 213, 149)',
                   'Ecuador': 'rgb(208, 237, 202)',
                   'French Guiana': 'rgb(211, 238, 205)',
                   'Peru': 'rgb(208, 237, 202)',
                   'Suriname': 'rgb(206, 236, 200)',
                   'Venezuela': 'rgb(202, 234, 196)',
                   'Puerto Rico': 'rgb(201, 235, 199)',
                   'Argentina': 'rgb(203, 225, 185)',
                   'Bolivia': 'rgb(217, 197, 165)',
                   'Paraguay': 'rgb(154, 217, 195)',
                   'Chile': 'rgb(231, 85, 168)',
                   'Aruba': 'rgb(191, 35, 128)'
                   }

    SAsia_color = {'Singapore': '#0000EE',
                   'Vietnam': '#1E90FF',
                   'Malaysia': '#1E90AF',
                   'Philippines': '#1E90AE',
                   'Thailand': '#1E90AB',
                   'Myanmar': '#1E90AC',
                   'Cambodia': '#1E90AA',
                   'Indonesia': '#1E90AA',
                   'Brunei': '#1E90BA',
                   'Laos': '#1E90BF'
                   }

    pl_SAsia = [[0.0, '#1E90FF'], [0.5, '#1E90FF'], [0.5, '#0000EE'], [1.0, '#0000EE']]

    Oceania_color = {'American Samoa': 'rgb(209,95,238)',
                     'Fiji': 'rgb(238,130, 238)',
                     'French Polynesia': 'rgb(148,0,211)',
                     'Tonga': 'rgb(238,130, 238)',
                     'Australia': 'rgb(233,125, 235)',
                     'Micronesia': 'rgb(231,123, 235)',
                     'New Caledonia': 'rgb(229,119, 233)',
                     'Marshall Islands': 'rgb(227,117, 231)',
                     'Guam': 'rgb(267,137, 251)',
                     'Papua New Guinea': 'rgb(277,187, 291)',
                     'Solomon Islands': 'rgb(22,167, 251)',
                     'Cook Islands': 'rgb(20,187, 211)',
                     'Samoa': 'rgb(50,127, 221)',
                     'Nauru': 'rgb(34, 92, 98)',
                     'Palau': 'rgb(214, 132, 238)',
                     'Vanuatu': 'rgb(252, 42, 128)',
                     'Niue': 'rgb(272, 52, 158)',
                     'New Zealand': 'rgb(242, 71, 133)'
                     }

    China_color = {'China': 'rgb(255,185,15'}

    JapanKorea_color = {'Japan': '#fcdd04'}

    SubsaharanAfrica_color = {'Guinea': 'rgb(209,95,238)',
                              'Liberia': 'rgb(238,130, 238)',
                              'Sierra Leone': 'rgb(148,0,211)',
                              'Cote D Ivoire': 'rgb(145,0,209)',
                              'Angola': 'rgb(143,0,207)',
                              'Seychelles': 'rgb(145,10,217)',
                              'Comoros': 'rgb(141,5,203)',
                              'Madagascar': 'rgb(233,60, 281)',
                              'Eritrea': 'rgb(202, 122, 118)',
                              'Somalia': 'rgb(115,51,222)',
                              'Djibouti': 'rgb(203,57, 211)',
                              'Burkina Faso': 'rgb(141,21,239)',
                              'Ghana': 'rgb(102,57,232)',
                              'Tanzania': 'rgb(217,37, 291)',
                              'Mozambique': 'rgb(213,17, 231)',
                              'Senegal': 'rgb(231,133, 219)',
                              'Togo': 'rgb(121,21,198)',
                              }

    Africa_color = {'Sudan': 'rgb(209,95,238)',
                    'Gambia': 'rgb(238,130, 238)',
                    'Nigeria': 'rgb(235,135, 233)',
                    'Mali': 'rgb(235,131, 229)',
                    'Senegal': 'rgb(231,133, 219)',
                    'Cote D Ivoire': 'rgb(145,0,209)',
                    'Burkina Faso': 'rgb(141,21,239)',
                    'Seychelles': 'rgb(145,10,217)',
                    'Somalia': 'rgb(115,51,222)',
                    'Ghana': 'rgb(102,57,232)',
                    'Tanzania': 'rgb(217,37, 291)',
                    'Mozambique': 'rgb(213,17, 231)',
                    'Djibouti': 'rgb(203,57, 211)',
                    'Madagascar': 'rgb(233,60, 281)',
                    'Comoros': 'rgb(141,5,203)',
                    'Togo': 'rgb(121,21,198)',
                    'Angola': 'rgb(212, 92, 138)',
                    'Eritrea': 'rgb(202, 122, 118)',
                    'Guinea': 'rgb(209,95,238)',
                    'Sierra Leone': 'rgb(148,0,211)',
                    'Liberia': 'rgb(238,130, 238)',
                    'Tunisia': 'rgb(228,99, 298)',
                    'Cameroon': 'rgb(207,78, 199)',
                    'South Africa': 'rgb(222,22, 222)',
                    'Congo': 'rgb(231,41, 172)',
                    'Algeria': 'rgb(237,35, 168)',
                    'Morocco': 'rgb(223,27, 165)',
                    'Zambia': 'rgb(218,62, 265)',
                    'Kenya': 'rgb(118,2, 215)',
                    'Uganda': 'rgb(128,35, 265)',
                    'Egypt': 'rgb(143,52, 265)',
                    'Ethiopia': 'rgb(206,36, 265)',
                    'Niger': 'rgb(121,52, 187)',
                    'Mayotte': 'rgb(101,32,165)',
                    'Rwanda': 'rgb(325,25,144)',
                    'Gabon': 'rgb(319,7,197)'
                    }

    Europe_color = {'France': 'rgb(209,95,238)',
                    'Germany': 'rgb(238,130, 238)',
                    'Italy': 'rgb(238,130, 238)',
                    'United Kingdom': 'rgb(238,130, 238)',
                    'Netherlands': 'rgb(148,0,211)',
                    'Spain': 'rgb(141,7,221)',
                    'Portugal': 'rgb(139,11,219)',
                    'Ireland': 'rgb(128,15,279)',
                    'Slovakia': 'rgb(121,25,209)',
                    'Romania': 'rgb(171,45,197)',
                    'Sweden': 'rgb(135,96,208)',
                    'Norway': 'rgb(138,56,213)',
                    'Slovenia': 'rgb(138,45,265)',
                    'Denmark': 'rgb(258,25,265)',
                    'Iceland': 'rgb(138,7,185)',
                    'Ukraine': 'rgb(298,65,265)',
                    'Czech Republic': 'rgb(226,96,128)',
                    'Albania': 'rgb(111,20,201)',
                    'Greece': 'rgb(108,63,265)',
                    'Latvia': 'rgb(121,35,299)'
                    }

    country = []
    region = []
    color = [intermediate_node_color] * len(X)

    for k, strain in enumerate(df['Strain']):

        i = text.index(strain)

        # Split journal title if gt 50 characters
        new_title_journal = split_at_n_caracter(df.loc[k, 'Journal'], 50)

        text[i] = text[i] + '<br>Country: ' + '{:s}'.format(df.loc[k, 'Country']) + '<br>Region: ' + '{:s}'.format(
            df.loc[k, 'Region']) + \
                  '<br>Collection date: ' + '{:s}'.format(df.loc[k, 'Date']) + \
                  '<br>Journal: ' + '{:s}'.format(new_title_journal) + '<br>Authors: ' + '{:s}'.format(
            df.loc[k, 'Authors'])
        country.append(df.loc[k, 'Country'])
        region.append(df.loc[k, 'Region'])
        if df.loc[k, 'Region'] == 'North America':
            color[i] = NA_color[df.loc[k, 'Country']]
        elif df.loc[k, 'Region'] == 'South America':
            color[i] = SAmer_color[df.loc[k, 'Country']]
        elif df.loc[k, 'Region'] == 'Southeast Asia':
            color[i] = SAsia_color[df.loc[k, 'Country']]
        elif df.loc[k, 'Region'] == 'Oceania':
            color[i] = Oceania_color[df.loc[k, 'Country']]
        elif df.loc[k, 'Region'] == 'China':
            color[i] = '#fecc00'
        elif df.loc[k, 'Region'] == 'Japan Korea':
            color[i] = '#dc7928'
        elif df.loc[k, 'Region'] == 'Subsaharan Africa':
            color[i] = SubsaharanAfrica_color[df.loc[k, 'Country']]
        elif df.loc[k, 'Region'] == 'Africa':
            color[i] = Africa_color[df.loc[k, 'Country']]
        elif df.loc[k, 'Region'] == 'Europe':
            color[i] = Europe_color[df.loc[k, 'Country']]
        else:
            pass

    axis = dict(showline=False,
                zeroline=False,
                showgrid=False,
                showticklabels=False,
                title=''  # y title
                )

    label_legend = set(list(df[ord_by]))
    nodes = []

    for elt in label_legend:
        node = dict(type='scatter',
                    x=X,
                    y=Y,
                    mode='markers',
                    marker=dict(color=color,
                                size=5),
                    text=text,  # vignet information of each node
                    hoverinfo='',
                    name=elt
                    )
        nodes.append(node)

    layout = dict(title=graph_title,
                  paper_bgcolor='rgba(0,0,0,0)',
                  dragmode="select",
                  font=dict(family='Balto', size=14),
                  # width=1000,
                  height=1000,
                  autosize=True,
                  showlegend=True,
                  xaxis=dict(showline=True,
                             zeroline=False,
                             showgrid=True,  # To visualize the vertical lines
                             ticklen=4,
                             showticklabels=True,
                             title='branch length'),
                  yaxis=axis,
                  hovermode='closest',
                  shapes=line_shapes,
                  plot_bgcolor='rgb(250,250,250)',
                  legend={'x': 0, 'y': 1},
                  )

    fig = dict(data=nodes, layout=layout)
    return fig


def split_at_n_caracter(title, n):
    sentences = "<br>".join([title[i:i + n] for i in range(0, len(title), n)])
    return sentences


# TO DO validation file and directory exist
def create_paths_file(virus_name, level1="", level2="", level3=""):
    dir = "data/" + virus_name + "/"
    if level1 == "" and level2 == "" and level3 == "":
        tree_file = dir + "nextstrain_" + virus_name + "_tree.new"
        metadata_file = dir + "nextstrain_" + virus_name + "_metadata.csv"
        stat_file = dir + "stat_year_nextstrain_" + virus_name + "_metadata.csv"
        return tree_file, metadata_file, stat_file
    elif level2 == "" and level3 == "":
        dir = dir + "/" + level1 + "/"
        tree_file = dir + "nextstrain_" + virus_name + "_" + level1 + "_tree.new"
        metadata_file = dir + "nextstrain_" + virus_name + "_" + level1 + "_metadata.csv"
        stat_file = dir + "stat_year_nextstrain_" + virus_name + "_" + level1 + "_metadata.csv"
        return tree_file, metadata_file, stat_file
    elif level3 == "":
        dir = dir + "/" + level1 + "/" + level2 + "/"
        tree_file = dir + "nextstrain_" + virus_name + "_" + level1 + "_" + level2 + "_tree.new"
        metadata_file = dir + "nextstrain_" + virus_name + "_" + level1 + "_" + level2 + "_metadata.csv"
        stat_file = dir + "stat_year_nextstrain_" + virus_name + "_" + level1 + "_" + level2 + "_metadata.csv"
        return tree_file, metadata_file, stat_file
    else:
        dir = dir + "/" + level1 + "/" + level2 + "/" + level3 + "/"
        tree_file = dir + "nextstrain_" + virus_name + "_" + level1 + "_" + level2 + "_" + level3 + "_tree.new"
        metadata_file = dir + "nextstrain_" + virus_name + "_" + level1 + "_" + level2 + "_" + level3 + "_metadata.csv"
        stat_file = dir + "stat_year_nextstrain_" + virus_name + "_" + level1 + "_" + level2 + "_" + level3 + "_metadata.csv"
        return tree_file, metadata_file, stat_file


def date_4_number(date):
    if 0 <= date <= 18:
        date = date + 2000
    else:
        date = date + 1900
    return date


def slicer(min_date, max_date):
    step = 1
    if 5 < max_date - min_date <= 10:
        step = 2
    elif 10 < max_date - min_date <= 50:
        step = 5
    elif max_date - min_date > 50:
        step = 10
    marks_data = {}
    for i in range(int(min_date), int(max_date) + 1, step):
        if i > int(max_date):
            marks_data[i] = str(int(max_date))
        else:
            marks_data[i] = str(i)

    if i < int(max_date):
        marks_data[int(max_date)] = str(max_date)

    return marks_data


def min_max_date(df):
    min_date = df['Year'].min()
    max_date = df['Year'].max()
    if min_date > max_date:
        tmp = min_date
        min_date = max_date
        max_date = tmp
    return min_date, max_date


tree_file, metadata_file, metadata_file_stat = create_paths_file(virus_name, level1="", level2="", level3="")
# To know the minimum and the maximum values of date for slicer
df_stat_metadata = pd.read_csv(metadata_file_stat)
min_date, max_date = min_max_date(df_stat_metadata)
# create the dictionary of slider
marks_data = slicer(min_date, max_date)
min_max_date_value = [min_date, max_date]

fig = create_tree(virus_name, tree_file, metadata_file, "Country")
tree_fig[tree_file] = fig

fig_map_bubble = create_map_bubble_year(virus_name, metadata_file_stat, 2, min_date, max_date)
fig_curve_line = create_curve_line(df_stat_metadata, virus_name, min_date, max_date)


######################################### MAIN APP #########################################
app.layout = html.Div([
    # Banner display
    html.Div([
        html.H2(
            'Phylogeny trees and global spread of 6 viruses',
            id='title'
        ),
        html.Img(
            src="https://s3-us-west-1.amazonaws.com/plotly-tutorials/logo/new-branding/dash-logo-by-plotly-stripe-inverted.png"
        )
    ],
        className="banner"
    ),

    # Body
    html.Div(
        [
            html.Div(
                [
                    html.Div(
                        [
                            html.Div(
                                [
                                html.H6(children='Dataset'),
                                dcc.Dropdown(
                                    id='d_virus-name',
                                    options=[{'label': species[i], 'value': species[i]} for i in
                                                range(len(species))],
                                    value='Measles',
                                ),
                                html.Div(id='output-container'),

                                html.Div(id='controls-container_mumps', children=[
                                    dcc.Dropdown(
                                        id='d_mumps',
                                        options=[{'label': i, 'value': i} for i in ['global', 'na']],
                                        value='global',
                                    ),
                                ]),

                                html.Div(id='controls-container_dengue', children=[
                                    dcc.Dropdown(
                                        id='d_dengue',
                                        options=[{'label': i, 'value': i} for i in
                                                    ['all', 'denv1', 'denv2', 'denv3', 'denv4']],
                                        value='all',
                                    ),
                                ]),

                                html.Div(id='controls-container_lassa', children=[
                                    dcc.Dropdown(
                                        id='d_lassa',
                                        options=[{'label': i, 'value': i} for i in ['s', 'l']],
                                        value='s',
                                    ),
                                ]),

                                html.Div(id='controls-container_avian', children=[
                                    dcc.Dropdown(
                                        id='d_avian_opt1',
                                        options=[{'label': i, 'value': i} for i in ['h7n9']],
                                        value='h7n9',
                                    ),
                                    dcc.Dropdown(
                                        id='d_avian_opt2',
                                        options=[{'label': i, 'value': i} for i in
                                                    ['ha', 'mp', 'na', 'ns', 'np', 'pa', 'pb2', 'pb1']],
                                        value='ha',
                                    ),
                                ]),

                                html.Div(id='controls-container_flu', children=[
                                    dcc.Dropdown(
                                        id='d_flu_opt1',
                                        options=[{'label': i, 'value': i} for i in
                                                    ['h3n2', 'h1n1pdm', 'vic', 'yam']],
                                        value='h3n2',
                                    ),
                                    dcc.Dropdown(
                                        id='d_flu_opt2',
                                        options=[{'label': i, 'value': i} for i in
                                                    ['ha', 'na']],
                                        value='ha',
                                    ),
                                    dcc.Dropdown(
                                        id='d_flu_opt3',
                                        options=[{'label': i, 'value': i} for i in
                                                    ['2y', '3y', '6y', '12y']],
                                        value='3y',
                                    ),
                                ]),
                                html.H6(children='Data Range'),
                                html.Div(id='id-slicer', children=[
                                    dcc.RangeSlider(
                                        id='id-year',
                                        min=min_date,
                                        max=max_date,
                                        step=1,
                                        marks=marks_data,
                                        value=min_max_date_value
                                    ),
                                ]),
                                html.Div(id='output-container-range-slider'),
                                dcc.Graph(
                                    id='curve-line-graph',
                                    figure=fig_curve_line,
                                    style={'height': 700}
                                ),
                            ]),
                        ],
                        className="four columns",
                        style={'margin-top': '10'}
                    ),
                    html.Div(
                        className="eight columns",
                        style={'margin-top': '10'},
                        children=html.Div([
                            html.Div(id='right-top-graph')
                        ])
                    )
                ], className="row"),

            html.Div(
                [
                    html.Div(
                        [
                            dcc.Graph(
                                id='graph_map',
                                figure=fig_map_bubble
                            )
                        ],
                        className='six columns',
                        style={'margin-top': '10'}
                    ),
                    html.Div(
                        [
                            html.Div(id="id-histo")
                        ],
                        className='six columns',
                        style={'margin-top': '10'}
                    ),
                ], className="row"
            )
        ],
        className="container"
    )
])


######################################### UPDATING FIGURES #########################################
@app.callback(
    Output('output-container', 'children'),
    [Input('d_virus-name', 'value')])
def _update_legend_gene(virus_name):
    return 'You have selected "{}" virus'.format(virus_name)


@app.callback(
    Output('controls-container_mumps', 'style'),
    [Input('d_virus-name', 'value')])
def _update_mumps_option(virus_name):
    if virus_name == "Mumps":
        return {'display': 'block'}
    else:
        return {'display': 'none'}


@app.callback(
    Output('controls-container_dengue', 'style'),
    [Input('d_virus-name', 'value')])
def _update_dengue_option(virus_name):
    if virus_name == "Dengue":
        return {'display': 'block'}
    else:
        return {'display': 'none'}


@app.callback(
    Output('controls-container_lassa', 'style'),
    [Input('d_virus-name', 'value')])
def _update_lassa_option(virus_name):
    if virus_name == "Lassa":
        return {'display': 'block'}
    else:
        return {'display': 'none'}


@app.callback(
    Output('controls-container_avian', 'style'),
    [Input('d_virus-name', 'value')])
def _update_avian_option(virus_name):
    if virus_name == "Avian":
        return {'display': 'block'}
    else:
        return {'display': 'none'}


@app.callback(
    Output('controls-container_flu', 'style'),
    [Input('d_virus-name', 'value')])
def _update_flu_option(virus_name):
    if virus_name == "Flu":
        return {'display': 'block'}
    else:
        return {'display': 'none'}


@app.callback(
    Output('right-top-graph', 'children'),
    [Input('d_virus-name', 'value'),
     Input('d_mumps', 'value'),
     Input('d_dengue', 'value'),
     Input('d_lassa', 'value'),
     Input('d_avian_opt1', 'value'), Input('d_avian_opt2', 'value'),
     Input('d_flu_opt1', 'value'), Input('d_flu_opt2', 'value'),
     Input('d_flu_opt3', 'value')])
def _update_pĥylogentic_tree(virus_name, mumps, dengue, lassa, avian_opt1, avian_opt2, flu_opt1, flu_opt2, flu_opt3):
    virus_name = virus_name.lower()
    ord_by_elt = "Country"
    data_virus_info = {}
    if virus_name == "ebola" or virus_name == "zika" or virus_name == "measles":
        data_virus_info = {"virus_name": virus_name, "level1": "", "level2": "", "level3": ""}
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name, level1="",
                                                                                                 level2="", level3="")
    elif virus_name == "mumps":
        data_virus_info = {"virus_name": virus_name, "level1": mumps, "level2": "", "level3": ""}
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=mumps,
                                                                                                 level2="", level3="")
    elif virus_name == "dengue":
        data_virus_info = {"virus_name": virus_name, "level1": dengue, "level2": "", "level3": ""}
        tree_file_filtred, metadata_file_filtred, metadata_file_stat = create_paths_file(virus_name, level1=dengue,
                                                                                         level2="", level3="")
    elif virus_name == "lassa":
        data_virus_info = {"virus_name": virus_name, "level1": lassa, "level2": "", "level3": ""}
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=lassa,
                                                                                                 level2="", level3="")
    elif virus_name == "avian":
        data_virus_info = {"virus_name": virus_name, "level1": avian_opt1, "level2": avian_opt2, "level3": ""}
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=avian_opt1,
                                                                                                 level2=avian_opt2,
                                                                                                 level3="")
    elif virus_name == "flu":
        data_virus_info = {"virus_name": virus_name, "level1": flu_opt1, "level2": flu_opt2, "level3": flu_opt3}
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=flu_opt1,
                                                                                                 level2=flu_opt2,
                                                                                                 level3=flu_opt3)

    # tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(data_virus_info)
    if tree_file_filtred in tree_fig:
        fig = tree_fig[tree_file_filtred]
    else:
        if ord_by_elt == "Country" or ord_by_elt == "Division" or ord_by_elt == "Date":
            fig = create_tree(virus_name, tree_file_filtred, metadata_file_filtred, ord_by_elt)

        tree_fig[tree_file_filtred] = fig
    return dcc.Graph(
        id='top-graph',
        figure=fig
    ),


@app.callback(
    Output('graph_map', 'figure'),
    [Input('d_virus-name', 'value'),
     Input('d_mumps', 'value'),
     Input('d_dengue', 'value'),
     Input('d_lassa', 'value'),
     Input('d_avian_opt1', 'value'), Input('d_avian_opt2', 'value'),
     Input('d_flu_opt1', 'value'), Input('d_flu_opt2', 'value'),
     Input('d_flu_opt3', 'value'),
     Input('id-year', 'value')])
def _update_map(virus_name, mumps, dengue, lassa, avian_opt1, avian_opt2, flu_opt1, flu_opt2, flu_opt3, id_year):
    virus_name = virus_name.lower()
    if virus_name == "ebola" or virus_name == "zika" or virus_name == "measles":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name, level1="",
                                                                                                 level2="", level3="")
    elif virus_name == "mumps":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=mumps,
                                                                                                 level2="", level3="")
    elif virus_name == "dengue":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=dengue,
                                                                                                 level2="", level3="")
    elif virus_name == "lassa":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=lassa,
                                                                                                 level2="", level3="")
    elif virus_name == "avian":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=avian_opt1,
                                                                                                 level2=avian_opt2,
                                                                                                 level3="")
    elif virus_name == "flu":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=flu_opt1,
                                                                                                 level2=flu_opt2,
                                                                                                 level3=flu_opt3)
    df = pd.read_csv(metadata_file_stat_filtred)

    min_date, max_date = id_year
    # To select only the data between min_date and max_date
    df = df[df["Year"] >= min_date]
    df = df[df["Year"] <= max_date]
    return create_map_bubble_year(virus_name, metadata_file_stat_filtred, 2, min_date, max_date)


@app.callback(
    Output('id-slicer', 'children'),
    [Input('d_virus-name', 'value'),
     Input('d_mumps', 'value'),
     Input('d_dengue', 'value'),
     Input('d_lassa', 'value'),
     Input('d_avian_opt1', 'value'), Input('d_avian_opt2', 'value'),
     Input('d_flu_opt1', 'value'), Input('d_flu_opt2', 'value'),
     Input('d_flu_opt3', 'value')])
def _update_slicer(virus_name, mumps, dengue, lassa, avian_opt1, avian_opt2, flu_opt1, flu_opt2, flu_opt3):
    virus_name = virus_name.lower()
    if virus_name == "ebola" or virus_name == "zika" or virus_name == "measles":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name, level1="",
                                                                                                 level2="", level3="")
    elif virus_name == "mumps":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=mumps,
                                                                                                 level2="", level3="")
    elif virus_name == "dengue":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=dengue,
                                                                                                 level2="", level3="")
    elif virus_name == "lassa":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=lassa,
                                                                                                 level2="", level3="")
    elif virus_name == "avian":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=avian_opt1,
                                                                                                 level2=avian_opt2,
                                                                                                 level3="")
    elif virus_name == "flu":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=flu_opt1,
                                                                                                 level2=flu_opt2,
                                                                                                 level3=flu_opt3)
    df = pd.read_csv(metadata_file_stat_filtred)
    min_date, max_date = min_max_date(df)
    # create the dictionary of slider
    marks_data = slicer(min_date, max_date)
    min_max_date_value = [min_date, max_date]

    # To select only the data between min_date and max_date
    df = df[df["Year"] >= min_date]
    df = df[df["Year"] <= max_date]
    return dcc.RangeSlider(
                id='id-year',
                min=min_date,
                max=max_date,
                step=1,
                marks=marks_data,
                value=min_max_date_value
            )


@app.callback(
    Output('curve-line-graph', 'figure'),
    [Input('d_virus-name', 'value'),
     Input('d_mumps', 'value'),
     Input('d_dengue', 'value'),
     Input('d_lassa', 'value'),
     Input('d_avian_opt1', 'value'), Input('d_avian_opt2', 'value'),
     Input('d_flu_opt1', 'value'), Input('d_flu_opt2', 'value'),
     Input('d_flu_opt3', 'value'),
     Input('id-year', 'value')])
def _update_curve(virus_name, mumps, dengue, lassa, avian_opt1, avian_opt2, flu_opt1, flu_opt2, flu_opt3, id_year):
    virus_name = virus_name.lower()
    if virus_name == "ebola" or virus_name == "zika" or virus_name == "measles":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name, level1="",
                                                                                                 level2="", level3="")
    elif virus_name == "mumps":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=mumps,
                                                                                                 level2="", level3="")
    elif virus_name == "dengue":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=dengue,
                                                                                                 level2="", level3="")
    elif virus_name == "lassa":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=lassa,
                                                                                                 level2="", level3="")
    elif virus_name == "avian":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=avian_opt1,
                                                                                                 level2=avian_opt2,
                                                                                                 level3="")
    elif virus_name == "flu":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=flu_opt1,
                                                                                                 level2=flu_opt2,
                                                                                                 level3=flu_opt3)
    df = pd.read_csv(metadata_file_stat_filtred)
    min_date, max_date = id_year

    # To select only the data between min_date and max_date
    df = df[df["Year"] >= min_date]
    df = df[df["Year"] <= max_date]

    return create_curve_line(df, virus_name, min_date, max_date)


@app.callback(
    Output('id-histo', 'children'),
    [Input('d_virus-name', 'value'),
     Input('d_mumps', 'value'),
     Input('d_dengue', 'value'),
     Input('d_lassa', 'value'),
     Input('d_avian_opt1', 'value'), Input('d_avian_opt2', 'value'),
     Input('d_flu_opt1', 'value'), Input('d_flu_opt2', 'value'),
     Input('d_flu_opt3', 'value'),
     Input('id-year', 'value')])
def _update_histo(virus_name, mumps, dengue, lassa, avian_opt1, avian_opt2, flu_opt1, flu_opt2, flu_opt3, id_year):
    virus_name = virus_name.lower()
    if virus_name == "ebola" or virus_name == "zika" or virus_name == "measles":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name, level1="",
                                                                                                 level2="", level3="")
    elif virus_name == "mumps":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=mumps,
                                                                                                 level2="", level3="")
    elif virus_name == "dengue":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=dengue,
                                                                                                 level2="", level3="")
    elif virus_name == "lassa":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=lassa,
                                                                                                 level2="", level3="")
    elif virus_name == "avian":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=avian_opt1,
                                                                                                 level2=avian_opt2,
                                                                                                 level3="")
    elif virus_name == "flu":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name,
                                                                                                 level1=flu_opt1,
                                                                                                 level2=flu_opt2,
                                                                                                 level3=flu_opt3)
    df = pd.read_csv(metadata_file_stat_filtred)
    min_date, max_date = id_year

    # To select only the data between min_date and max_date
    df = df[df["Year"] >= min_date]
    df = df[df["Year"] <= max_date]

    # Count the number of viruses by Country
    df_group_by_country = df.groupby(['Country'])['Value'].sum()
    # Translate groupby object in dataframe
    df_group_by_country = df_group_by_country.to_frame()
    # Rename the first column in Value
    df_group_by_country.columns = ['Value']
    # Move the index values (i.e. Country column) in column and reset index
    df_group_by_country = df_group_by_country.reset_index()

    return dcc.Graph(
        id='right-bottom-histo',
        figure={
            'data': [{
                'x': df_group_by_country['Country'],
                'y': df_group_by_country['Value'],
                'type': 'bar'
            }],
            'layout': {
                'paper_bgcolor': 'rgba(0,0,0,0)',
                'plot_bgcolor': 'rgba(0,0,0,0)',
                'autosize': True,
                'margin': '0px 0px 0px 0px',
                'title': '<br>Distribution of {} <br>Between {} and {}'.format(virus_name.title(), min_date, max_date)
            }
        }
    )


######################################### CSS #########################################
external_css = [
    "https://cdnjs.cloudflare.com/ajax/libs/normalize/7.0.0/normalize.min.css",  # Normalize the CSS
    "https://fonts.googleapis.com/css?family=Open+Sans|Roboto"  # Fonts
    "https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css",
    "https://cdn.rawgit.com/TahiriNadia/styles/faf8c1c3/stylesheet.css",
    "https://cdn.rawgit.com/TahiriNadia/styles/b1026938/custum-styles_phyloapp.css"
]

for css in external_css:
    app.css.append_css({"external_url": css})



# Running the server
if __name__ == '__main__':
    app.run_server(debug=True, port=8051)

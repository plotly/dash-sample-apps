library(plotly)

library(dash)
library(dashHtmlComponents)
library(dashCoreComponents)
library(dashTable)
library(rlist)
library("RColorBrewer")

source("datasource.R")
source("word_cloud.R")
source("lda.R")


dbcCard <- function(title, widget) {
  return(
    htmlDiv(
      className = "card",
      children = list(
        htmlDiv(
          className = "card-header",
          children = list(htmlH5(title))
        ),
        htmlDiv(
          className = "card-body",
          children = list(widget)
        )
      )
    )
  )
}


top_10_widget <- function(top_10) {
  return(dccGraph(figure = list(
    data = list(
      list(
        x = top_10$company,
        y = top_10$n,
        type = 'bar',
        #name='SF',
        text = top_10$company,

        #marker = list(
        #  color = 'rgb(158,202,225)',
        #  line = list(color = 'rgb(8,48,107)',
        #              width = 3)
        #),
        text = ~company,
        textposition = 'auto',
        insidetextfont = list(size = 30, color = 'white')
      )
    ),
    layout = list(
      title = "Top 10 banks by number of complaints",
      xaxis = list(
        title = "",
        showlegend = FALSE,
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
      ),
      yaxis = list(title = "")
    )
    ))
  )
}


HEADER <- htmlDiv(
  className = "banner",
  id = "header",
  style = list(
    width = "100%",
    background = "#262B3D",
    color = "#E2EFFA"
  ),
  children = list(htmlDiv(
    className = "container scalable",
    children = list(
      htmlA(htmlImg(
        src = "assets/dash-logo.png",
        style = list(
          "height" = "60px",
          "width" = "auto",
          "margin-bottom" = "30",
          'top' = '40'
        )
      ),
            href = "https://plot.ly/products/dash/"),
      htmlH2(
        htmlA(
          "Top 10 banks by number of complaints",
          className = "title-left",
          style = list(textDecoration = "none",
                       color = "inherit")
        )
      )

    )
  ))
)


LEFT_COLUMN <-
  htmlDiv(
    className = "jumbotron",
    children = list(
      htmlH4(children = "Select bank & dataset size", className = "display-6"),
      htmlHr(className = "my-2"),
      htmlLabel("Select percentage of dataset", className = "lead"),
      htmlP(
        "(Lower is faster. Higher is more precise)",
        style = list("fontSize" = 10, "font-weight" = "lighter")
      ),
      dccSlider(
        id = "n-selection-slider",
        min = 0,
        max = 100,
        step = 10,
        marks = as.list(
          setNames(
            as.character(c("0%", "", "20%", "", "40%", "", "60%", "", "80%", "", "100%")),
            as.character(seq(0, 100, by = 10))
          )
        ),
        value = 20
      ),
      htmlHr(className = "my-3"),
      htmlDiv(list(
        htmlLabel(
          "Select a bank",
          style = list("marginTop" = 59),
          className = "lead"
        )
      )),
      htmlP(
        "(You can use the dropdown or click the barchart on the right)",
        style = list("fontSize" = 10, "font-weight" = "lighter")
      ),
      dccDropdown(
        id = 'bank-selection',
        options = list(),
        value = top_N(get_sample_data(10), 10)$company[1]
      ),

      htmlLabel("Select time frame", className = "lead"),
      dccRangeSlider(
        id = "time-slider",
        min = min(DATASOURCE$date),
        max = max(DATASOURCE$date),
        step = (max(DATASOURCE$date) - min(DATASOURCE$date)) / 10,

        marks = as.list(
          setNames(
            as.character(c("2015 Q1", "Q2", "Q3", "Q4", "2016 Q1", "Q2", "Q3", "Q4", "2017", "Q2", "Q3", "Q4")),
            as.character(c(1420066800, 1427839200, 1435701600, 1443650400, 1451602800, 1459461600, 1467324000, 1475272800, 1483225200, 1490997600, 1498860000, 1506808800))
          )
        ),
        value = list(1420066800, 1506808800)
      )
    )
  )

SELECTION <- htmlDiv(
  className = 'row',
  children = list(
    htmlDiv(children = LEFT_COLUMN, className = "col-md-5"),
    htmlDiv(
      className = "col-md-7",
      dbcCard(
      "Top 10 banks by number of complaints",
        dccLoading(id = 'top_10_widget_container', style = list(width="100%"))
      )
    )
  )
)

TEST <- htmlDiv(children = list(htmlH1(id = 'test-id')))

LDA_PLACEHOLDER <- htmlDiv(
  className = 'row',
  children = list(
    htmlDiv(
      className = "col-md-12",
      children = list(
        dccLoading(
          children = list(
            dccGraph(
              id = "lda"
            )
          )
        )
      )
    )
  )
)

#### TAB
get_tabs <- function(fig_treemap, fig_wordcloud) {
  if (is.na(fig_treemap)) {
    fig_treemap <- dccLoading(
      id = "loading-treemap",
      children = list(dccGraph(id = "bank-treemap")),
      type = "default"
    )
  }
  if (is.na(fig_wordcloud)) {
    fig_wordcloud <- dccLoading(
      id = "loading-wordcloud",
      children = list(dccGraph(id = "bank-wordcloud")),
      type = "default"
    )
  }
  return(
    htmlDiv(
      id = 'tabs',
      className = "row",
      children = list(
        htmlDiv(
          id = 'wordfreq',
          children = dbcCard("Word frequency", fig_treemap),
          className = "col-md-6"
        ),
        htmlDiv(
          id = 'Wordcloud',
          children = dbcCard("Wordcloud", fig_wordcloud),
          className = "col-md-6"
        )
      )
    )
  )
}


get_words_histogram <- function (sample_pct, company, after, before) {
    filtered_df <- filtered_data(sample_pct, company, after, before)
    stoplist <- stoplist_funct(company)
    return(word_cloud(filtered_df$complaint, stoplist))
}


# TODO - move to a separate file
words_histogram <- function(df) {
  # Returns a graph for words frequency

  df <- df[1:25,]

  f <- dccGraph(figure = list(
    data = list(
      list(
        y = df$word,
        x = df$freq,
        type = 'bar',
        #name='SF',
        text = df$word,
        orientation= "h",

        #marker = list(
        #  color = 'rgb(158,202,225)',
        #  line = list(color = 'rgb(8,48,107)',
        #              width = 3)
        #),
        text = ~company,
        textposition = 'auto',
        insidetextfont = list(size = 30, color = 'white')
      )
    ),
    layout = list(
      #title = "Word frequency",
      xaxis = list(
        title = "",
        showlegend = FALSE,
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
      ),
      yaxis = list(
        title = "",
        autorange = "reversed"
      )
    )
  ))

  return(f)
}


words_treemap <- function(df) {
    # Returns a graph for words frequency

    df <- df[1:10,]

    f <- dccGraph(figure = list(
        data = list(
            list(
                # labels=word_list_top, parents=[""] * len(word_list_top), values=freq_list_top
                  #labels = df$word,
                 # text = df$word,
                # values = rep(25, nrow(df)), #df$freq * 100,
                labels  = list("AAA", "BBB"),
                parents = list("", ""),
                type = 'treemap',

                textinfo = "label+value+percent parent+percent entry",
                outsidetextfont = list(size= 20, color= "darkblue"),
                marker = list(line= list(width= 2))
                #,
                #pathbar = list(visible= False}

              )
          )
        #,
          # layout = list(
          #   title = "Word frequency",
          #   xaxis = list(
          #     title = "",
          #     showlegend = FALSE,
          #     zeroline = FALSE,
          #     showline = FALSE,
          #     showticklabels = FALSE,
          #     showgrid = FALSE
          #   ),
          #   yaxis = list(title = "")
          # )
        ))

      return(f)
}


draw_wordcloud <- function(df) {
  # Returns a graph for words frequency

  df <- df[1:25,]

  # append columns with coordinates
  df$size <- 50 * df$freq / max(df$freq)
  df$size <- as.integer(df$size)
  df$i <- seq(1:nrow(df))
  # TODO - take into account size and word length when calculating X?
  df$x <- sapply(df$i, function(i) { return(sqrt(i) * cos(8*pi * sqrt(i) / sqrt(nrow(df)))) } )
  df$y <- sapply(df$i, function(i) { return(sqrt(i) * sin(8*pi * sqrt(i) / sqrt(nrow(df)))) } )


  f <- dccGraph(figure = list(
    data = list(
      list(
        y = df$x,
        x = df$y,
        type = 'scatter',
        mode = 'markers+text',
        text = df$word,

        marker = list(
          color = 'rgb(255,255,255)',
          size = rep(0, nrow(df))
        ),

        textfont = list(
          size = df$size
        )
      )
    ),

      layout = list(
        #title = "Word cloud",
        xaxis = list(
          title = "",
          showlegend = FALSE,
          zeroline = FALSE,
          showline = FALSE,
          showticklabels = FALSE,
          showgrid = FALSE
        ),
        yaxis = list(
          title = "",
          showlegend = FALSE,
          zeroline = FALSE,
          showline = FALSE,
          showticklabels = FALSE,
          showgrid = FALSE
        )
      )
    )
  )

  return(f)
}

app <- Dash$new(
  external_stylesheets=list("https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css")
)

app$layout(
  htmlDiv(
    className = "container",
    children = list(
      HEADER,
      SELECTION,
      TEST,
      LDA_PLACEHOLDER
    )
  )
)

## -- callbacks --

app$callback(
  output('top_10_widget_container', 'children'),
  params = list(
    input('n-selection-slider', 'value')
  ),
  function(n_selected) {
    sample <- get_sample_data(n_selected)
    top10 <- top_N(sample, 10)
    return(
      top_10_widget(top10)
    )
  }
)

app$callback(
  output('bank-selection', 'options'),
  params = list(
    input('n-selection-slider', 'value')
  ),
  function(n_selected) {
    sample <- get_sample_data(n_selected)
    top10 <- top_N(sample, 10)

    available_indicators <- unique(top10$company)

    option_indicator <- lapply(available_indicators,
                           function(available_indicator) {
                             list(label = available_indicator,
                                  value = available_indicator)
                           })


    return(
      option_indicator
    )
  }
)


app$callback(
  output('test-id', 'children'),
  params = list(
    input('n-selection-slider', 'value'),
    input('bank-selection', 'value'),
    input('time-slider', 'value')
  ),
  function(sample_pct, company, selected_time) {
    time <- selected_time
    after <- as.integer(selected_time[1])
    before <- as.integer(selected_time[2])

    hist_data <- get_words_histogram(as.integer(sample_pct), company, after, before)
    return(
        get_tabs(
            words_histogram(hist_data),
            draw_wordcloud(hist_data)
        )
    )
  }
)



app$callback(
  output('lda', 'figure'),
  params = list(
    input('n-selection-slider', 'value'),
    input('bank-selection', 'value'),
    input('time-slider', 'value')
  ),
  function(sample_pct, company, selected_time) {
    time <- selected_time
    after <- as.integer(selected_time[1])
    before <- as.integer(selected_time[2])

    df <- filtered_data(ceiling(sample_pct / 10), company, after, before)
    # df <- df[1:25,]

    lda_df <- build_lda_df(df)
    
    num_topics <- max(tsne_df$topic)
    tops <- top_terms(tsne_df)
    top_terms_by_topic <- aggregate(term~topic, tops, paste, collapse = ',')
    colors <- brewer.pal(n = num_topics, name = "RdBu")
    
    graphs <- list()
    for (i in seq(1, num_topics)) {
      df_i <- lda_df[lda_df$topic == i, ]
      g <- list(
        y = df_i$V1,
        x = df_i$V2,
        type = 'scatter',
        mode = 'markers',
        #text = df$word,
        name = top_terms_by_topic[top_terms_by_topic$topic == i, ]$term,
        marker = list(
          color = colors[i]
          # size = rep(0, nrow(lda_df))
        )
      )
      
      graphs <- list.append(graphs, g)
    }
    
    lda_scatter_figure <- list(
      data = graphs,
      
      layout = list(
        title = "LDA",
        xaxis = list(
          title = "",
          showlegend = FALSE,
          zeroline = FALSE,
          showline = FALSE,
          showticklabels = FALSE,
          showgrid = TRUE
        ),
        yaxis = list(
          title = "",
          showlegend = FALSE,
          zeroline = FALSE,
          showline = FALSE,
          showticklabels = FALSE,
          showgrid = TRUE
        )
      )
    )


    return(lda_scatter_figure)
  }
)


app$run_server(showcase = FALSE)

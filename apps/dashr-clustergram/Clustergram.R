library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(dashBio)
library(heatmaply)
library(data.table)


app <- Dash$new()

options_tabs <- htmlDiv(id = 'clustergram-body', className = 'app-body', children = list(
  dccLoading(className = 'dashbio-loading', children = htmlDiv(
    children = dccGraph(id = 'clustergram', style = list('display' = 'none'))
  )),
  
  htmlDiv(id = 'clustergram-control-tabs', className = 'control-tabs', children = list(
    dccTabs(id = 'clustergram-tabs', value = 'what-is', children = list(
      dccTab(
        label = 'About',
        value = 'what-is',
        children = htmlDiv(className = 'control-tab', children = list(
          htmlH4(className = 'what-is', children = "What is Clustergram?"),
          htmlP('Clusterman is a combination of a heatmap and dendrograms that
                allows you to display hierarchical clustering data. Clusters on
                the dendrogram are highlighted in one color if they comprise
                 data points that share some minimal level of correlation.'),
          htmlP('In the "Data" tab, you can select a preloaded dataset or, 
                alternatively, upload one of your own. A sample dataset is also
                available for download in the tab.'),
          htmlP('In the "Graph" tab, you can choose the dimension(s) along which
                clustering will be performed (row or column). You can also change
                the threshold that determines the point at which clusters are 
                highlighted for the tow and column dendrograms, and choose which rows
                and columns are used to compute the clustering.'),
          htmlP('In addition, you can highlight specific clusters by adding annotations
                to the clustergram, and choose whether to show or hide the labels for
                the rows and/or columns.')
        ))
      ),
      dccTab(
        label = 'Data',
        value = 'datasets',
        children = htmlDiv(className = 'control-tab', children = list(
          
          htmlDiv(
            'Preloaded dataset',
            title = 'Choose from some pre-loaded datasets to
            view them on the heatmap.',
            className = 'fullwidth-app-controls-name'
          ),
          
          dccDropdown(
            id = 'clustergram-datasets',
            options = list(
              list('label' = 'Andersons Iris Data', 'value' = 'iris'),
              list('label' = 'mtcars', 'value' = 'mtcars'),
              list('label' = 'Prostate Cancer', 'value' = 'prostatecancer'),
              list('label' = 'Lung Cancer Subtypes', 'value' = 'lungcancer')
            ),
            value = 'prostatecancer'
          ),
          
          a
        ))
      )
    ))
  ))
))













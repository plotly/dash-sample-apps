library(dplyr)
library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(data.table)
library(plotly)
library(rlist)
library(glue)
library(parallel)
library(compiler)
library(dtplyr)
library(knitr)
library(png)
library(magrittr)
library(Rtsne)
library(stringr)
library(jsonlite)
library(tictoc)

appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}

#Following the python notation
IMAGE_DATASETS <- c("mnist_3000", "cifar_gray_3000", "fashion_3000")
WORD_EMBEDDINGS <- c("wikipedia_3000", "twitter_3000", "crawler_3000")


readme <- htmlDiv(children=list(
  htmlButton("Close", id='close-button', n_clicks=0)
  ,
  dccMarkdown("# t-SNE Explorer

This is a demo of the Dash for R framework developed by Plotly [Plotly](https://plot.ly/).

Dash abstracts away all of the technologies and protocols required to build an interactive web-based application and is a simple and effective way to bind a user interface around your R code. To learn more check out our [documentation](https://plot.ly/dash). 

For an introductory and extensive explanation of t-SNE how to use it properly, please check out the [demo app](https://dash-tsne.plot.ly/).
")
  ,
  htmlImg(id='gif1', src='assets/animated1.gif')
  ,
  htmlImg(id='gif2', src='assets/animated2.gif')
  ,
  dccMarkdown("
## Getting Started
### Using the demo
To get started, choose a dataset you want to visualize. When the scatter plot appears on the graph, you can see the original image by clicking on a data point. 

Alternatively, you can explore the GloVe Word Vectors datasets, which are encoded vectors of large collection of texts from Wikipedia, Twitter, and acquired through Web Crawlers. Upon clicking a data point, you will be able to see the 5 closest neighbors of the word you clicked.

### Running the app locally



Clone the git repo, then install the requirements 
```
git clone https://github.com/plotly/dash-sample-apps/dashr-tsne.git
cd dash-sample-apps/dashr-tsne

```
Download the script and run init.R

```
Rscript init.R
```
Run the app
```
Rscript app.R
```
or using ctrl+shift+enter (cmd+shift+enter)
### How to use the local version
To train your own t-SNE algorithm, input a high-dimensional dataset with only numerical values, and the corresponding labels inside the upload fields. For convenience, small sample datasets are included inside the data directory. The training can take a lot of time depending on the size of the dataset (the complete MNIST dataset could take 15-30 min), and it is not advised to refresh the webpage when you are doing so.

### Generating data
Will be included

which will create the csv file with the corresponding parameters. At the moment, we have the following datasets:
* MNIST
* CIFAR10
* Fashion_MNIST


## About the app
### What is t-SNE?
t-distributed stochastic neighbor embedding, created by van der Maaten and Hinton in 2008, is a visualization algorithm that reduce a high-dimensional space (e.g. an image or a word embedding) into two or three dimensions, facilitating visualization of the data distribution. 

A classical example is MNIST, a dataset of 60,000 handwritten digits, 28x28 grayscale. Upon reducing the set of images using t-SNE, you can see all the digit clustered together, with few outliers caused by poor calligraphy. [You can read a detailed explanation of the algorithm on van der Maaten's personal blog.](https://lvdmaaten.github.io/tsne/)


## Built With
* [DashR](https://dash.plot.ly/) - Main server and interactive components
* [Plotly R](https://plot.ly/r/) - Used to create the interactive plots
* [data.table](https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html)
* [Rtsne](https://cran.r-project.org/web/packages/tsne/tsne.pdf) - Run the t-SNE algorithm
* [dplyr](https://www.tidyverse.org/)

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

* **Xing Han Lu** - *Python version* - [@xhlulu](https://github.com/xhlulu)
* **Dan Yunheum Seol** -  *Initial work* - [@dan-seol](https://github.com/dan-seol)

See also the list of [contributors](https://github.com/plotly/dash-svm/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

##Screenshots
")
  ,
  htmlImg(id='png1', src='assets/demo_image1.png')
  ,
  htmlImg(id='png2', src='assets/fashion_mnist_example.png')
  ,
  htmlImg(id='png3', src='assets/screenshot1.png')
  ,
  htmlImg(id='png4', src='assets/screenshot2.png')
))

genMarks <- function(n){
  l <- list(glue('{n}'))
  names(l) <- 'label'
  m <- list(l)
  names(m) <- glue('{n}')
  return(m)
}

input_field <- function(title, state_id, state_value, state_max, state_min){
  #takes (title, state_id, state_value, state_max, state_min), returns htmlDiv(htmlP(title), dccInput(...))
  return(htmlDiv(list(
    htmlP(title), 
    dccInput(id=state_id, type='number', value=state_value, max=state_max, min=state_min, size='7')
    
  ))
  
  )}

genMark <- function(n){
  l <- list(glue('{n}'))
  names(l) <- 'label'
  return(l)
}

genMarks <- function(min, max, by){
  s <- seq(from=min, to=max, by)
  l <- lapply(s, genMark)
  names(l) <- s
  return(l)
}

namedSlider <- function(name, short, min, max,  val, marks=NULL){
  if(is.null(marks)){
    marks = genMarks(min, max, step)
  } else {
    step=NULL
    
  }
  
  return(htmlDiv(
    
    style=list(margin='25px 5px 30px 0px')
    ,
    children = list(
      
      glue('{name}:')
      ,
      htmlDiv(
        style=list('margin-left'='5px'),
        children = dccSlider(id=glue('slider-{short}'), min=min, max=max, marks=marks, value=val)
      )
    )))
}

namedRadioItems <- function(name, short, options, val, ...){
  kwargs <- list(...)
  kwargs$display <- 'inline-block'
  return(htmlDiv(
    id = glue('div-{short}')
    , 
    style = kwargs
    , 
    children=list(glue({name}), 
                  dccRadioItems(
                    id= glue('radio-{short}'), 
                    options=options, 
                    value=val, 
                    labelStyle=list(
                      display = 'inline-block', 
                      'font-weight'=300
                    ),
                    style = list(
                      display='inline-block'
                      
                    )
                  ))
  ))
}

Card <- function(kids, ...){
  kwargs <- list(...)
  kwargs$padding<-20
  kwargs$margin<-5
  kwargs$borderRadius<-5
  kwargs$border<-'thin lightgrey solid'
  kwargs[['user-select']] <-  'none'
  kwargs[['-moz-user-select']] <-  'none'
  kwargs[['-webkit-user-select']] <-  'none'
  kwargs[['-ms-user-select']] <- 'none'
  return(htmlSection(children = kids, style= kwargs))
  
}

#Generate the default 3D scatter plot
TSNE_DT <- fread('./data/tsne_3d.csv')

colnames(TSNE_DT) <- c("Digit", "x", "y", "z")

defaultPlot <- plot_ly(
  TSNE_DT, 
  type='scatter3d', 
  x=~x, 
  y=~y, 
  z=~z, 
  color = ~as.factor(Digit), 
  mode='markers', 
  marker=list(symbol='circle', size=2.5)
) %>% plotly::layout(.,width=850, height=600)


#Initialize the app

app <- Dash$new()
##############################################################################################################################

###########LAYOUT BEGINS########################################################################################################


app$layout(
  htmlDiv( className='row, background', children=list(
    htmlDiv(id='footer', children=readme)
    ,
    htmlDiv(className='row', children=list(
      
      htmlH2('t-SNE Explorer', className='title', id='app-title')
      ,
      htmlImg(id='plotly-image', src="assets/logo.png", style=list(height='90px', float='right', 'margin-top'='10px', 'margin-right'='30px'))
      
    )
    ,
    style=list('background-color'='rgb(255,255,255)')
    )
    ,
    htmlDiv(className='row, background', list(
      
      htmlDiv(className='row, background', 
              
              children=list(dccMarkdown(
                "The Scatter plot above is the result of running the t-SNE algorithm on the MNIST digits, resulting in a 3D visualization of the image dataset.   
For demo purposes, all the data were pre-generated using limited number of input parameters, a subset of 3000 samples, and displayed instantly.   
You can run the t-SNE algorithm with your custom dataset and label using a dropdown at the top right. To learn more about how the t-SNE Explorer works, click on 'Learn More' below.  
This is the R version of the t-SNE explorer. To view the source code, please visit the [GitHub Repository](https://github.com/plotly/dashr-tsne)"),
                htmlButton('Learn More', id='learn-button', n_clicks=0)
                
              )
              ,
              style=list('margin-left'='5px')
      )
      ,
      dccDropdown(
        id = 'dropdown-mode-choice',
        options=list(
          list(label = 'Pre-generated data', value='demo'),
          list(label = 'Custom data', value = 'custom')
          
        ),
        value = 'demo',
        clearable = FALSE,
        style = list(width='45%', float='right')
      )
      ,
      
      htmlDiv(className='row, background', list(   
        
        htmlH4( children='', id='tsne_h4', style=list(float='left')),
        htmlDiv(className= 'three columns', id='control-panel', children = list(), style=list(width='25vh'))
        
      ),  style=list(float='left', display='block', width='27vw'))
      ,
      
      htmlDiv(
        id='demo-graph', 
        className='five columns', 
        children= dccGraph(id='tsne-3d-plot', figure=defaultPlot),
        
      )
      ,
      htmlDiv(
        id='custom-graph',
        className='five columns', 
        children= dccGraph(id='tsne-3d-plot-custom', figure=defaultPlot ), 
        style=list(display='none'))
      ,
      htmlDiv(list(htmlDiv(id='KL-div',style=list(display='none') )
                   ,
                   htmlDiv(id='end-time', style=list(display='none'))
                   , 
                   htmlDiv(id='error-message', style=list(display='none'))
                   
                   
      ), id='plot-div', style=list(display='none'))
      ,
      htmlDiv(id='palette', className='three columns' )
      ,
      htmlDiv(id='custom-container', className='three columns')
    ))
    
  )
  ,
  style=list('max-width'='100%', 'font-size'='1.5rem', padding='0px 0px') )
)
################################################################LAYOUT DONE################################################

################################################################CALLBACK BEGINS############################################
#

app$callback(
  
  output=list(id='footer', property='style')
  ,
  params=list(
    input(id='close-button', property='n_clicks'), 
    input(id='learn-button', property='n_clicks')
  )
  ,
  function(c, l){
    if(l>c){
      return(list(display='block', backgroundColor = '#F9F9F9'))
    } else{
      return(list(display='none'))
    }
  }
  
)

demoPanel <- list(   
  dccDropdown(
    id = 'dropdown-dataset',
    options=list(
      list(label = 'MNIST Digits', value='mnist_3000'),
      list(label = 'Twitter (Glove)', value = 'twitter_3000'),
      list(label = 'Wikipedia (Glove)', value = 'wikipedia_3000')
      
    )
    ,
    placeholder="Select a dataset"
    ,
    value = ''
    ,
    clearable=F
    ,
    style = list(width='100%')
  )
  ,
  namedSlider("Number of Iterations",'iterations', 250, 1000,  500, genMarks(250, 1000, 250))
  ,
  namedSlider("Perplexity",'perplexity', 3, 100,  30, list('3'=list(label='3'), '10'=list(label='10'), '30'=list(label='30'), '50'=list(label='50'), '100'=list(label='100')))
  ,
  namedSlider("Initial PCA dimensions",'pca-dimensions', 25, 100,  50, list('25'=list(label='25'), '50'=list(label='50'), '100'=list(label='100')))
  ,
  namedSlider("Learning Rate",'learning-rate', 10, 200,  100, list('10'=list(label='10'), '50'=list(label='50'), '100'=list(label='100'), '200'=list(label='200')))
  ,
  htmlDiv(
    id='div-wordemb-controls'
    ,
    style=list(display='none', float='left')
    ,
    children=list(
      
      namedRadioItems(
        name="Display",
        short='wordemb-display-mode',
        options= list(
          list(label=' Regular', value='regular'), 
          list(label=' Top-100 Neighbors', value='neighbors')
        )
        ,
        val='regular'
      )
      ,
      htmlDiv(dccDropdown(
        id='dropdown-word-selected', 
        placeholder='Select a word to display its neighbors', value=''), 
        style=list(width='25vh'))
    )
  )
  
)
###########################################demo PANEL ends
#
#############################################custom PANEL begins


#
###############################################################################################################################
app$callback(
  output = list(id='control-panel', property='children'),
  params = list(input(id='dropdown-mode-choice', property='value')),
  function(choice){
    
    if(choice == 'demo'){
      return(demoPanel)
    } else {
      return(
        list(
          htmlDiv(id='data-df-and-message', children = '', style=list(display='none')),
          htmlDiv(id='label-df-and-message', children = '', style=list(display='none')),
          input_field('Number of Iterations:', "n-iter-state", 400, 1000, 250),
          input_field('Perplexity:','perplexity-state', 20, 50, 5),
          input_field("Learning Rate:", "lr-state", 200, 1000, 10),
          input_field("Initial PCA dimensions", 'pca-state', 30, 10000, 3),
          htmlButton(
            id='tsne-train-button', 
            n_clicks=0, 
            children='Start Training t-SNE'
          )
          ,
          dccUpload(
            id='upload-data'
            , 
            children=htmlA('Upload your input data here.')
            , 
            style=list(
              height= '45px',
              'line-height'= '45px',
              'border-width'= '1px',
              'border-style'= 'dashed',
              'border-radius'= '5px',
              'text-align'= 'center',
              'margin-top'= '5px',
              'margin-bottom'= '5 px'
            )
            , 
            multiple=F
            ,
            max_size = -1
          )
          ,
          dccUpload(
            id='upload-label'
            ,
            children=htmlA('Upload your labels here.')
            ,
            style=list(
              height= '45px',
              'line-height'= '45px',
              'border-width'= '1px',
              'border-style'= 'dashed',
              'border-radius'= '5px',
              'text-align'= 'center',
              'margin-top'= '5px',
              'margin-bottom'= '5 px'
            )
            , 
            multiple=F
            ,
            max_size = -1
          )
          ,
          htmlDiv(list(
            
            htmlP(id='upload-data-message', style = list('margin-bottom'='0px'))
            ,
            htmlP(id='upload-label-message', style = list('margin-bottom'='0px'))
            ,
            htmlDiv(
              id ='training-status-message', 
              style=list('margin-bottom'='0px', 'margin-top'='0px')),
            htmlP(id='error-status-message')
            
          )
          ,
          id='output-messages',
          style=list('margin-bottom'='2px', 'margin-top'='2px')
          )
        ))
    }
  }
  
)
app$callback(
  output = list(id='palette', property='children'),
  params=list(input(id='dropdown-mode-choice', property='value')),
  function(choice)
  {
    if(choice=='demo'){
      return(list(
        
        htmlDiv(
          id='div-plot-click-message'
          ,
          style=list(
            'text-align'='center', 
            'margin-bottom'='7px', 
            'font-weight'='bold'
          )
        ),
        htmlDiv(id='div-plot-click-image')
        ,
        htmlDiv(id='div-plot-click-wordemb', style=list(width='40vh', height='30vh'))
      ))
    } else{
      return("")
    }
  }
  
)
app$callback(
  output = list(id='demo-graph', property='style'),
  params = list(input(id='dropdown-mode-choice', property='value')),
  function(value){
    if(value=='demo'){
      return(NULL)
      
    } else{
      return(list(display='none'))
    }
  }
)
app$callback(
  output = list(id='custom-graph', property='style'),
  params = list(input(id='dropdown-mode-choice', property='value')),
  function(value){
    if(value=='custom'){
      return(NULL)
    } else{
      return(list(display='none'))
    }
  }
)

## Global varlable for the data storage required.
#
#
datLst <<- list(mnist_3000 = fread("data/mnist_3000_input.csv.gz", header=T), 
                fashion_3000= fread("data/fashion_3000_input.csv.gz", header=T),
                cifar_gray_3000 = fread("data/cifar_gray_3000_input.csv.gz", header=T),
                wikipedia_3000 = fread("data/wikipedia_3000.csv.gz", header=T),
                crawler_3000 = fread("data/crawler_3000.csv.gz", header=T),
                twitter_3000 = fread('data/twitter_3000.csv.gz', encoding='Latin-1', header=T)
)

datLabels <<- list(
  mnist_3000 = unique(fread("data/mnist_3000_labels.csv.gz", header=T)), 
  fashion_3000= unique( fread("data/fashion_3000_labels.csv.gz", header=T)),
  cifar_gray_3000 = unique(fread("data/cifar_gray_3000_labels.csv.gz", header=T)),
  wikipedia_3000 = unique(fread("data/wikipedia_3000.csv.gz", header=T)[, 1]),
  crawler_3000 = unique(fread("data/crawler_3000.csv.gz", header=T)[, 1]),
  twitter_3000 =  unique(fread('data/twitter_3000.csv.gz', encoding='Latin-1', header=T)[,1])
)

app$callback(
  output = list(id='div-wordemb-controls', property='style')
  ,
  params = list(input(id='dropdown-dataset', property='value'))
  ,
  function(dataset){
    if(is.element(dataset, WORD_EMBEDDINGS)){
      return(list(width='20vh', display='block', float='left'))
    } else {
      return(list(display='none'))
    }
  }
)
#
app$callback(
  output = list(id='dropdown-word-selected', property='disabled')
  ,
  params = list(input(id='radio-wordemb-display-mode', property='value'))
  ,
  function(mode){
    if (mode == 'neighbors'){
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
  
)
#
#helper function
labeledList <- function(n){
  l <- list()
  l$label <- n
  l$value <-n
  return(l)
}

#
app$callback(
  output = list(id='dropdown-word-selected', property='options')
  ,
  params = list(input(id='dropdown-dataset', property='value'))
  ,
  function(dataset){
    
    D <- datLabels[[dataset]]
    L <- lapply(unlist(D), labeledList)
    
    return(unname(L))
  }
)

#############CALLBACK that shows 3D plot#####


### Euclidean distance
eu.norm <- function(v1, v2){
  return(sqrt(sum((v1-v2)^2)))
}
#
### functions that enable us to load
#
generate_figure_word_vec <- function(embedding_df, wordemb_display_mode, selected_word, datTbl){
  
  colnames(embedding_df) <- c("label", 'x', 'y', 'z')
  
  #Regular displays the full scatter plot with only circles
  if(wordemb_display_mode == 'regular'){
    plot_mode <- 'markers'
    
    p <- plot_ly(
      embedding_df, 
      type='scatter3d', 
      x=~x, 
      y=~y, 
      z=~z, 
      color = '#ED9C69',
      text=~label, 
      mode='markers', 
      marker=list(symbol='circle', size=2.5))
    
    return(p)
  }
  # Nearest Neighbors display only the 100 nearest neighbors of the selected_word, in text rather than circles
  else if(wordemb_display_mode == 'neighbors'){
    if(selected_word == ''){
      return(plot_ly(type =  'scatter3d', x=0, y=0, z=0, mode='text', text='NA'))
    }
    plot_mode <- 'text'
    #Get the nearest neighbors incides using Euclidean distance
    
    selected_vector <- filter(datTbl, label==selected_word)
    mtx <- as.matrix(datTbl[, -1])
    sel_vec <- unlist(selected_vector[, -1])
    datTbl$distance <-  mtx %>% apply(., 1, . %>% eu.norm(., sel_vec))
    sorted_DT <- datTbl[order(distance), ]
    neighbors_label <- sorted_DT[1:101, 1]
    neighbors <- embedding_df[label %in% unlist(neighbors_label)]
    return(plot_ly(
      neighbors, 
      type='scatter3d', 
      x=~x, 
      y=~y, 
      z=~z,
      name=~label,
      text=~label, 
      textposition='middle-center', 
      showlegend=FALSE, 
      mode='text'
    ))
    
  }
}

app$callback(
  output = list(id='tsne-3d-plot', property='figure')
  ,
  params = list(input(id='dropdown-dataset', property='value'),
                input(id='slider-iterations',  property='value'),
                input(id='slider-perplexity',  property='value'),
                input(id='slider-pca-dimensions',  property='value'),
                input(id='slider-learning-rate',  property='value'),
                input(id='dropdown-word-selected',  property='value'),
                input(id='radio-wordemb-display-mode',  property='value')
  )
  ,
  function(dataset, iterations, perplexity, pca_dim, learning_rate, selected_word, wordemb_display_mode){
    if(dataset==''|| is.null(dataset)){
      return(defaultPlot)
    } else {
      axes <- list(title='', showgrid=T, zeroline=F, showticklabels=F)
      path <- as.character(glue('demo_embeddings/{dataset}/iterations_{iterations}/perplexity_{perplexity}/pca_{pca_dim}/learning_rate_{learning_rate}/data.csv'))
      embedding_DT <- fread(path, encoding='Latin-1', header=T)
      datTbl <- datLst[[dataset]]
      colnames(datTbl)[1] <- 'label'
      if(is.element(dataset, IMAGE_DATASETS)){
        
        colnames(embedding_DT) <- c('label', 'x','y','z')
        p <- plot_ly(
          embedding_DT, 
          type='scatter3d', 
          x=~x, 
          y=~y, 
          z=~z, 
          color=~as.factor(label), 
          mode='markers', 
          marker=list(symbol='circle', size=2.5))
        
        return(p %>% plotly::layout(.,width=850, height=600))
      } else if(is.element(dataset, WORD_EMBEDDINGS)){
        figure <- generate_figure_word_vec(embedding_DT, wordemb_display_mode, selected_word, datTbl)
        
        return(figure %>% plotly::layout(.,width=850, height=600))
      }else{
        return(plot_ly(type =  'scatter3d'))}
      
    }
  }
)
#
#### Returning the 5 nearest neighborhood for the point clicked

app$callback(
  output =  list(id = 'div-plot-click-wordemb', property='children')
  ,
  params = list(
    input(id = 'tsne-3d-plot', property='clickData')
    ,
    input(id='dropdown-dataset', property='value')
  )
  ,
  function(clickData, dataset){
    if(dataset %in% WORD_EMBEDDINGS){
      clickDat <- clickData[[1]] 
      # contains a list(character(6))
      selected_word <- clickDat[[1]][6]
      datTbl <- datLst[[dataset]]
      colnames(datTbl)[1] = 'label'
      selected_vector <- filter(datTbl, label==selected_word)
      mtx <- as.matrix(datTbl[, -1])
      sel_vec <- unlist(selected_vector[, -1])
      datTbl$distance <- mtx %>% apply(., 1, . %>% eu.norm(., sel_vec))
      sorted_DT <- datTbl[order(distance), ]
      neighbors_label <- sorted_DT[2:6, 1]
      p <- plot_ly(type='bar', y=neighbors_label$label, orientation='h' ) %>% layout( xaxis= list(title='Euclidean Distance'))
      if(selected_word==''||is.na(selected_word)){toplabel <- "Default 5-NN graph"} 
      else{toplabel <- paste("The 5 nearest neighbors of", selected_word)}
      return(list(htmlH5(toplabel, style=list(float='center')), dccGraph(
        id='graph-bar-nearest-neighbors-word',
        figure=p,
        style=list(height= '25vh'),
        config=list('displayModeBar'= F)) %>% htmlDiv(., 
                                                      style=list(
                                                        height='25vh', 
                                                        display='block', 
                                                        margin='auto'))))
    } else {
      return(list())
    }
  }
  
)

#
####
app$callback(
  output = list(id='div-plot-click-image',property='children')
  ,
  params = list(
    input(id='tsne-3d-plot', property='clickData')
    ,
    input(id='dropdown-dataset', property='value')
    ,
    input(id='slider-iterations', property='value')
    ,
    input(id='slider-perplexity', property='value')
    ,
    input(id='slider-pca-dimensions', property='value')
    ,
    input(id='slider-learning-rate', property='value')
    
  )
  ,
  function(clickData, dataset, iterations, perplexity, pca_dim, learning_rate){
    #clickData = list(list(x=0, y=0, z=0, curveNumber=1, pointNumber=291))
    if(dataset %in% IMAGE_DATASETS){
      datTbl <- datLst[[dataset]]
      path = as.character(glue('demo_embeddings/{dataset}/iterations_{iterations}/perplexity_{perplexity}/pca_{pca_dim}/learning_rate_{learning_rate}/data.csv'))
      
      embedding_DT <- fread(path, encoding='Latin-1', header=T)
      
      clickDat <- clickData[[1]]
      embedding_DT$x = embedding_DT$x %>% round(., 4)
      embedding_DT$y = embedding_DT$y %>% round(., 4)
      embedding_DT$z = embedding_DT$z %>% round(., 4)
      clickPoint = c(clickDat[[1]][1], clickDat[[1]][2], clickDat[[1]][3]) %>% as.numeric(.) %>% round(., 4)
      imageIndex <- which(embedding_DT$x== clickPoint[1] & embedding_DT$y==clickPoint[2] & embedding_DT$z==clickPoint[3])
      
      imageVec<- as.numeric(datTbl[imageIndex, ])
      #
      if(dataset=='cifar_gray_3000'){ imageMtx <- matrix(imageVec, nrow=32)} else {imageMtx <- matrix(imageVec, nrow=28)}
      DIGIT <- writePNG(t(imageMtx))
      DIGIT_B64 <- base64_enc(DIGIT)
      
      return(htmlImg(
        src =  glue('data:image/png;base64, ', DIGIT_B64), 
        style =  list(height='25vh', display='block', margin='auto', float='right')
      ))
      
    } else{
      return(list())
    }
  }
  
)
#########CALLBACKS for DEMO ends here
#
#
#########CALLBACKS for custom data begins
#

#helper function
parse_contents <- function(contents, filename){
  contents_parsed <- stringr::str_split(contents, ',')
  decoded <- jsonlite::base64_dec(contents_parsed[[1]][2]) %>% rawToChar(.)
  if(str_detect(filename, '.csv')){
    DT <- fread(decoded)
    
  } else if(str_detect(filename, '.xls')){
    DT <- read_xls(decoded)
  } else{
    print("Invalid file")
    return(list(NULL, message='INVALID FILE'))
  }
  
  return(list(dt=DT, message=glue('Data {filename} has been successfully uploaded')))
}


# main callback
# Store the uploaded data to the hidden data htmlDiv container
# 
app$callback(
  
  output = list(id='data-df-and-message', property='children')
  ,
  params = list(
    
    input(id='upload-data', property='contents')
    ,
    input(id='upload-data', property='filename')
    
  )
  ,
  function(contents, filename){
    contents_parsed <- parse_contents(contents, filename)
    if(is.na(filename)){return(NULL)}
    if(is.na(contents_parsed)){return(NULL)}
    dt <-contents_parsed$dt %>% toJSON(., force=TRUE)
    message <- contents_parsed$message %>% toJSON(., force=TRUE)
    
    C <-c(dt, message) %>% toJSON(., force=TRUE)
    return(C)
  }
)


app$callback(
  
  output = list(id='upload-data-message', property='children')
  ,
  params = list(
    
    input(id='data-df-and-message', property='children')
  )
  ,
  function(JSON){
    if(is.na(JSON)){return("NO DATA LOADED")}
    
    box <- fromJSON(JSON)
    
    return(box[2] %>% fromJSON)
  }
  
)


app$callback(
  
  output = list(id='label-df-and-message', property='children')
  ,
  params = list(
    
    input(id='upload-label', property='contents')
    ,
    input(id='upload-label', property='filename')
    
  )
  ,
  function(contents, filename){
    contents_parsed <- parse_contents(contents, filename)
    if(is.na(filename)){return(NULL)}
    if(is.na(contents_parsed)){return(NULL)}
    dt <-contents_parsed$dt %>% toJSON(., force=TRUE)
    message <- contents_parsed$message %>% toJSON(., force=TRUE)
    
    C <-c(dt, message) %>% toJSON(., force=TRUE)
    return(C)
  }
)


app$callback(
  
  output = list(id='upload-label-message', property='children')
  ,
  params = list(
    input(id='label-df-and-message', property='children')
  )
  ,
  function(JSON){
    if(is.na(JSON)){return("NO DATA LOADED")}
    box <- fromJSON(JSON)
    return(box[2] %>% fromJSON)
  }
  
)
## pre- t-SNE ends

####t-SNE process

# Store the uploaded data to the hidden data htmlDiv container

softBounds <- function(c, a=c, b=c){
  if(c<a){
    return(a)
  } else if(b<c){
    return(b)
  } else{
    return(c)
  }
}
is.empty <- function(x){
  return(is.null(x)||is.na(x))
}

app$callback(
  
  output = list(id='custom-container', property='children')
  ,
  params = list(
    
    input(id='tsne-train-button', property= 'n_clicks')
    ,
    state(id='perplexity-state',property='value' )
    ,
    state(id='n-iter-state',property='value' )
    ,
    state(id='lr-state',property='value' )
    ,
    state(id='pca-state',property='value' )
    ,
    state(id='data-df-and-message',property='children' )
    ,
    state(id='label-df-and-message',property='children' )
  )
  ,
  function(n_clicks, perplexity, n_iter, learning_rate, pca_dim, data_div, label_div){
    
    #RUN THE t-SNE ALGORITHM UPLON CLICKING THE TRAINING BUTTON
    error_message = "No error"
    # Fix up for startup post
    if((n_clicks <= 0)|| is.empty(data_div) || is.empty(label_div) ){
      kl_js <- toJSON(NULL, null='null', force=TRUE)
      time_js <- toJSON('0 seconds elapsed')
      error_js <- toJSON('Empty')
      C <- c(kl_js, time_js, error_js, 'null') %>% toJSON(., raw='base64')
      return(C)
      
    } else{
      Data_DT <- fromJSON(fromJSON(data_div)[1])
      label_DT <- fromJSON(label_div)[1] %>% fromJSON()
      niter <- softBounds(n_iter, 250, 1000)
      perp <- softBounds(perplexity, 5, 50)
      lr <- softBounds(learning_rate, 10, 1000)
      initialdim <- softBounds(pca_dim, 3, dim(Data_DT)[2])
      tic()
      TSNE <- Rtsne(
        Data_DT, 
        dims=3, 
        initial_dims = initialdim, 
        perplexity=perp, 
        eta=lr, 
        max_iter = niter, 
        num_threads=4)
      KLDiv <- TSNE$itercosts %>% .[length(.)]
      OUTPUT <- data.table(cbind(label_DT,TSNE$Y))
      colnames(OUTPUT) <- c('label','x', 'y', 'z')
      toc(log=T)
      L <- list(kl=KLDiv, time=tic.log()[[1]], error=error_message, fig=OUTPUT)
      V <- mapply(toJSON, x=L, force=TRUE, dataframe='rows')
      tic.clearlog()
      return(V %>% toJSON(., raw='base64', force=TRUE))
    }
    
  }
)

app$callback(
  output = list(id='KL-div', property='children')
  ,
  params =  list(
    input(id='custom-container', property='children')
  )
  ,
  function(JSON){
    L <- fromJSON(JSON) 
    return(fromJSON(L[1]))
  }
)

app$callback(
  output = list(id='end-time', property='children')
  ,
  params =  list(
    input(id='custom-container', property='children')
  )
  ,
  function(JSON){
    L <- fromJSON(JSON) 
    return({fromJSON(L[2])})
  }
)
app$callback(
  output = list(id='error-message', property='children')
  ,
  params =  list(
    input(id='custom-container', property='children')
  )
  ,
  function(JSON){
    L <- fromJSON(JSON) 
    return(fromJSON(L[3]))
  }
)
app$callback(
  output = list(id='tsne-3d-plot-custom', property='figure')
  ,
  params =  list(
    input(id='custom-container', property='children')
  )
  ,
  function(JSON){
    L <- fromJSON(JSON) 
    OUTPUT <- fromJSON(L[4])
    p <- plotly::plot_ly(
      OUTPUT, 
      type='scatter3d', 
      x=~x, 
      y=~y, 
      z=~z, 
      color = ~as.factor(label), 
      mode='markers', 
      marker=list(symbol='circle', size=2.5)
    )%>% plotly::layout(.,width=850, height=600)
    return(p)
  }
)

#
########CALLBACKS for custom data ends



# 
# app$run_server(port="8893")

app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))


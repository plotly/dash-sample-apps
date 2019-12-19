appName <- Sys.getenv("DASH_APP_NAME")

if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)

  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}



#Source assets
source("assets/FinancialFunctions.R")

# Load Necessary Packages
library('dash')
library('dashCoreComponents')
library('dashHtmlComponents')
library('plyr')
library('plotly')
library('lubridate')

#################################################################################

app <- Dash$new()


grey_line <- htmlDiv(list(
  htmlHr(className = "greyline")
))


overview <- htmlDiv(list(
  Header()
))


#################################################################################

#Intro Paragraph and Overview

paragraph <- htmlDiv(list(
  htmlDiv(list(
    htmlH4('Product Summary', style = list("color" = "#ffffff")),
    
    htmlP("\
                            As the industry's first index fund for individual investors, \
                            the Calibre Index Fund is a low-cost way to gain diversified exposure \
                            to the U.S. equity market. The fund offers exposure to 500 of the \
                            largest U.S. companies, which span many different industries and \
                            account for about three-fourths of the U.S. stock market's value. \
                            The key risk for the fund is the volatility that comes with its full \
                            exposure to the stock market. Because the Calibre Index Fund is broadly \
                            diversified within the large-capitalization market, it may be \
                            considered a core equity holding in a portfolio.")), className = 'rowrow2'
  )), className = "product"
)


#################################################################################


#Code for the Average Annual Performance Bar Charts



ax <- list(
  title = ""
)

years <- c("1 year", "3 year", "5 year", "10 year", "41 year")
SP_Index <- c(21.83, 11.41, 15.79, 8.50, 0)
Index_Fund <- c(21.67, 11.26, 15.62, 8.37, 11.11)


performance_data <- data.frame(years, SP_Index, Index_Fund)

performance_data$years <-
  factor(performance_data$years, levels = performance_data$years[order(c(1, 2, 3, 4, 5))])

performance <-
  plot_ly(
    performance_data,
    x = ~ years,
    y = ~ Index_Fund,
    type = "bar",
    name = 'Calibre Index Fund',
    width = 340,
    height = 200
  ) %>%
  add_trace(y = ~ SP_Index, name = "S&P Index Fund") %>%
  layout(
    yaxis = ax,
    xaxis = ax,
    colorway = c('#98151B', '#DCDCDC'),
    legend = list(
      x = 1.0,
      y =  0.95,
      orientation = 'h',
      yanchor = "top",
      font = list(size = 9)
    ),
    autosize = FALSE,
    bargap = 0.35,
    hovermode = "closest",
    margin = list(
      r = 0,
      t = 20,
      b = 10,
      l = 10
    )
  )


performance_graph <- htmlDiv(
  list(
    htmlH6(list('Average annual performance'), className = 'subtitle'),
    htmlBr(),
    dccGraph(id= "Annual Performance", figure=performance)
  ), className = "six columns"
)

#################################################################################


#Code for the Hypothetical Growth Line Graphs

hypothetical_yaxis <- list(
  title = "",
  autotick = FALSE,
  tickmode = "array",
  tickvals = c(0, 10000, 20000, 30000),
  showline = TRUE,
  color = 'gray',
  linecolor = toRGB('lightgray')
)


nogrid <- list(
  title = "" ,
  showgrid = FALSE,
  autotick = FALSE,
  tickmode = "array",
  tickvals = c(2010, 2015),
  showline = TRUE,
  color = 'gray',
  linecolor = toRGB('lightgray')
)

x = c(2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018)
y = c(10000, 7500, 9000, 10000, 10500, 11000, 14000, 18000, 19000, 20500, 30000)


hypothetical_data <- data.frame(x, y)



hypothetical_growth <-
  plot_ly(
    hypothetical_data,
    x = x,
    y = y,
    type = "scatter",
    mode = 'lines',
    line = list(color = '#98151B', width = 3),
    name = "Calibre Index Fund Inv",
    height = 200,
    width = 340
  ) %>%
  layout(
    yaxis = hypothetical_yaxis,
    xaxis = nogrid,
    legend = list(orientation = 'h', y = -0.2),
    showlegend = TRUE,
    autosize = FALSE
    ,
    margin = list(
      r = 0,
      t = 20,
      b = 10,
      l = 10
    )
  )



hypothetical_graph <- htmlDiv(
  list(
    htmlH6(list('Hypothetical Growth of $10,000'), className = 'subtitle'),
    htmlBr(),
    dccGraph(id= "Hypothetical Graph", figure=hypothetical_growth)
  ) , className = 'six columns'
)

#################################################################################

#Code for the performance line graphs.


df_graph <- read.csv('data/df_graph.csv')


performance_yaxis <- list(
  title = "",
  autotick = FALSE,
  tickmode = "array",
  tickvals = c(100,200),
  showline = TRUE,
  color = 'gray',
  linecolor = toRGB('lightgray')
)


performance_xaxis <- list(
  title = "" ,
  showgrid = TRUE,
  autotick = FALSE,
  tickmode = "array",
  tickvals = c(2008,2010,2012,2014,2016,2018),
  showline = TRUE,
  color = 'gray',
  linecolor = toRGB('lightgray'),
  type="date",
  
  rangeselector = (list(
    buttons = list(
      list(
        count = 1,
        label = "1Y",
        step = "year",
        stepmode = "backward"
      ),
      list(
        count = 3,
        label = "3Y",
        step = "year",
        stepmode = "backward"
      ),
      list(
        count = 5,
        label = "5Y",
        step = "year",
        stepmode = "backward"
      ),
      list(
        count = 10,
        label = "10Y",
        step = "year",
        stepmode = "backward"
      )
    )
  ))
  
  
)


df_graph$Date <- ymd(as.character(df_graph$Date))

df_graph$Date <- floor_date(df_graph$Date, "month")

df_graph <-
  ddply(
    df_graph,
    "Date",
    summarise,
    Vanguard.500.Index.Fund = mean(Vanguard.500.Index.Fund),
    MSCI.EAFE.Index.Fund..ETF. = mean(MSCI.EAFE.Index.Fund..ETF.)
  )






performance_lines <-
  plot_ly(
    df_graph,
    x = df_graph$Date,
    y = df_graph$Vanguard.500.Index.Fund,
    type = "scatter",
    mode = 'lines',
    line = list(color = '#98151B', width = 3),
    name = "Calibre Index Fund",
    height = 200,
    width = 750
  ) %>%
  add_trace(
    y = df_graph$MSCI.EAFE.Index.Fund..ETF.,
    line = list(color = '#DCDCDC', width = 3),
    name = "MSCI EAFE Index Fund (ETF)"
  ) %>%
  layout(
    yaxis = performance_yaxis,
    xaxis = performance_xaxis,
    autosize = FALSE,
    margin = list(
      r = 40,
      t = 40,
      b = 30,
      l = 40
    )
  )

performance_lines_graph <- htmlDiv(list(
  htmlH6(list('Performance'), className = 'subtitle'),
  htmlBr(),
  dccGraph(id= "Performance_Lines", figure=performance_lines)
) , className = 'twelve columns'
)


#################################################################################

stock_style_graph <- plot_ly(x = list("1"), y = list("1"), hoverinfo = "none",
                             marker = list("opacity" = 0), mode = "markers", name = "B", type = "scatter",
                             width = 200,
                             height = 150) %>%
  layout(title = "",
         annotations = list(
           list(
             "x" = 0.990130093458,
             "y" = 1.00181709504,
             "align" = "left",
             "font" = list("size" = 9),
             "showarrow" = FALSE,
             "text" = "<b>Market<br>Cap</b>",
             "xref" = "x",
             "yref" = "y"
           ),
           list(
             "x"= 1.00001816013,
             "y"= 1.35907755794e-16,
             "font" = list("size" = 9),
             "showarrow" = FALSE,
             "text" = "<b>Style</b>",
             "xref" = "x",
             "yref" = "y"
           )
         ),
         
         hovermode = "closest",
         margin = list("r" = 30, "t" = 20, "b" = 20, "l" = 30),
         
         shapes = list(
           list(
             "fillcolor" = "rgb(127,127,127)",
             "line" = list("color"="rgb(0,0,0)", "width"=2),
             "opacity" = 0.3,
             "type" = "rect",
             "x0" = 0,
             "x1" = 0.33,
             "xref" = "paper",
             "y0" = 0,
             "y1" = 0.33,
             "yref" = "paper"
           ),
           list(
             "fillcolor" = "rgb(127,127,127)",
             "line" = list("color"="rgb(0,0,0)", "dash" ="solid", "width"=2),
             "opacity" = 0.3,
             "type" = "rect",
             "x0" = 0.33,
             "x1" = 0.66,
             "xref" = "paper",
             "y0" = 0,
             "y1" = 0.33,
             "yref" = "paper"
           ),
           list(
             "fillcolor" = "rgb(127,127,127)",
             "line" = list("color"="rgb(0,0,0)", "width"=2),
             "opacity" = 0.3,
             "type" = "rect",
             "x0" = 0.66,
             "x1" = 0.99,
             "xref" = "paper",
             "y0" = 0,
             "y1" = 0.33,
             "yref" = "paper"
           ),
           list(
             "fillcolor" = "rgb(127,127,127)",
             "line" = list("color"="rgb(0,0,0)", "width"=2),
             "opacity" = 0.3,
             "type" = "rect",
             "x0" = 0,
             "x1" = 0.33,
             "xref" = "paper",
             "y0" = 0.33,
             "y1" = 0.66,
             "yref" = "paper"
           ),
           list(
             "fillcolor" = "rgb(127,127,127)",
             "line" = list("color"="rgb(0,0,0)", "width"=2),
             "opacity" = 0.3,
             "type" = "rect",
             "x0" = 0.33,
             "x1" = 0.66,
             "xref" = "paper",
             "y0" = 0.33,
             "y1" = 0.66,
             "yref" = "paper"
           ),
           list(
             "fillcolor" = "rgb(127,127,127)",
             "line" = list("color"="rgb(0,0,0)", "width"=2),
             "opacity" = 0.3,
             "type" = "rect",
             "x0" = 0.66,
             "x1" = 0.99,
             "xref" = "paper",
             "y0" = 0.33,
             "y1" = 0.66,
             "yref" = "paper"
           ),
           list(
             "fillcolor" = "rgb(127,127,127)",
             "line" = list("color"="rgb(0,0,0)", "width"=2),
             "opacity" = 0.3,
             "type" = "rect",
             "x0" = 0,
             "x1" = 0.33,
             "xref" = "paper",
             "y0" = 0.66,
             "y1" = 0.99,
             "yref" = "paper"
           ),
           list(
             "fillcolor" = "#98151B",
             "line" = list("color"="rgb(0,0,0)", "width"=1),
             "opacity" = 0.9,
             "type" = "rect",
             "x0" = 0.33,
             "x1" = 0.66,
             "xref" = "paper",
             "y0" = 0.66,
             "y1" = 0.99,
             "yref" = "paper"
           ),
           list(
             "fillcolor" = "rgb(127,127,127)",
             "line" = list("color"="rgb(0,0,0)", "width"=1),
             "opacity" = 0.9,
             "type" = "rect",
             "x0" = 0.66,
             "x1" = 0.99,
             "xref" = "paper",
             "y0" = 0.66,
             "y1" = 0.99,
             "yref" = "paper"
           )
         ),
         xaxis = list(
           "autorange" = TRUE,
           "range" = list(0.989694747864, 1.00064057995),
           "showgrid" = FALSE,
           "showline" = FALSE,
           "showticklabels" = FALSE,
           "title"  = "<br>",
           "type" = "linear",
           "zeroline" = FALSE
         ),
         yaxis = list(
           "autorange" = TRUE,
           "range" = list(-0.0358637178721, 1.06395696354),
           "showgrid" = FALSE,
           "showline" = FALSE,
           "showticklabels" = FALSE,
           "title"  = "<br>",
           "type" = "linear",
           "zeroline" = FALSE
         )
         
  )

stock_graph <- htmlDiv(list(
  htmlH6(list('Portfolio'), className = 'subtitle'),
  htmlBr(),
  dccGraph(id= "Stock_Graph", figure=stock_style_graph, className = "four columns")))

stock_graph_noheader <- htmlDiv(list(
  dccGraph(id= "Stock_Graph", figure=stock_style_graph, className = "four columns")))

stock_text <- htmlDiv(list(
  htmlLi("Calibre Index Fund seeks to track the performance of\
                     a benchmark index that meaures the investment return of large-capitalization stocks."),
  htmlLi("Learn more about this portfolio's investment strategy and policy.")
) , className="eight columns middle-aligned")



# Fees and Expenses


fees_graph <- plot_ly(x = list("Category Average", "This fund"), y = list("2242", "329"),
                      marker = list("color" = "#98151B"), name = "A", type = "bar",  height = 150,
                      width = 340) %>%
  add_trace(x = list("This Fund"), y = list("1913"), marker = list("color"="#D3D3D3"),
            name = "B", type = "bar") %>%
  layout(annotations = list(
    list(
      "x" = "-0.0111111111111",
      "y" = "2381.92771084",
      "font" = list("color"="rgb(0,0,0)", "size" = 10),
      "showarrow" = FALSE,
      "text" = "$2,242",
      "xref" = "x",
      "yref" = "y"
    ),
    list(
      "x" = "0.995555555556",
      "y" = "509.638554217",
      "font" = list("color"="rgb(0,0,0)", "size" = 10),
      "showarrow" = FALSE,
      "text" = "$329",
      "xref" = "x",
      "yref" = "y"
    ),
    list(
      "x" = "0.995551020408",
      "y" = "1730.32432432",
      "font" = list("color"="rgb(0,0,0)", "size" = 10),
      "showarrow" = FALSE,
      "text" = "You save<br><b>$1,913</b>",
      "xref" = "x",
      "yref" = "y"
    )
  ),
  autosize = FALSE,
  bargap = 0.5,
  barmode = "stack",
  hovermode = "closest",
  margin = list("r" = 40, "t" = 20, "b" = 20, "l" = 40),
  showlegend = FALSE,
  title = "",
  xaxis = list(
    "autorange" = TRUE,
    "mirror" = FALSE,
    "showline" = TRUE,
    "tickfont" = list("size" = 10),
    "title" = "",
    "type" = "category",
    "zeroline" = FALSE
  ),
  yaxis = list(
    "autorange" = FALSE,
    "mirror" = FALSE,
    "nticks" = 3,
    "range" = list(0, 3000),
    "showline" = TRUE,
    "tickfont" = list("size" = 10),
    "title" = "",
    "tickprefix" = "$",
    "type" = "linear",
    "zeroline" = FALSE
  )
  )

fees_bars <- htmlDiv(list(
  htmlStrong("Fees on $10,000 invested over 10 years"),
  htmlBr(),
  dccGraph(id= "Fees", figure=fees_graph)
), className = "six columns"
)




#################################################################################

news <- htmlDiv(list(
  htmlH6(list('Calibre News'), className = 'subtitle'),
  htmlBr(),
  htmlLi('10/25/16    The rise of indexing and the fall of costs'),
  htmlBr(),
  htmlLi("08/31/16    It's the index mutual fund's 40th anniversary: Let the low-cost, passive party begin")
), className = "six columns")

reviews <- htmlDiv(list(
  htmlH6(list('Reviews'), className = 'subtitle'),
  htmlBr(),
  htmlLi('Launched in 1976.'),
  htmlLi('On average, has historically produced returns that have far outpaced the rate of inflation.*'),
  htmlLi("Calibre Quantitative Equity Group, the fund's advisor, is among the world's largest equity index managers."),
  htmlBr(),
  htmlP("Did you know? The fund launched in 1976 as Calibre First Index Investment 
         Trust-the nation's first index fund available to individual investors."),
  htmlBr(),
  htmlP("* The performance of an index is not an exact representation of any particular investment, as you cannot invest directly in an index."),
  htmlBr(),
  htmlP("Past performance is no guarantee of future returns. See performance data current to the most recent month-end.")
), className = "six columns")


fees_text <- htmlDiv(list(
  htmlDiv(list(
    htmlH6(list("Fees"), className = "subtitle"),
    
    htmlBr(list()),
    
    htmlDiv(list(
      
      htmlDiv(list(
        htmlStrong(list("Purchase fee"))
      ), className = "three columns right-aligned"),
      
      htmlDiv(list(
        htmlP(list("None"))
      ), className = "nine columns")
    ), className = "row"),
    
    htmlDiv(list(
      htmlDiv(list(
        htmlStrong(list("Redemption fee"))
      ), className = "three columns right-aligned"),
      
      htmlDiv(list(
        htmlP(list("None"))
      ), className = "nine columns")
    ), className = "row"),
    
    htmlDiv(list(
      htmlDiv(list(
        htmlStrong(list("12b-1 fee"))
      ), className = "three columns right-aligned"),
      
      htmlDiv(list(
        htmlP(list("None"))
      ), className = "nine columns")
    ), className = "row"),
    
    htmlDiv(list(
      htmlDiv(list(
        htmlStrong(list("Account service fee"))
      ), className = "three columns right-aligned"),
      
      htmlDiv(list(
        htmlStrong(list("Nonretirement accounts, traditional IRAs, Roth IRAs, 
                        UGMAs/UTMAs, SEP-IRAs, and education savings accounts (ESAs)")),
        htmlP(list("We charge a $20 annual account service fee for each Calibre Brokerage Account, 
                   as well as each individual Calibre mutual fund holding with a balance of less than 
                   $10,000 in an account. This fee does not apply if you sign up for account access on 
                   Calibre.com and choose electronic delivery of statements, confirmations, and 
                   Calibre fund reports and prospectuses. This fee also does not apply to members 
                   of Flagship SelectT, Flagship®, Voyager Select®, and Voyager® Services.")),
        htmlBr(list()),
        htmlStrong(list("SIMPLE IRAs")),
        htmlP(list("We charge participants a $25 annual account service fee for each fund they hold in their Calibre SIMPLE IRA. 
                   This fee does not apply to members of Flagship Select, Flagship, Voyager Select, and Voyager Services.")),
        htmlBr(list()),
        htmlStrong(list("403(b)(7) plans")),
        htmlP(list("We charge participants a $15 annual account service fee for each fund they hold in their Calibre 403(b)(7) account.
                   This fee does not apply to members of Flagship Select, Flagship, Voyager Select, and Voyager Services.")),
        htmlBr(list()),
        htmlStrong(list("Individual 401(k) plans")),
        htmlP(list("We charge participants a $20 annual account service fee for each fund they hold in their 
                   Calibre Individual 401(k) account. This fee will be waived for all participants in the plan 
                   if at least 1 participant qualifies for Flagship Select, Flagship, Voyager Select, and Voyager Services")),
        htmlBr(list())
        
      ), className = "nine columns")
    ), className = "row")
  ), className = "twelve columns")
), className = "row")

#################################################################################
#Tables


#Fund Facts Table

fund_facts_data <- read.csv("data/df_fund_facts.csv")

fund_facts_table <- generate_table(as.matrix(fund_facts_data))

funds <- htmlDiv(list(
  htmlH6(list('Fund Facts'), className = 'subtitle'),
  htmlBr(),
  fund_facts_table
), className = "six columns")



#Price Performance Table

price_perf_data <- read.csv("data/df_price_perf.csv")

price_perf_table <- generate_table(as.matrix(price_perf_data))

price_perf <- htmlDiv(list(
  htmlH6(list('Price and Performance (%)'), className = 'subtitle'),
  htmlBr(),
  price_perf_table
), className = "six columns")



#Current Prices Table


current_prices_data <- read.csv("data/df_current_prices.csv")

current_prices_table <- generate_table(as.matrix(current_prices_data))

current_prices <- htmlDiv(list(
  htmlH6(list('Current Prices'), className = 'subtitle'),
  htmlBr(),
  current_prices_table
), className = "six columns")



#Historical Prices Table

historical_prices_data <- read.csv("data/df_hist_prices.csv")

historical_prices_table <- generate_table(as.matrix(historical_prices_data))

historical_prices <- htmlDiv(list(
  htmlH6(list('Historical Prices'), className = 'subtitle'),
  htmlBr(),
  historical_prices_table
), className = "six columns")


# AVerage annual returns table
avg_returns_data <- read.csv("data/df_avg_returns.csv")

avg_returns_data[1,1] <- NA

avg_returns_data[3,6] <- NA

avg_returns_table <- generate_table(as.matrix(avg_returns_data))

avg_returns <- htmlDiv(list(
  htmlH6(list('Average annual returns--updated monthly as of 02/28/2018'), className = 'subtitle'),
  htmlBr(),
  avg_returns_table
), className = "twelve columns")

# After Tax Returns Table

after_tax_data <- read.csv("data/df_after_tax.csv")

after_tax_data[1,1] <- NA
after_tax_data[2, 2:6] <- NA
after_tax_data[6, 2:6] <- NA
after_tax_data[8, 2:6] <- NA
after_tax_data[9, 2:6] <- NA
after_tax_data[4,6] <- NA
after_tax_data[5,6] <- NA
after_tax_data[7,6] <- NA

after_tax_table <- generate_table(as.matrix(after_tax_data))

after_tax <- htmlDiv(list(
  htmlH6(list('After-tax returns--updated quarterly as of 12/31/2017'), className = 'subtitle'),
  htmlBr(),
  after_tax_table
), className = "twelve columns")

# Dividends Table

dividends_data <- read.csv("data/df_dividend.csv")



dividends_data$Distribution.Yield <- "-"

dividends_data[2:6,8] <- "-"



dividends_table <- generate_table(remove.factors(dividends_data))


dividends <- htmlDiv(list(
  htmlH6(list('Dividend and capital gains distributions'), className = 'subtitle'),
  htmlBr(),
  dividends_table
), className = "row")


# Realized and Unrealized Gains Table

realized_data <- read.csv("data/df_realized.csv")

realized_table <- generate_table(as.matrix(realized_data))

realized <- htmlDiv(list(
  htmlBr(),
  htmlH6(list('Realized/unrealized gains as of 01/31/2018'), className = 'subtitle'),
  htmlBr(),
  realized_table
), className = "six columns")

# Unrealized appreciation table

unrealized_data <- read.csv("data/df_unrealized.csv")

unrealized_table <- generate_table(as.matrix(unrealized_data))

unrealized <- htmlDiv(list(
  htmlBr(),
  htmlH6(list('Unrealized appreciation/depreciation'), className = 'subtitle'),
  htmlBr(),
  unrealized_table
), className = "six columns")

# Recent Investment Returns

recent_returns_data <- read.csv("data/df_recent_returns.csv")

recent_returns_table <- generate_table(as.matrix(recent_returns_data))

recent_returns <- htmlDiv(list(
  htmlH6(list('Recent Investment Returns'), className = 'subtitle'),
  htmlBr(),
  recent_returns_table
), className = "twelve columns")


# Equity Sector Table

equity_characteristics_data <- read.csv("data/df_equity_char.csv")

equity_characteristics_data[8:12, 3] <- NA


equity_characteristics_table <- generate_table(as.matrix(equity_characteristics_data))


equity_characteristics <- htmlDiv(list(
  htmlH6(list("Equity characteristics as of 01/31/2018"), className = 'subtitle'),
  htmlBr(),
  equity_characteristics_table
), className = "twelve columns")

#Equity Sector Diversification Table

equity_div_data <- read.csv('data/df_equity_diver.csv')

equity_div_table <- generate_table(as.matrix(equity_div_data))


equity_div <- htmlDiv(list(
  htmlH6(list("Equity sector diversification"), className = 'subtitle'),
  htmlBr(),
  equity_div_table
), className = "twelve columns")

# Expenses Table


expenses_data <- read.csv("data/df_expenses.csv")

expenses_table <- generate_table(as.matrix(expenses_data))

expenses <- htmlDiv(list(
  htmlH6(list('Expenses'), className = 'subtitle'),
  htmlBr(),
  expenses_table
), className = "six columns")


# Minimums Table

minimums_data <- read.csv("data/df_minimums.csv")

minimums_table <- generate_table(as.matrix(minimums_data))

minimums <- htmlDiv(list(
  htmlH6(list('Minimums'), className = 'subtitle'),
  htmlBr(),
  minimums_table
), className = "six columns")


######################################################################################

#Second Page First Row

secondpage_firstrow <- htmlDiv(list(
  current_prices,
  historical_prices
), className = "row")



#First Page First Row

firstpage_firstrow <- htmlDiv(list(
  funds,
  performance_graph
), className = "row")

#First Page Second Row

firstpage_secondrow <- htmlDiv(list(
  # htmlBr(),
  hypothetical_graph,
  price_perf
), className = "row")

#News and Reviews
news_reviews <- htmlDiv(list(
  htmlBr(),
  news,
  reviews
), className = "row")


#News and Reviews
news_reviews <- htmlDiv(list(
  htmlBr(),
  news,
  reviews
), className = "row")


# Distributions Row Two
distributions <- htmlDiv(list(
  htmlBr(),
  realized,
  unrealized
), className = "row")

# Portfolio&Stock Page First Row

portfolio_firstrow <- htmlDiv(list(
  stock_graph,
  stock_text
), className = "row")


# RiskReward Image and Layout

risk_reward <- htmlDiv(list(
  htmlH6(list("Risk Potential"), className = 'subtitle'),
  htmlBr(),
  htmlImg(src= 'assets/risk_reward.png', height = "140", width = "332")
), className = "six columns")

white_space <- htmlDiv(list(
  htmlBr()
), className = "six columns")

risk_reward_proper <- htmlDiv(list(
  white_space,
  risk_reward
), className = "row")

# Expenses Row

expenses_row <- htmlDiv(list(
  expenses
), className = "row")

# Fees and Minimums Row

fees_mins_row <- htmlDiv(list(
  minimums,
  fees_bars
), className = "row")


after_tax_row <- htmlDiv(list(
  after_tax
), className = "row")


avg_returns_row <- htmlDiv(list(
  avg_returns
), className = "row")

performance_lines_row <- htmlDiv(list(
  performance_lines_graph
), className = "row")


equity_characteristics_row <- htmlDiv(list(
  equity_characteristics
), className = "row")


equity_div_row <- htmlDiv(list(
  equity_div
), className = "row")


dividends_row <- htmlDiv(list(
  dividends
), className = "row")

#################################################################################

#Define pages to be selected with URL calls.


app$layout(htmlDiv(list(
  
  #URL
  dccLocation(id = 'url', refresh=FALSE),
  
  #Content
  htmlDiv(id='page-content')
)))

index_page <- (htmlDiv(list(
  overview,
  grey_line,
  paragraph,
  htmlDiv(list(
    firstpage_firstrow,
    htmlBr(),
    htmlBr(),
    firstpage_secondrow,
    risk_reward_proper,
    secondpage_firstrow,
    performance_lines_row,
    avg_returns_row,
    after_tax_row,
    stock_graph_noheader,
    stock_text,
    expenses_row,
    fees_mins_row,
    fees_text,
    dividends,
    distributions,
    news_reviews
  ), className='subpagefull')
), className = 'page'))

page_1_layout <- (htmlDiv(list(
  overview,
  grey_line,
  paragraph,
  htmlDiv(list(
    firstpage_firstrow,
    firstpage_secondrow,
    risk_reward_proper,
    htmlDiv(id='page-1-content'),
    htmlBr()
  ), className='subpageone')
), className = 'page'))

page_2_layout <- (htmlDiv(list(
  overview,
  grey_line,
  htmlDiv(list(
    secondpage_firstrow,
    performance_lines_row,
    avg_returns_row,
    after_tax_row,
    recent_returns,
    htmlDiv(id='page-1-content')
  ), className='subpagetwo')
), className = 'page'))

news_layout <- (htmlDiv(list(
  overview,
  grey_line,
  htmlDiv(list(
    news_reviews,
    htmlBr(),
    htmlDiv(id='page-1-content'),
    htmlBr()
  ), className='subpagefive')
), className = 'pagetwo'))


distributions_layout <- (htmlDiv(list(
  overview,
  grey_line,
  htmlDiv(list(
    dividends_row,
    htmlBr(),
    htmlBr(),
    htmlBr(),
    distributions,
    htmlDiv(id='page-1-content'),
    htmlBr()
  ), className='subpagefive')
), className = 'pagetwo'))

portfolio_layout <- (htmlDiv(list(
  overview,
  grey_line,
  htmlDiv(list(
    portfolio_firstrow,
    equity_characteristics_row,
    equity_div_row
  ), className='subpagetwo')
), className = 'page'))


fees_minimums_layout <- htmlDiv(list(
  overview,
  grey_line,
  htmlDiv(list(
    expenses_row,
    fees_mins_row,
    fees_text
  ), className='subpagefour')
), className = 'page')
#################################################################################

app$callback(output = list(id='page-content', property = 'children'),
             params = list(input(id='url', property = 'pathname')),
             display_page <- function(pathname) {
               if (pathname == '/dashr-financial-report/overview') {
                 return(page_1_layout)
               }
               else if (pathname == '/') {
                 return(page_1_layout)
               }
               else if (pathname == '/dashr-financial-report/price-performance') {
                 return(page_2_layout)
               }
               else if (pathname == '/dashr-financial-report/news-and-reviews') {
                 return(news_layout)
               }
               else if (pathname == '/dashr-financial-report/portfolio-management') {
                 return(portfolio_layout)
               }
               else if (pathname == '/dashr-financial-report/distributions') {
                 return(distributions_layout)
               }
               else if (pathname == '/dashr-financial-report/fees') {
                 return (fees_minimums_layout)
               }
               else {
                 return(index_page)
               }
             }
)

if (appName != "") {
  
  app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
  
} else {
  
  app$run_server()
  
}

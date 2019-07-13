

# Generate a Table Function 
generate_table <- function(df, nrows=20)
  
  # function generates a html table from a supplied data frame (df)
  
  #  and number of rows (nrows) to display
{
  
  n <- min(nrows, nrow(df))
  
  rows <- lapply(seq(1, n), function(i) {
    
    htmlTr(children = lapply(as.character(df[i,]), htmlTd))
    
  })
  
  header <- htmlTr(children = lapply(names(df), htmlTh))
  
  htmlTable(
    
    children = c(list(header), rows)
  )
}

# Remove.Factors

remove.factors <- function(x) as.character(levels(x))[x]


#Create a HTML Header:

Header <- function() {
  return(
    htmlDiv(list(
      getLogo(),
      htmlBr(),
      getHeader(),
      htmlDiv(list(
        grey_line
      )),
      htmlBr(),
      getMenu()
    ))
  )
}


getLogo <- function() {
  logo = htmlDiv(list(
    
    htmlDiv(list(
      (htmlImg(src = 'assets/Logo.png', height='40', width='220')
      )), className = 'ten columns padded'),
    
    
    htmlDiv(list(htmlA( href = "/dashr-financial-report/fullversion",
                        htmlButton(list("Full Version"), className = "redbutton"
                        )))),
    
    htmlDiv(list(htmlA( href = 'https://github.com/plotly/dash-sample-apps/dashr-financial-report',
                        htmlButton(list("Learn More"), className = "learnbutton"
                        )))),
    
    htmlDiv(list(htmlA(
      htmlImg(src='assets/dashLogo.png', 
              height='50', width='190'),
      href='https://plot.ly/products/dash/', className = 'logo')))
    
    
  ), className = 'row')
  return(logo)
}

getHeader <- function() {
  header = htmlDiv(list(
    
    htmlDiv(list(
      htmlP(
        'Calibre Financial Index Fund Investor Shares')
    ),
    className = 'maintitle')
  ))
  return(header)
}

getMenu <- function() {
  menu = htmlDiv(list(
    dccLink('Overview   ', href='/dashr-financial-report/overview', className="tab"),
    
    dccLink('Price Performance   ', href='/dashr-financial-report/price-performance', className="tab"),
    
    dccLink('Portfolio & Management   ', href='/dashr-financial-report/portfolio-management', className="tab"),
    
    dccLink('Fees & Minimums   ', href='/dashr-financial-report/fees', className="tab"),
    
    dccLink('Distributions   ', href='/dashr-financial-report/distributions', className="tab"),
    
    dccLink('News & Reviews   ', href='/dashr-financial-report/news-and-reviews', className="tablast")
  ), className = "rowrow")
  return (menu)
}


# 








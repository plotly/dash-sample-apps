appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
  setwd(sprintf("/app/apps/%s", appName))
}

library(dash)
library(dashDaq)
library(dashHtmlComponents)
library(dashCoreComponents)
library(magrittr)
library(purrr)
library(rlang)
library(compiler)
library(plotly)


app <- Dash$new()
font_color <- list(dark='#ffffff', light='#222')
background_color <- list(dark='#2a3f5f', light='#ffffff')
axis_color <- list(dark='#EBF0F8', light='#506784')
marker_color <- list(dark='#f2f5fa', light='#2a3f5f')

theme <- list(
  dark=F,
  primary = '#447EFF',
  secondary = '#D3D3D3',
  detail = '#D3D3D3'
)

init_input <- list(list(
    function_generator = T,
    oscilloscope = T,
    frequency_input = 1E6,
    amplitude_input = 1,
    offset_input = 0,
    function_type = 'SIN'
  ))

##HELPERS

sawtooth <- cmpfun(function(t, width=1){
  w <- width+ numeric(length(t))
  y <- numeric(length(t))
  mask1 <- (w>1 || w<0)
  y[mask1] <- NA
  tmod <- t %%(2*pi)
  mask2 <- (!mask1) && (tmod< 2*w*pi)
  tsub <- tmod[mask2]
  wsub <- w[mask2]
  y[mask2] <- tsub/(pi*wsub)-1
  
  mask3 <- (!mask1) && (!mask2)
  tsub_ <- tmod[mask3]
  wsub_ <- w[mask3]
  y[mask3] <- (pi*(wsub_+1)-tsub_)/(pi*(1-wsub_))
  return(y)
})

header <- cmpfun(function(){
  return(htmlDiv(
    id='header',
    className='banner',
    style=list(backgroundColor='#6682C0'),
    children=list(
      htmlH2(
        "Dash DAQ: Function Generator & Oscilloscope Control Panel",
        style=list(
          color='white', 
          marginLeft='40px',
          display='inline-block',
          'text-align'='center'
          )),
      htmlImg(
        #src="https://s3-us-west-1.amazonaws.com/plotly-tutorials/excel/dash-daq/dash-daq-logo-by-plotly-stripe+copy.png",
        src='assets/dash-logo.png',
        style=list(
          position='relative',
          float='right',
          right='10px',
          height='75px'
        )
        )
    )))
})

knobs <- cmpfun(function(cur_input, cur_tab){
  return(htmlDiv(list(
    ###FREQ INPUT
    daqKnob(
      value=cur_input[[cur_tab]][['frequency_input']],
      id='frequency-input',
      label='Frequency (Hz)',
      labelPosition='bottom',
      size=75,
      color=theme[['primary']],
      scale=list(interval=100000),
      max=2500000,
      min=100000,
      className='four columns'
    ),
    ###AMP INPUT
    daqKnob(
      value=cur_input[[cur_tab]][['amplitude_input']],
      id='amplitude-input',
      label='Amplitude (mV)',
      labelPosition='bottom',
      size=75,
      color=theme[['primary']],
      scale=list(labelInterval=10),
      max=10,
      min=0,
      className='four columns'
    ),
    ###OFFSET
    daqKnob(
      value=cur_input[[cur_tab]][['offset_input']],
      id='offset-input',
      label='Offset (mV)',
      labelPosition='bottom',
      size=75,
      color=theme[['primary']],
      scale=list(labelInterval=10),
      max=10,
      min=0,
      className='four columns'
    )
  ),
  style=list(marginLeft='20%', textAlign='center')))
})
led_displays <- cmpfun(function(cur_input, cur_tab){
  return(htmlDiv(list(
    daqLEDDisplay(
      id='frequency-display',
      size=10,
      value=cur_input[[cur_tab]][['frequency_input']],
      label="Frequency (Hz)",
      labelPosition="bottom",
      color=theme[['primary']],
      style=list(marginBottom= '30px'),
      className='four columns'),
    
    daqLEDDisplay(
      id='amplitude-display',
      size=10,
      value=cur_input[[cur_tab]][['amplitude_input']],
      label="Amplitude (mV)",
      labelPosition="bottom",
      color=theme[['primary']],
      className='four columns'),
    
    daqLEDDisplay(
      id='offset-display',
      size=10,
      value=cur_input[[cur_tab]][['offset_input']],
      label="Offset (mV)",
      labelPosition="bottom",
      color=theme[['primary']],
      className='four columns')
    ), style=list(marginLeft = '20%', textAlign = 'center')))
})

radioitem <- cmpfun(function(cur_input, cur_tab){
  return(dccRadioItems(
    id='function-type',
    options=list(
      list(label= 'Sine', value= 'SIN'),
      list(label= 'Square', value= 'SQUARE'),
      list(label= 'Ramp', value= 'RAMP')
      ),
    value=cur_input[[cur_tab]][['function_type']],
    labelStyle=list(display= 'inline-block'),
    style=list(
      margin= '30px auto 0px auto',
      display = 'flex',
      width = '80%',
      alignItems = 'center',
      justifyContent = 'space-between'
      )))
  })


power_setting_div <- cmpfun(function(cur_inputs, cur_tab){
  if(is.na(cur_inputs)||rlang::is_empty(cur_inputs)){
  cur_inputs = init_input
  }
return(htmlDiv(
  className='row power-settings-tab',
  children=list(
    htmlDiv(
      className='Title',
      children=htmlH3("Power", 
                      id='power-title', 
                      style=list(color=theme[['primary']]))),
    htmlDiv(
      # power-controllers
      list(
        htmlDiv(
            daqPowerButton(
              id='function-generator',
              on=cur_inputs[[cur_tab]][['function_generator']],
              label="Function Generator",
              labelPosition='bottom',
              color=theme[['primary']]),
          className='six columns',
          style=list('margin-bottom'='15px')),
        htmlDiv(
            daqPowerButton(
              id='oscilloscope',
              on=cur_inputs[[cur_tab]][['oscilloscope']],
              label="Oscilloscope",
              labelPosition='bottom',
              color=theme[['primary']]
            ),
          className='six columns',
          style=list('margin-bottom'='15px'))
        ),style=list(margin='15px 0'))
    )))
})


function_setting_div <- cmpfun(function(cur_input, cur_tab){
  if(is.na(cur_input)||rlang::is_empty(cur_input)){
    cur_input = init_input
  }
  return(htmlDiv(
    className='row power-settings-tab',
    children=list(
      htmlDiv(
        id='function-title',
        className='Title',
        style=list(color=theme[['primary']]),
        children=htmlH3("Function")),
      # Knobs
      knobs(cur_input, cur_tab),
      # LED Displays
      led_displays(cur_input, cur_tab),
      # # RadioItems
      radioitem(cur_input, cur_tab)
    )))
})

app$layout(htmlDiv(
  id='main-page',
  className='container',
  children=list(
    #toggle
    htmlDiv(
      id='toggleDiv',
    style=list(
      width='fit-content',
      margin='0 auto'
      ),
    children=list(
      daqToggleSwitch(
        id='toggleTheme',
        style=list(
          # position='absolute',
          # transform='translate(-50%, 20%)',
          # 'z-index'='9999',
          float = 'center'
        ),
        size=30,
        value=F
        )
    )),
    header(),
    htmlDiv(
      children=htmlDiv(
        children=list(
          htmlDiv(
            className='five columns left-panel',
            children=list(
              htmlDiv(
                id='dark-theme-components',
                children=daqDarkThemeProvider(
                  theme=theme,
                  children=list(
                    power_setting_div(NULL, 1),
                    function_setting_div(NULL, 1)
                  )
                )
              ),
              daqColorPicker(
                id='color-picker',
                label='Color Picker',
                value=list(hex='#6682C0'),
                size=164,
                style=list(marginTop = '20px', backgroundColor = 'inherit')
              )
            )
          ),
          #Oscillator Panel - Right
          htmlDiv(
            className='seven columns right-panel',
            children=list(
              htmlDiv(htmlH3("Graph", id='graph-title'),
                      style=list(color=theme[['primary']]),
                      className='Title'),
              dccTabs(
                id='tabs',
                children=list(dccTab(
                  label='Run #1',
                  value='1')),
                value='1',
                className='oscillator-tabs',
                colors=list(
                  border=theme[['primary']],
                  primary=theme[['primary']],
                  background='#f2f2f2'
                )
              ),
              
              htmlDiv(
                className='row oscope-info',
                children=list(
                  htmlDiv(
                    htmlDiv(
                      htmlDiv(
                        id='graph-info',
                        children='-',
                        style=list(
                          border=sprintf('1px solid %s', theme[['primary']])
                          )
                        ),
                    className='row graph-param'), className='six columns'),
                  htmlButton(
                    '+',
                    id='new-tab',
                    n_clicks=0,
                    type='submit',
                    style=list(
                      height='20px',
                      width='20px',
                      padding='2px',
                      lineHeight='10px',
                      float='right',
                      color='inherit'
                    )),
                  htmlButton(
                    '-',
                    id='del-tab',
                    n_clicks=0,
                    type='submit',
                    style=list(
                      height='20px',
                      width='20px',
                      padding='2px',
                      lineHeight='10px',
                      float='right',
                      color='inherit'
                    ))
                  )
              ),
              htmlHr(),
              dccGraph(id='oscope-graph')
            )
          )
        )
      )
    ),
    dccStore(id='control-inputs', data=list())
    # {tabs_number: {value1=x, value2=x}}
)))

#Does not support multiple outputs - need to be updated
app$callback(
  output=list(id='oscilloscope', property='on'),
  params=list(input(id='tabs', property='value'), 
              state(id='control-inputs', property='value'),
              state(id='oscilloscope', property='on'),
              state(id='function-generator', property='on')),
  cmpfun(function(tab_index, cur_inputs, osci_on, func_gen){
    if(!tab_index %in% names(cur_inputs)){
      return(osci_on)
    }
    td <- cur_inputs[[tab_index]]
    return(td[['oscilloscope']])
  })
)

app$callback(
  output=list(id='function-generator', property='on'),
  params=list(input(id='tabs', property='value'), 
              state(id='control-inputs', property='value'),
              state(id='oscilloscope', property='on'),
              state(id='function-generator', property='on')),
  cmpfun(function(tab_index, cur_inputs, osci_on, func_gen){
    if(!tab_index %in% names(cur_inputs)){
      return(func_gen)
    }
    td <- cur_inputs[[tab_index]]
    return(td[['function_generator']])
  })
)


app$callback(
  output=list(id='frequency-input', property='value'),
  params=list(input(id='tabs', property='value'), 
              state(id='control-inputs', property='value'),
              state(id='oscilloscope', property='on'),
              state(id='function-generator', property='on')),
  cmpfun(function(tab_index, cur_inputs, osci_on, func_gen){
    if(!tab_index %in% names(cur_inputs)){
      return(1000000)
    }
    td <- cur_inputs[[tab_index]]
    return(td[['frequency_input']])
  })
)

app$callback(
  output=list(id='amplitude-input', property='value'),
  params=list(input(id='tabs', property='value'), 
              state(id='control-inputs', property='value'),
              state(id='oscilloscope', property='on'),
              state(id='function-generator', property='on')),
  cmpfun(function(tab_index, cur_inputs, osci_on, func_gen){
    if(!tab_index %in% names(cur_inputs)){
      return(1)
    }
    td <- cur_inputs[[tab_index]]
    return(td[['amplitude_input']])
  })
)

app$callback(
  output=list(id='offset-input', property='value'),
  params=list(input(id='tabs', property='value'), 
              state(id='control-inputs', property='value'),
              state(id='oscilloscope', property='on'),
              state(id='function-generator', property='on')),
  cmpfun(function(tab_index, cur_inputs, osci_on, func_gen){
    if(!tab_index %in% names(cur_inputs)){
      return(0)
    }
    td <- cur_inputs[[tab_index]]
    return(td[['offset_input']])
  })
)

app$callback(
  output=list(id='function-type', property='value'),
  params=list(input(id='tabs', property='value'), 
              state(id='control-inputs', property='value'),
              state(id='oscilloscope', property='on'),
              state(id='function-generator', property='on')),
  cmpfun(function(tab_index, cur_inputs, osci_on, func_gen){
    if(!tab_index %in% names(cur_inputs)){
      return('SIN')
    }
    td <- cur_inputs[[tab_index]]
    return(td[['function_type']])
  })
)

app$callback(
  output=list(id='control-inputs', property='data'),
  params=list(
    input(id='oscilloscope',property='on'),
    input(id='function-generator',property='on'),
    input(id='frequency-input',property='value'),
    input(id='amplitude-input',property='value'),
    input(id='offset-input',property='value'),
    input(id='function-type',property='value'),
    state(id='tabs', property='value'),
    state(id='control-inputs', property='data')),
  cmpfun(function(osc_on, fnct_on, freq, amp, offset, wave, sel_tab, cur_inputs){
    cur_inputs[[sel_tab]] = list(oscilloscope=osc_on, 
                                 function_generator=fnct_on, 
                                 frequency_input=freq,
                                 amplitude_input=amp,
                                 offset_input=offset,
                                 function_type=wave)
    return(cur_inputs)
  })
)

# new tab created, not saved to store unless control inputs changes 

app$callback(
  output=list(id='oscope-graph', property='figure'),
  params=list(input(id='control-inputs', property='data'),
              input(id='toggleTheme', property='value'),
              state(id='tabs', property='value')),
  cmpfun(function(cur_inputs, theme_value, tab_index){
    theme_select <- ifelse(theme_value, 'dark', 'light')
    axis <- axis_color[[theme_select]]
    marker <- marker_color[[theme_select]]
    time <- seq(-0.000045, 0.000045, length.out = 1000)
    base_figure <- plot_ly(
                type='scatter',
                x=time,
                y=numeric(length(time)),
                marker=list(color=marker),
                mode='lines')
      
    base_layout <-  purrr::partial(plotly::layout, xaxis=list(title='s',
                             color=axis,
                             titlefont=list(family='Dosis'),
                             size=13),
                  yaxis=list(title='Voltage (mV)',
                             color=axis,
                             range=c(-10, 10),
                             titlefont=list(family='Dosis'),
                             size=13),
                  margin=list(l=40, b=40, t=20, r=50),
                  plot_bgcolor='rgba(0,0,0,0)',
                  paper_bgcolor='rgba(0,0,0,0)'
                  )
    
    if(!(tab_index %in% names(cur_inputs))){
      return(base_figure %>% base_layout)
    }
    tab_data <- cur_inputs[[tab_index]]

    if(rlang::is_empty(tab_data[['oscilloscope']])||is.na(tab_data[['oscilloscope']])){
      base_figure = plot_ly(type='scatter')
      base_layout  %<>% partial(., showticklabels=F, showline=F, zeroline=F)
      return(base_figure %>% base_layout)
    }

    if(rlang::is_empty(tab_data[['function_generator']])||is.na(tab_data[['function_generator']])){
      return(base_figure %>% base_layout)
    }

    if(tab_data[['function_type']]=='SIN'){
      Y <- unlist(lapply(time, function(n){
          return(tab_data[['offset_input']]+tab_data[['amplitude_input']]*sin(2*pi*0.1*tab_data[['frequency_input']]*n))
          }))
      base_figure <- plot_ly(
        type='scatter',
        x=time,
        y=Y,
        line=list(color=marker),
        mode='lines')
      return(base_figure %>% base_layout)
    } else if(tab_data[['function_type']]=='SQUARE'){
      Y <- unlist(lapply(time, function(n){
        return(tab_data[['offset_input']]+tab_data[['amplitude_input']]*sign(sin(2*pi*tab_data[['frequency_input']]*0.1*n)))
      }))
      base_figure <- plot_ly(
        type='scatter',
        x=time,
        y=Y,
        line=list(color=marker),
        mode='lines')
      return(base_figure %>% base_layout)
    } else if(tab_data[['function_type']]=='RAMP'){
      Y <- abs(sawtooth(2*pi*tab_data[['frequency_input']]/10 * time))*tab_data[['amplitude_input']]
      Y <- tab_data[['offset_input']] + 2 * Y - tab_data[['amplitude_input']]
      base_figure <- plot_ly(
        type='scatter',
        x=time,
        y=Y,
        line=list(color=marker),
        mode='lines')
      return(base_figure %>% base_layout)
    } else {
      return(base_figure %>% base_layout)
    }
  })
)

app$callback(
  output=list(id='graph-info', property='children'),
  params=list(input(id='control-inputs', property='data'),
              input(id='toggleTheme', property='value'),
              state(id='tabs', property='value')),
  cmpfun(function(cur_inputs, theme_value, tab_index){
    tab_data <- cur_inputs[[tab_index]]
    return(sprintf('%s | %.3f Hz | %.3f mV | %.3f mV', tab_data[['function_type']], tab_data[['frequency_input']], tab_data[['amplitude_input']], tab_data[['offset_input']]))
  })
)

app$callback(
  output=list(id='frequency-display', property='value'),
  params=list(input(id='frequency-input', property='value')),
  cmpfun(function(val){return(round(val, 3))})
)

app$callback(
  output=list(id='amplitude-display', property='value'),
  params=list(input(id='amplitude-input', property='value')),
  cmpfun(function(val){return(round(val, 3))})
)

app$callback(
  output=list(id='offset-display', property='value'),
  params=list(input(id='offset-input', property='value')),
  cmpfun(function(val){return(round(val, 3))})
)

app$callback(
  output=list(id='dark-theme-components', property='children'),
  params=list(
    input(id='toggleTheme', property='value'),
    input(id='color-picker', property='value'),
    state(id='control-inputs', property='data'),
    state(id='tabs', property='value')
    ),
  cmpfun(function(turn_dark, color_pick, cur_inputs, cur_tab_value){
    theme[['dark']] <- turn_dark
    if(rlang::is_empty(color_pick)||is.na(color_pick)){
      theme[['primary']] <- color_pick[['hex']]
    }
    return(
    daqDarkThemeProvider(
      theme=theme, 
      children=list(
        power_setting_div(cur_inputs, cur_tab_value),
        function_setting_div(cur_inputs, cur_tab_value)
      ))
    )
  })
)

app$callback(
  output=list(id='power-title', property='style'),
  params=list(input(id='color-picker', property='value')),
   cmpfun(function(color){
    return(list('color'=color[['hex']]))
  })
)
app$callback(
  output=list(id='function-title', property='style'),
  params=list(input(id='color-picker', property='value')),
  cmpfun(function(color){
    return(list('color'=color[['hex']]))
  })
)
app$callback(
  output=list(id='graph-title', property='style'),
  params=list(input(id='color-picker', property='value')),
  cmpfun(function(color){
    return(list('color'=color[['hex']]))
  })
)

app$callback(
  output=list(id='graph-info', property='style'),
  params=list(input(id='color-picker', property='value')),
  cmpfun(function(color){
    return(list(border=sprintf('1px solid %s', color[['hex']]),'color'='inherit'))
  })
)

app$callback(
  output=list(id='tabs', property='colors'),
  params=list(input(id='color-picker', property='value')),
  cmpfun(function(color){
    return(list('border'=color[['hex']], primary=color[['hex']], background='#f2f2f2'))
  })
)

app$callback(
  output=list(id='header', property='style'),
  params=list(input(id='color-picker', property='value')),
  cmpfun(function(color){
    return(list('backgroundColor'=color[['hex']]))
  })
)

app$callback(
  output=list(id='function-generator', property='color'),
  params=list(input(id='color-picker', property='value')),
  cmpfun(function(color){
    return(color[['hex']])
  })
)

app$callback(
  output=list(id='oscilloscope', property='color'),
  params=list(input(id='color-picker', property='value')),
  cmpfun(function(color){
    return(color[['hex']])
  })
)

app$callback(
  output=list(id='frequency-input', property='color'),
  params=list(input(id='color-picker', property='value')),
  cmpfun(function(color){
    return(color[['hex']])
  })
)

app$callback(
  output=list(id='amplitude-input', property='color'),
  params=list(input(id='color-picker', property='value')),
  cmpfun(function(color){
    return(color[['hex']])
  })
)

app$callback(
  output=list(id='offset-input', property='color'),
  params=list(input(id='color-picker', property='value')),
  cmpfun(function(color){
    return(color[['hex']])
  })
)

app$callback(
  output=list(id='frequency-display', property='color'),
  params=list(input(id='color-picker', property='value')),
  cmpfun(function(color){
    return(color[['hex']])
  })
)

app$callback(
  output=list(id='amplitude-display', property='color'),
  params=list(input(id='color-picker', property='value')),
  cmpfun(function(color){
    return(color[['hex']])
  })
)

app$callback(
  output=list(id='offset-display', property='color'),
  params=list(input(id='color-picker', property='value')),
  cmpfun(function(color){
    return(color[['hex']])
  })
)

############################ DARK THEME CALLBACKS
app$callback(
  output=list(id='new-tab', property='style'),
  params=list(input(id='toggleTheme', property='value'),
              state(id='new-tab', property='style')),
  cmpfun(function(turn_dark, cur_style){
    if(turn_dark){
      cur_style[['backgroundColor']] <- background_color[['dark']]
      return(cur_style)
    }else{
      cur_style[['backgroundColor']] <- background_color[['light']]
      return(cur_style)
    }
  })
)
app$callback(
output=list(id='del-tab', property='style'),
params=list(input(id='toggleTheme', property='value'),
            state(id='del-tab', property='style')),
cmpfun(function(turn_dark, cur_style){
  if(turn_dark){
    cur_style[['backgroundColor']] <- background_color[['dark']]
    return(cur_style)
  }else{
    cur_style[['backgroundColor']] <- background_color[['light']]
    return(cur_style)
  }
})
)

app$callback(
  output = list(id='main-page', property='style'),
  params = list(input(id='toggleTheme', property='value')),
  cmpfun(function(turn_dark){
    if(turn_dark){
      return(list(backgroundColor=background_color[['dark']], color=font_color[['dark']]))
    }else{
      return(list(backgroundColor=background_color[['light']], color=font_color[['light']]))
    }
  })
)

####### GENERATING TABS
app$callback(
  output=list(id='tabs', property='children'),
  params=list(input(id='new-tab', property='n_clicks'),
              state(id='control-inputs', 'data')),
  cmpfun(function(n_clicks, cur_inputs){
    return(lapply(1:(n_clicks), function(n){return(dccTab(label=sprintf('Run #%d', n), value=sprintf('%d',n)))}))
  })
)
####### DELETING TABS
app$callback(
  output=list(id='new-tab', property='n_clicks'),
  params=list(input(id='del-tab', property='n_clicks'),
              state(id='new-tab', property='n_clicks')),
  cmpfun(function(d, n){
    return(max(n-1, 1))
  })
)

####### LAUNCH THE APP
app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050))
# Reusable components to avoid too many repetitive inline style parameters

card <- function(children, ...){
  htmlSection(
    children,
    style = list(
      padding = "20",
      margin = "5",
      color = "#a5b1cd",
      "user-select" = "none",
      "-moz-user-select" = "none",
      "-webkit-user-select" = "none",
      "-ms-user-select" = "none"
    )
  )
}

namedDropdown <- function(name, ...){
  htmlDiv(
    style = list(margin = "10px 0px"),
    children = list(
      htmlP(
        children = sprintf("%s:", name),
        style = list(marginLeft = "3px")
      ),
      dccDropdown(...)
    )
  )
}

namedSlider <- function(name, ...){
  htmlDiv(
    style = list(padding = "5px 10px 25px"),
    children = list(
      htmlP(sprintf("%s:", name)),
      htmlDiv(
        style = list(marginLeft = "6px"),
        children = dccSlider(...)
      )
    )
  )
}

formattedSlider <- function(...){
  htmlDiv(
    style = list(padding = "5px 10px 25px"),
    children = dccSlider(...)
  )
}

namedRadioItems <- function(name, ...){
  htmlDiv(
    style = list(padding = "20px 10px 25px 4px"),
    children = list(
      htmlP(children = sprintf("%s:", name)),
      dccRadioItems(...)
    )
  )
}

# AUTO GENERATED FILE - DO NOT EDIT

'DashNgl'DashNgl <- function(id=NULL, viewportStyle=NULL, stageParameters=NULL, data=NULL) {
    
    props <- list(id=id, viewportStyle=viewportStyle, stageParameters=stageParameters, data=data)
    if (length(props) > 0) {
        props <- props[!vapply(props, is.null, logical(1))]
    }
    component <- list(
        props = props,
        type = 'DashNgl',
        namespace = 'dash_ngl',
        propNames = c('id', 'viewportStyle', 'stageParameters', 'data'),
        package = 'dashNgl'
        )

    structure(component, class = c('dash_component', 'list'))
}

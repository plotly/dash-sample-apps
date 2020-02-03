library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(dashBio)
library(dashDaq)
library(bio3d)
library(jsonlite)

appName <- Sys.getenv("DASH_APP_NAME")
if (appName != ""){
  pathPrefix <- sprintf("/%s/", appName)
  
  Sys.setenv(DASH_ROUTES_PATHNAME_PREFIX = pathPrefix,
             DASH_REQUESTS_PATHNAME_PREFIX = pathPrefix)
}

############################# set color list #############################
ATOM_COLOR_DICT <- list(
  "C" = "#c8c8c8",
  "H" = "#ffffff",
  "N" = "#8f8fff",
  "S" = "#ffc832",
  "O" = "#f00000",
  "F" = "#ffff00",
  "P" = "#ffa500",
  "K" = "#42f4ee",
  "G" = "#3f3f3f"
)

CHAIN_COLOR_DICT <- list(
  "A" = "#320000",
  "B" = "#8a2be2",
  "C" = "#ff4500",
  "D" = "#00bfff",
  "E" = "#ff00ff",
  "F" = "#ffff00",
  "G" = "#4682b4",
  "H" = "#ffb6c1",
  "I" = "#a52aaa",
  "J" = "#ee82ee",
  "K" = "#75FF33",
  "L" = "#FFBD33",
  "M" = "#400040",
  "N" = "#004000",
  "O" = "#008080",
  "P" = "#008080",
  "x" = "#9c6677",
  "Y" = "#b7c5c8"
)

RESIDUE_COLOR_DICT <- list(
  'ALA' = '#C8C8C8',
  'ARG' = '#145AFF',
  'ASN' = '#00DCDC',
  'ASP' = '#E60A0A',
  'CYS' = '#E6E600',
  'GLN' = '#00DCDC',
  'GLU' = '#E60A0A',
  'GLY' = '#EBEBEB',
  'HIS' = '#8282D2',
  'ILE' = '#0F820F',
  'LEU' = '#0F820F',
  'LYS' = '#145AFF',
  'MET' = '#E6E600',
  'PHE' = '#3232AA',
  'PRO' = '#DC9682',
  'SER' = '#FA9600',
  'THR' = '#FA9600',
  'TRP' = '#B45AB4',
  'TYR' = '#3232AA',
  'VAL' = '#0F820F',
  'ASX' = '#FF69B4',
  'GLX' = '#FF69B4',
  'A' = '#A0A0FF',
  'DA' = '#A0A0FF',
  'G' = '#FF7070',
  'DG' = '#FF7070',
  'I' = '#80FFFF',
  'C' = '#FF8C4B',
  'DC' = '#FF8C4B',
  'T' = '#A0FFA0',
  'DT' = '#A0FFA0',
  'U' = '#FF8080'
)

tmp <- setNames(
  list(  # pylint: disable=invalid-name
    c('GLY', 'ALA', 'LEU', 'ILE', 'VAL', 'MET', 'PRO'),
    c('ASN', 'GLN', 'SER', 'THR', 'CYS'),
    c('ASP', 'GLU'),
    c('LYS', 'ARG', 'HIS'),
    c('TRP', 'TYR', 'PHE'),
    c('A', 'G', 'DA', 'DG'),
    c('DT', 'DC', 'U', 'I', 'C')
  ),
  c(toupper('hydrophobic'), toupper('polar'), toupper('acidic'), 
    toupper('basic'), toupper('aromatic'), toupper('purine'),
    toupper('pyrimidine'))
)
tmp_vec <- unlist(tmp)

RESIDUE_TYPE_COLOR_DICT <- setNames(
  lapply(names(tmp_vec),
         function(name) {
           if(grepl('hydrophobic', name, ignore.case = TRUE)) {
             '#00ff80'
           } else if(grepl('polar', name, ignore.case = TRUE)) {
             '#ff00bf'
           } else if(grepl('acidic', name, ignore.case = TRUE)) {
             '#ff4000'
           } else if(grepl('basic', name, ignore.case = TRUE)) {
             '#0040ff'
           } else if(grepl('aromatic', name, ignore.case = TRUE)) {
             '#ffff00'
           } else if(grepl('purine', name, ignore.case = TRUE)) {
             '#A00042'
           } else if(grepl('pyrimidine', name, ignore.case = TRUE)) {
             '#4F4600'
           } else NULL
         }),
  tmp_vec
)

############################# helper functions #############################
set_custom_colors <- function(custom_colors, custom_colors_name, custom_colors_list, ...) {
  name <- gsub('[[:digit:]]+', '', names(custom_colors[[custom_colors_name]]))
  for(i in 1:length(name)) {
    if(custom_colors_name == 'residue_type') {
      args <- list(...)
      name_i <- args$tmp[[name[i]]]
    } else name_i <- name[i]
    custom_colors_list[
      match(name_i, names(custom_colors_list))
      ] <- custom_colors[[custom_colors_name]][[i]]
  }
  custom_colors_list
}

set_styles <- function(atoms, types, mol_style, custom_colors_list, 
                       default_ATOM_color = '#BEA06E',
                       default_HETATM_color = '#330000',
                       ATOM_COLOR_DICT = ATOM_COLOR_DICT) {
  
  len <- length(atoms)
  lapply(1:len,
         function(i){
           
           atom <- atoms[[i]]
           type <- types[i]
           
           if(type == 'ATOM') {
             
             residue_name <- gsub('[[:digit:]]+', '', atom$residue_name)
             id <- match(residue_name, names(custom_colors_list))
             color <- if(length(id) > 0 && !is.na(id)) {
               custom_colors_list[[id]][1]
             } else default_color
             
             list(
               color = color, 
               visualization_type = mol_style
             )
           } else if(type == 'HETATM') {
             
             id <- match(atom$element, names(ATOM_COLOR_DICT))
             color <- if(length(id) > 0 && !is.na(id)) {
               ATOM_COLOR_DICT[[id]][1]
             } else default_HETATM_color
             
             list(
               color = color, 
               visualization_type = 'stick'
             )
           } else NULL
         })
}

create_style <- function(mol_style, color_style, custom_colors, atoms, types) {
  
  if(length(custom_colors) > 0) {
    switch(names(custom_colors), 
           "residue" = {
             RESIDUE_COLOR_DICT <- set_custom_colors(custom_colors, 
                                                     custom_colors_name = 'residue', 
                                                     custom_colors_list = RESIDUE_COLOR_DICT)
           },
           "residue_type" = {
             RESIDUE_TYPE_COLOR_DICT <- set_custom_colors(custom_colors, 
                                                          custom_colors_name = 'residue_type', 
                                                          custom_colors_list =  RESIDUE_TYPE_COLOR_DICT,
                                                          tmp = tmp)
           },
           "atom" = {
             ATOM_COLOR_DICT <- set_custom_colors(custom_colors, 
                                                  custom_colors_name = 'atom', 
                                                  custom_colors_list =  ATOM_COLOR_DICT)
           },
           "chain" = {
             CHAIN_COLOR_DICT <- set_custom_colors(custom_colors, 
                                                   custom_colors_name = 'chain', 
                                                   custom_colors_list =  CHAIN_COLOR_DICT)
           }
    )
  } 
  
  styles <- switch(
    color_style,
    "residue" = set_styles(atoms = atoms, 
                           types = types, 
                           mol_style = mol_style, 
                           custom_colors_list = RESIDUE_COLOR_DICT,
                           ATOM_COLOR_DICT = ATOM_COLOR_DICT),
    "residue_type" = set_styles(atoms = atoms, 
                                types = types, 
                                mol_style = mol_style, 
                                custom_colors_list = RESIDUE_TYPE_COLOR_DICT,
                                ATOM_COLOR_DICT = ATOM_COLOR_DICT),
    "chain" = set_styles(atoms = atoms, 
                         types = types, 
                         mol_style = mol_style, 
                         custom_colors_list = CHAIN_COLOR_DICT,
                         ATOM_COLOR_DICT = ATOM_COLOR_DICT),
    "atom" = set_styles(atoms = atoms, 
                        types = types, 
                        mol_style = mol_style, 
                        custom_colors_list = ATOM_COLOR_DICT,
                        ATOM_COLOR_DICT = ATOM_COLOR_DICT)
  )
  # remove NULL
  Filter(Negate(is.null), styles)
}
############################# dash apps #############################
header_colors <- function() {
  list(
    bg_color = '#e7625f',
    font_color = 'white'
  )
}

DATAPATH <- './sample_data/molecule3d_'

data_info <- setNames(
  list(
    list(
      name = 'Measles Nucleocapsid',
      description = dccMarkdown("
The measles nucleoprotein forms a large helical complex with
RNA... It is thought to chaperone the process of replication and
transcription by providing a site ready for binding of the
polymerase/phosphoprotein complex while a new RNA chain is being
built.
The structure includes the stable core domain of the
nucleoprotein and a strand of RNA, but the flexible tail of
nucleoprotein was removed in this study.
"),
      link = 'http://pdb101.rcsb.org/motm/231'
    ),
    list(
      name = 'a-cobratoxin-AChBP complex',
      description = dccMarkdown("
The crystal structure of the snake long alpha-neurotoxin,
alpha-cobratoxin, bound to the pentameric
acetylcholine-binding protein (AChBP) from Lymnaea
stagnalis...
The structure unambiguously reveals the positions and
orientations of all five three-fingered toxin molecules
inserted at the AChBP subunit interfaces and the
conformational changes associated with toxin binding.
"),
      link = 'https://www.rcsb.org/structure/1yi5'
    ),
    list(
      name = 'Calcium ATPase',
      description = dccMarkdown("The calcium pump allows muscles to relax after muscle 
contraction. The pump is found in the membrane of the
sarcoplasmic reticulum. In some cases, it is so plentiful that
it may make up 90% of the protein there. Powered by ATP, it
pumps calcium ions back into the sarcoplasmic reticulum,
reducing the calcium level around the actin and myosin
filaments and allowing the muscle to relax.

The structure has a big domain poking out on the outside
of the sarcoplasmic reticulum, and a region that is embedded
in the membrane, forming a tunnel to the other side.
"),
      link = 'http://pdb101.rcsb.org/motm/51'
    ),
    list(
      name = 'DNA',
      description = dccMarkdown("DNA is read-only memory, archived safely inside cells. Genetic
information is stored in an orderly manner in strands of
DNA. DNA is composed of a long linear strand of millions of
nucleotides, and is most often found paired with a partner
strand. These strands wrap around each other in the familiar
double helix...
"),
      link = 'http://pdb101.rcsb.org/motm/23'
    ),
    list(
      name = 'T7 RNA Polymerase',
      description = dccMarkdown("RNA polymerase is a huge factory with many moving parts. 
The constituent proteins form a machine that surrounds DNA
strands, unwinds them, and builds an RNA strand based on the
information held inside the DNA. Once the enzyme gets started,
RNA polymerase marches confidently along the DNA copying RNA
strands thousands of nucleotides long.

...

This structure includes a very small RNA polymerase that is
made by the bacteriophage T7... a small transcription
bubble, composed of two DNA strands and an RNA strand, is
bound in the active site. "),
      link = 'http://pdb101.rcsb.org/motm/40'
    )
  ), c(paste0(DATAPATH, "4uft.pdb"), paste0(DATAPATH, "1yi5.pdb"),
       paste0(DATAPATH, "1su4.pdb"), paste0(DATAPATH, "1bna.pdb"),
       paste0(DATAPATH, "1msw.pdb"))
)

##### main code ####
app <- Dash$new()

tab1 <-   dccTab(
  label = 'About',
  value = 'what-is',
  children = htmlDiv(
    className = 'control-tab',
    children = list(
      htmlH4(className = 'what-is', children = 'What is Molecule3D?'),
      htmlP('Molecule3D is a visualizer that allows you to 
            view biomolecules in multiple representations: 
            sticks, spheres, and cartoons.'),
      htmlP('You can select a preloaded structure, or upload your own, 
            in the "Data" tab. A sample structure is also 
            available to download.'),
      htmlP('In the "View" tab, you can change the style and 
            coloring of the various components of your molecule.')
    )
  )
)

tab2 <- dccTab(
  label = 'Data',
  value = 'upload-download',
  children = htmlDiv(
    className = 'control-tab',
    children = list(
      htmlDiv(
        title = 'download a sample data file to view',
        children = list()
      ),
      htmlDiv(
        title = 'Select molecule to view',
        # className = "app-controls-block",
        children = list(
          htmlDiv(className = 'app-controls-name',
                  children = 'Select structure',
                  style = list(width = '120px')),
          htmlDiv(className = 'app-controls-name',
                  children = ' ',
                  style = list(width = '40px')),
          dccDropdown(
            id = 'dropdown-demostr',
            options = list(
              list(
                label = 'Measles nucleocapsid',
                value = paste0(DATAPATH, "4uft.pdb")
              ),
              list(
                label = 'a-cobratoxin-AChBP complex',
                value = paste0(DATAPATH, "1yi5.pdb")
              ),
              list(
                label = 'Calcium ATPase',
                value = paste0(DATAPATH, "1su4.pdb")
              ),
              list(
                label = 'DNA',
                value =  paste0(DATAPATH, "1bna.pdb")
              ),
              list(
                label = 'T7 RNA polymerase',
                value =  paste0(DATAPATH, "1msw.pdb")
              )
            ),
            value = paste0(DATAPATH, "1bna.pdb"),
            style = list(width = '130px')
          )
        )
      ),
      htmlBr(),
      htmlDiv(
        title = 'Upload biomolecule to view here',
        # className = 'app-controls-block',
        id = 'mol3d-upload-container', 
        children = list(
          dccUpload(
            id = 'mol3d-upload-data',
            className = 'control-upload',
            children = htmlDiv(
              list('Drag and drop or click to upload a file.')
            ),
            style = list(
              width = '100%', 
              height = '80px', 
              lineHeight = '80px',
              borderWidth = '1px', 
              borderStyle = 'dashed',
              borderRadius = '5px', 
              textAlign = 'center'
            )
          ),
          htmlBr(),
          htmlA(
            htmlButton(
              "Download sample structure",
              id = "mol3d-download-sample-data",
              className = 'control-download',
              style = list(
                width = '100%'
              )
            ),
            href = 'assets/sample_data/molecule3d_2mru.pdb',
            download = '2mru.pdb'
          ),
          htmlBr(),
          htmlBr(),
          htmlDiv(id = 'mol3d-data-info')
        )
      )
    )
  )
)

tab3 <- dccTab(
  label = 'View',
  value = 'view-options',
  children = htmlDiv(
    className = 'control-tab', 
    children = list(
      # Textarea container to display the selected atoms
      htmlDiv(
        title = 'view information about selected atoms of biomolecule',
        # className = "app-controls-block",
        id = "mol3d-selection-display",
        children = list(
          htmlP(
            "Selection",
            style= list(
              'font-weight' = 'bold',
              'margin-bottom' = '10px'
            )
          ),
          htmlDiv(id='mol3d-selection-output')
        )
      ),
      
      htmlBr(),
      # Dropdown to select chain representation
      # (sticks, cartoon, sphere)
      htmlDiv(
        title = 'select style for molecule representation',
        # className = "app-controls-block",
        id = 'mol3d-style',
        children = list(
          htmlP(
            'Style',
            style = list(
              'font-weight' = 'bold',
              'margin-bottom' = '10px'
            )
          ),
          dccDropdown(
            id = 'dropdown-styles',
            options = list(
              list(label = 'Sticks', value = 'stick'),
              list(label = 'Cartoon', value = 'cartoon'),
              list(label = 'Spheres', value = 'sphere')
            ),
            value = 'cartoon'
          )
        )
      ),
      htmlBr(),
      # Dropdown to select color of representation
      htmlDiv(
        title = 'select color scheme for viewing biomolecule',
        # className = "app-controls-block",
        id = 'mol3d-style-color',
        children = list(
          htmlP(
            'Color',
            style = list(
              'font-weight' = 'bold',
              'margin-bottom' = '10px'
            )
          ),
          dccDropdown(
            id = 'dropdown-style-color',
            options = list(
              list(label = 'Atom',
                   value = 'atom'),
              list(label = 'Residue identity',
                   value = 'residue'),
              list(label = 'Residue type',
                   value = 'residue_type'),
              list(label = 'Chain',
                   value = 'chain')
            ),
            value = 'residue',
            style = list(width = "130px")
          ),
          dccDropdown(
            id = 'mol3d-coloring-key',
            options = list(),
            style = list(width = "100px")
          ) 
        )
      ),
      htmlBr(),
      htmlDiv(
        title = 'Customize molecule coloring.',
        # className = "app-controls-block",
        children = list(
          htmlP(
            id = 'mol3d-customize-coloring',
            style = list(
              'font-weight' = 'bold',
              'margin-bottom' = '10px'
            )
          ),
          daqColorPicker(
            id = 'mol3d-coloring-value',
            size = 315
          )
        )
      )
    )
  )
)

header <- htmlDiv(
  id = "app-page-header",
  style = list(
    width = "100%",
    background = header_colors()[["bg_color"]],
    color = header_colors()[["font_color"]]
  ),
  children = list(
    htmlA(
      id = "dashbio-logo",
      children = list(
        htmlImg(src='assets/plotly-dash-bio-logo.png', height = '36', width = '190',
                style = list('top' = '10', 'margin-left' = '10px'))
      ),
      href = "/Portal"
    ),
    htmlH2("Molecule3D"),
    htmlA(
      id = "gh-link",
      children = list("View on GitHub"),
      href = "https://github.com/plotly/dash-sample-apps/blob/master/apps/dashr-molecule3d/app.R",
      style = list(color = "white", border = "solid 1px white")
    ),
    htmlImg(
      src = "assets/GitHub-Mark-Light-64px.png"
    )
  )
)

#### dash app layout ####
app$layout(
  header,
  htmlBr(),
  htmlDiv(
    id = 'mol3d-body',
    className = 'app-body',
    children = list(
      htmlDiv(
        id = 'mol3d-control-tabs',
        className = 'control-tabs',
        children = list(
          dccTabs(
            id = 'mol3d-tabs',
            value = 'what-is',
            children = list(
              tab1,
              tab2,
              tab3
            )
          )
        )
      ),
      dccLoading(
        htmlDiv(
          id = 'mol3d-biomolecule-viewer',
          children = list()
        )
      ),
      dccStore(
        id = 'mol3d-color-storage',
        data = list()
      )
    )
  )
)

#### dash app callback ####
app$callback(
  output = list(id = 'mol3d-data-info', property = 'children'),
  params = list(input(id = 'dropdown-demostr', property = 'value')),
  function(molecule_selected) {

    if(molecule_selected %in% names(data_info) && !is.null(unlist(molecule_selected))) {
      mol <- data_info[[molecule_selected]]
      
      list(
        htmlH4(mol$name),
        mol$description,
        htmlA(
          '(source)',
          href = mol$link
        ) 
      )
    } else ""
  }
)

app$callback(
  output = list(id = 'dropdown-demostr', property = 'value'),
  params = list(input(id = 'mol3d-upload-data', property = 'contents'),
                state(id = 'dropdown-demostr', property = 'value')),
  function(upload_content, dem) {
    if(is.null(unlist(upload_content))) {
      dem
    } else list()
  }
)

atom_color <- names(ATOM_COLOR_DICT)
residue_color <- names(RESIDUE_COLOR_DICT)
residue_type_color <- names(tmp)
chain_color <- names(CHAIN_COLOR_DICT)

app$callback(
  output = list(id = 'mol3d-coloring-key', property = 'options'),
  params = list(input(id = 'dropdown-style-color', property = 'value')),
  function(mol_style) {
    switch(mol_style,
           "atom" = lapply(atom_color, function(col) list(label = col, value = col)),
           "residue" = lapply(residue_color, function(col) list(label = col, value = col)),
           "residue_type" = lapply(residue_type_color, function(col) list(label = col, value = col)),
           "chain" = lapply(chain_color, function(col) list(label = col, value = col))
    )
  }
)

app$callback(
  output = list(id = 'mol3d-color-storage', property = 'data'),
  params = list(input(id = 'mol3d-coloring-value', property = 'value'),
                input(id = 'dropdown-style-color', property = 'value'),
                state(id = 'mol3d-coloring-key', property = 'value'),
                state(id = 'mol3d-color-storage', property = 'data')),
  function(color_value, color_style, color_key, current) {
    
    if(is.null(unlist(color_style))) return(list())
    
    if(is.null(unlist(color_key)) || is.null(unlist(color_value))) return(current)
    
    if(!is.null(names(current))){
      if(color_style != names(current)) return(list())
    }
    current[[color_style]][[color_key]] <- color_value[['hex']]
    current[[color_style]] <- setNames(as.list(current[[color_style]]),
                                       names(current[[color_style]]))
    
    current
  }
)

get_pdb <- function(decoded_pdb, type = c("ATOM", "HETATM")) {
  data <- strsplit(decoded_pdb, "\n")[[1]]
  get_ATOM <- lapply(data,
                     function(x){
                       if(any(sapply(type, function(y) grepl(y, x)))) {
                  
                         info <- Filter(f = function(y) y != "", strsplit(x, " ")[[1]])
                         if(info[1] %in% type) {
                           info
                         } else NULL
                       } else NULL
                     })
  get_ATOM <- setNames(
    as.data.frame(
      do.call(
        rbind,
        Filter(Negate(is.null), get_ATOM)
      ),
      stringsAsFactors = FALSE
    ), c('type', 'eleno', 'elety', 'resid', 'chain', 
         'resno', 'x', 'y', 'z', 'o', 'b', 'elesy')
  )

  hide <- lapply(c('eleno', 'resno', 'x', 'y', 'z', 'o', 'b'), 
                 function(column) {
                   get_ATOM[[column]] <<- as.numeric(get_ATOM[[column]])
                 }
  )

  list(
    atom = get_ATOM
  )
}

app$callback(
  output = list(id = 'mol3d-biomolecule-viewer', property = 'children'),
  params = list(
    input(id = 'mol3d-upload-data', property = 'contents'),
    input(id = 'dropdown-demostr', property = 'value'),
    input(id = 'dropdown-styles', property = 'value'),
    input(id = 'dropdown-style-color', property = 'value'),
    input(id = 'mol3d-color-storage', property = 'modified_timestamp'),
    state(id = 'mol3d-color-storage', property = 'data')
  ),
  function(contents, demostr, mol_style, color_style, wt, custom_colors) {
    
    bonds <- list()

    if(length(unlist(demostr)) > 0) {
      
      pdbData <- bio3d::read.pdb(demostr)
      tryCatch({bonds <- jsonlite::read_json(gsub(".pdb", ".json" ,demostr))$bonds}, 
               error = function(e) {warning("no bonds JSON file", call. = FALSE)})
    } else if(length(unlist(contents)) > 0) {
      content_string <- unlist(strsplit(contents, ","))[2]
      decoded <- jsonlite::base64_dec(content_string)
      pdbData <- get_pdb(rawToChar(decoded))
    } else return('demostr and contents are none')

    atom_data <- pdbData$atom
    atoms <- lapply(1:dim(atom_data)[1],
                    function(i) {
                      atom <- atom_data[i, ]
                      
                      if(atom$type %in% c("ATOM", "HETATM")) {
                        list(
                          name = atom$elety,
                          chain = atom$chain,
                          positions = list(
                            atom$x,
                            atom$y,
                            atom$z
                          ),
                          residue_index = atom$resno,
                          element = atom$elesy ,
                          residue_name = paste0(atom$resid, atom$resno),
                          serial = i - 1
                        )
                      } else NULL
                    })
    # remove NULL
    atoms <- Filter(Negate(is.null), atoms)
    
    modelData <- list(atoms = atoms, 
                      bonds = bonds)
    styles <- create_style(mol_style = mol_style,
                           color_style = color_style,
                           custom_colors = custom_colors, 
                           atoms = atoms,
                           types = atom_data$type)
    
    dashbioMolecule3dViewer(
      id = 'mol-3d',
      selectionType = 'atom',
      modelData = modelData,
      styles = styles,
      selectedAtomIds = list(),
      backgroundOpacity = '0',
      atomLabelsShown = FALSE
    )
  }
)

app$callback(
  output = list(id = 'mol3d-selection-output', property = 'children'),
  params = list(
    input(id = 'mol-3d', property = 'selectedAtomIds'),
    input(id = 'mol-3d', property = 'modelData')),
  function(selected_atom_ids, model_data) {

    residue_summary <- if(length(selected_atom_ids) > 0) {
      lapply(selected_atom_ids,
             function(id) {
               res_info <- model_data$atoms[[id]]
               htmlP(
                 children = jsonlite::toJSON(
                   x = data.frame("residue" =  res_info$residue_name,
                                  "atom" = res_info$name, 
                                  "chain" = res_info$chain, 
                                  "xyz" = paste(unlist(res_info$positions), collapse = " ")),
                   pretty = TRUE
                 )
               )
             }
      )
    } else "No atoms have been selected. Click on an atom to select it."
    
    htmlDiv(residue_summary)
  }
)

#### dash app run ####
app$run_server(host="0.0.0.0", port=8050)

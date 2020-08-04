% Auto-generated: do not edit by hand
\name{DashNglDashNgl}

\alias{DashNglDashNgl}

\title{DashNgl component}

\description{
The NglMoleculeViewer is used to render schematic diagrams of biomolecules in ribbon-structure representations. Read more about the used WebGL protein viewer here: https://github.com/arose/ngl
}

\usage{
DashNglDashNgl(id=NULL, viewportStyle=NULL, stageParameters=NULL,
imageParameters=NULL, downloadImage=NULL, pdbString=NULL,
data=NULL, molStyles=NULL)
}

\arguments{
\item{id}{Character. The ID of this component, used to identify dash components in callbacks.
The ID needs to be unique across all of the components in an app.}

\item{viewportStyle}{Lists containing elements 'width', 'height'.
those elements have the following types:
  - width (character; optional)
  - height (character; optional). The height and the width (in px) of the container
in which the molecules will be displayed.
Default: width:1000px / height:500px
It should be in JSON format.}

\item{stageParameters}{Lists containing elements 'quality', 'backgroundcolor', 'cameratype'.
those elements have the following types:
  - quality (character; optional)
  - backgroundcolor (character; optional)
  - cameratype (character; optional). Parameters (in JSON format) for the stage object of ngl.
Currently implemented are the quality of the visualisation
and the background colorFor a full list see:
http://nglviewer.org/ngl/api/file/src/stage/stage.js.html}

\item{imageParameters}{Lists containing elements 'antialias', 'transparent', 'trim', 'defaultfilename'.
those elements have the following types:
  - antialias (logical; optional)
  - transparent (logical; optional)
  - trim (logical; optional)
  - defaultfilename (character; optional). Parameters (in JSON format) for exporting the image}

\item{downloadImage}{Logical. flag if download image was pressed}

\item{pdbString}{Character. Variable which defines how many molecules should be shown and/or which chain
The following format needs to be used:
pdbID1.chain:start-end@atom1,atom2_pdbID2.chain:start-end
. indicates that only one chain should be shown
: indicates that a specific amino acid range should be shown (e.g. 1-50)
@ indicates that chosen atoms should be highlighted (e.g. @50,100,150)
 _ indicates that more than one protein should be shown}

\item{data}{List of lists containing elements 'filename', 'ext', 'selectedvalue', 'chain', 'aarange', 'color', 'chosen', 'config', 'resetview', 'uploaded'.
those elements have the following types:
  - filename (character; required)
  - ext (character; required)
  - selectedvalue (character; required)
  - chain (character; required)
  - aarange (character; required)
  - color (character; required)
  - chosen (optional): . chosen has the following type: lists containing elements 'residues', 'atoms'.
those elements have the following types:
  - residues (character; required)
  - atoms (character; required)
  - config (optional): . config has the following type: lists containing elements 'type', 'input'.
those elements have the following types:
  - type (character; required)
  - input (character; required)
  - resetview (logical; required)
  - uploaded (logical; required)s. The data (in JSON format) that will be used to display the molecule
filename: name of the used pdb/cif file
ext: file extensions (pdb or cif)
selectedValue: pdbString
molStyles: selected molecule representation (cartoon, stick, sphere)
chain: selected chain
color: chain color
aaRange: selected residues range
config.input: content of the pdb file
config.type: format of config.input
resetView: bool if the view should be resettet
uploaded: bool if the file was uploaded}

\item{molStyles}{Lists containing elements 'representations', 'chosenatomscolor', 'chosenatomsradius', 'molspacingxaxis'.
those elements have the following types:
  - representations (list of characters; optional)
  - chosenatomscolor (character; required)
  - chosenatomsradius (numeric; required)
  - molspacingxaxis (numeric; required). The data (in JSON format) that will be used to style the displayed molecule
representations: one or multiple selected molecule representation
 - Possible molecule styles:
   'backbone,'ball+stick','cartoon', 'hyperball','licorice','line',
   'ribbon',''rope','spacefill','surface','trace','tube'
 - Possible additional representations:
   'axes','axes+box','helixorient','unitcell'
chosenAtomsColor: color of the 'ball+stick' representation of the chosen atoms}
}

\value{named list of JSON elements corresponding to React.js properties and their values}


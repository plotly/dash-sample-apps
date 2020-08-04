import React, {Component} from 'react';
import PropTypes from 'prop-types'; //exports a range of data validators
import {Stage, Selection, download } from 'ngl'; //https://github.com/arose/ngl/blob/master/src/stage/stage.ts
import { equals } from 'ramda';

/**
 * The NglMoleculeViewer is used to render schematic diagrams
 * of biomolecules in ribbon-structure representations.
 * Read more about the used WebGL protein viewer here:
 * https://github.com/arose/ngl
 */
export default class DashNgl extends Component {
  // constructor might not be needed anylonger:
  // https://hackernoon.com/the-constructor-is-dead-long-live-the-constructor-c10871bea599
  constructor (props) {
    super(props) // initiate the parent's constructor method and allows the component to inherit methods from its parent
    this.state = { stage: null, orientationMatrix: null, structuresList: [] } // initial values are set
    console.log(this.props)
    console.log(this.state)
  }

  // called after the component is rendered
  componentDidMount () {
    const { id, stageParameters, viewportStyle } = this.props
    const params = { ...stageParameters }
    const stage = new Stage(id, params)
    stage.setSize(viewportStyle.width, viewportStyle.height)
    this.setState({ stage })

    const orientationMatrix = stage.viewerControls.getOrientation();
    this.setState ( { orientationMatrix })

    console.log('component did mount')
  }

  componentWillUnmount() {
    console.log('component will unmount')
  }

  // triggered by any update of the DOM (e.g. new dropdown selection)
  shouldComponentUpdate (prevProps, nextProps) {
    const { stageParameters, data, downloadImage, molStyles } = this.props
    console.log(data)

    // check if data has changed
    if (data !== null && prevProps.data !== null) {
      console.log({ prevProps, nextProps })

      // wait for the first pdb selection after startup
      if (nextProps.data !== undefined) {
        console.log('first pdb selection')
        return true
      }

      // check if pdb selection has changed
      const oldSelection = prevProps.data[0].selectedValue
      const newSelection = data[0].selectedValue
      console.log({ oldSelection, newSelection })
      if (oldSelection !== newSelection) {
        console.log('pdb selection has changed')
        return true
      } 

      // check if view should be resetet
      const resetView = data[0].resetView
      if (oldSelection === newSelection && resetView === true){
        console.log('reset view')  
        return true
      }

      console.log(data)
      // check if structure has been uploaded
      if (data[0].uploaded === true){
        console.log('data has been uploaded')
        return true
      }
    } 

    // check if molStyles has been changed
    if (!equals(prevProps.molStyles,molStyles)){
        console.log("mol style has been changed")
        return true
      }

    // check for stage params changed
    if (!equals(prevProps.stageParameters, stageParameters)) {
            return true
        }

    // check if export image has been selected
    console.log(prevProps.downloadImage)
    console.log(downloadImage)
    if (prevProps.downloadImage !== downloadImage){
      return true
    }

    return false
  }

  // called only if shouldComponentUpdate evaluates to true
  componentDidUpdate () {
    console.log('updated')
    const { stageParameters, data, downloadImage} = this.props
    const { stage, structuresList } = this.state
    console.log(downloadImage)
    const newSelection = data[0].selectedValue

    stage.setParameters(stageParameters)
    
    console.log({ data, stageParameters })
    console.log(structuresList)

    if (
      downloadImage === undefined ||
      (downloadImage === false && newSelection !== 'placeholder') 
    ) { 
        // console.log(newSelection)
        stage.eachComponent(function (comp) {
          comp.removeAllRepresentations()
        })
        this.processDataFromBackend(data, structuresList)
    }

    if (downloadImage === true){
      this.generateImage()
    }
  }
  
  // helper function for generating an 
  generateImage() {
    const { imageParameters } = this.props
    const { stage } = this.state
    console.log('generate image')
    console.log(stage)
    
    stage.makeImage({
      factor: 1,
      antialias: imageParameters.antialias,
      trim: imageParameters.trim,
      transparent: imageParameters.transparent
    }).then(function( blob ){
        console.log(blob)
        download( blob, imageParameters.defaultFilename + ".png" );
    });
  //set downloadImage to false to prevent retriggering of the download handler
  this.props.setProps({'downloadImage': false})
  }

  //helper function for highlighting atoms
  highlightAtoms(args,sele,struc,chosenAtoms,chosenResidues){
    const repr = 'ball+stick'
    if (chosenAtoms!==''){
      args.sele = sele + ' and @' + chosenAtoms
      struc.addRepresentation(repr, args)
    }
    if (chosenResidues!== ''){
      args.sele = sele + '.CA and (' + chosenResidues.replace(/,/g,' or ') + ')'
      struc.addRepresentation(repr, args)
    }
  }

  //helper function for adding one or multiple molecular representations
  addMolStyle (struc, sele, chosen, color) {
    const { molStyles } = this.props
    const temp = molStyles.representations
    let reprs = [...temp]
    console.log(reprs)
    let chosenAtoms = chosen.atoms
    let chosenResidues = chosen.residues
    console.log(chosen)
    let args = {
      sele: sele,
      showBox: reprs.includes('axes+box')
    }
    if (sele !== ':'){
      args.color = color
    }

    console.log(reprs)
    reprs.forEach(e => {
      let repr = e
      if (repr === 'axes+box') {
        // This is not a ngl provided moleculuar representation
        // but a combination of repr: 'axes' and showBox = true 
        repr = 'axes'
      }
      struc.addRepresentation(repr, args)
    })

    // check if atoms should be selected
    if (chosenAtoms !=='' || chosenResidues !==''){
      args.radius = molStyles.chosenAtomsRadius
      args.color = molStyles.chosenAtomsColor
      this.highlightAtoms (args,sele,struc,chosenAtoms,chosenResidues)
    }
  }

  // helper functions which styles the output of loadStructure/loadData
  showStructure (stageObj, chain, aaRange, chosen, color, xOffset) {
    const { stage, orientationMatrix } = this.state
    let sele = ':'

    console.log(chosen)
    console.log("orientation Matrix")
    console.log(orientationMatrix)
    stage.viewerControls.orient(orientationMatrix);

    if (chain === 'ALL'){
      this.addMolStyle(stageObj,sele, chosen, color)
    } else {
      sele += chain
      if (aaRange !== 'ALL') {
        sele += '/0 and ' + aaRange
        console.log (sele)
      }

      const selection = new Selection(sele)
      const structure = stageObj.structure.getView(selection)
      const pa = structure.getPrincipalAxes()
      const struc = stage.addComponentFromObject(structure)
      const strucCenter = struc.getCenter()

      struc.setRotation(pa.getRotationQuaternion())
      struc.setPosition(
        [(0-strucCenter.x)-xOffset,
          0-strucCenter.y,
          0-strucCenter.z]
          )
      this.addMolStyle(struc, sele, chosen, color)
    } 
    
    //stage.animationControls.moveComponent(stageObj, stageObj.getCenter(), 1000)
    //const center = stage.getCenter()
    //stage.animationControls.zoom(newZoom, duration)
    //stage.animationControls.zoomMove(center, newZoom, duration)

    stage.autoView()
  }

  // loadStructure (stage, filename, molStyles, chain, range, color, xOffset) {
  //   console.log('load from browser')
  //   // console.log(filename)
  //   const stageObj = stage.getComponentsByName(filename).list[0]
  //   this.showStructure(stageObj, molStyles, chain, range, color, xOffset, stage)
  // }
  
  // If not load the structure from the backend
  processDataFromBackend (data, structuresList) {
    const { molStyles } = this.props
    const { stage } = this.state

    console.log('processDataFromBackend')
    
    // loop over list of structures:
    for (var i = 0; i < data.length; i++) {
      const filename = data[i].filename
      const xOffset = i * molStyles.molSpacingXaxis
      
      // check if already loaded
      if (structuresList.includes(filename)) {
        // If user has selected structure already just add the new Representation
        console.log('load from browser')
        this.showStructure(
          stage.getComponentsByName(filename).list[0],
          data[i].chain,
          data[i].aaRange,
          data[i].chosen,
          data[i].color,
          xOffset)
      } else { // load from backend
        this.loadData(data[i], xOffset)
      }
    }
    // const center = stage.getCenter()
    // const newZoom = -500
    // const duration = 1000
    // stage.animationControls.zoomMove(center, newZoom, duration)

    // let orientationMatrix = stage.viewerControls.getOrientation();
    // console.log(orientationMatrix)
    // stage.viewerControls.orient(orientationMatrix);

    // autoView = stage.animationControls.zoomMove(center, zoom, 0)
    // stage.animationControls.zoomMove(center, zoom, 0)

    // const center = stage.getCenter()
    // console.log(center)
    // let zoom = stage.getZoom()
    // zoom = zoom - 500
    // console.log(zoom)

    // // stage.autoView()
    // // change zoom depending on sub units
    // const zoom3 = -500
    // stage.animationControls.zoomMove(center, zoom3, 1000)
    // // stage.autoView(center,zoom,1000)

    // const center2 = stage.getCenter()
    // console.log(center2)
    // const zoom2 = stage.getZoom()
    // // zoom = zoom + 5
    // console.log(zoom2)

    // stage.viewer.zoom(distance, set)
    // let orientationMatrix2 = stage.viewerControls.getOrientation()
    // console.log(orientationMatrix2)
  }

  loadData (data, xOffset) {
    console.log('load from backend')
    const {stage} = this.state
    const stringBlob = new Blob([data.config.input], { type: data.config.type })
    
    stage.loadFile(stringBlob, { ext: data.ext, defaultRepresentation: false }).then(stageObj => {
      stageObj.name = data.filename
      this.showStructure(
        stageObj,
        data.chain,
        data.aaRange,
        data.chosen,
        data.color,
        xOffset
      )

      this.setState(state => ({
        structuresList: state.structuresList.concat([
          data.filename])
      }))
      console.log(this.state)
    })
  }

  render () {
    const { id } = this.props
    return (
      <div id={id}/>
    )
  }
}

const defaultViewportStyle = {
  width: '500x',
  height: '500px'
}

const defaultStageParameters = {
  quality: 'medium',
  backgroundColor: 'white',
  cameraType: 'perspective',
}

const defaultImageParameters = {
  antialias: true,
  transparent: true,
  trim: true,
  defaultFilename: 'dashNGL_output'
}

const defaultData = [{
  filename: 'placeholder',
  ext: '',
  selectedValue: 'placeholder',
  chain: 'ALL',
  aaRange: 'ALL',
  color: 'red',
  chosen: {
    'chosenAtoms':'',
    'chosenResidues':''
  },
  config: {
    type: 'text/plain',
    input: ''
  },
  uploaded: true,
  resetView: false,
}]

console.log(defaultData)

DashNgl.defaultProps = {
  data: defaultData,
  viewportStyle: defaultViewportStyle,
  stageParameters: defaultStageParameters,
  imageParameters: defaultImageParameters,
  downloadImage: false,
  molStyles:{
    representations:['cartoon','axes+box'],
    chosenAtomsColor:'#808080',
    chosenAtomsRadius: 1,
    molSpacingXaxis: 100,
  }
}

DashNgl.propTypes = {
  /**
   * The ID of this component, used to identify dash components in callbacks.
   * The ID needs to be unique across all of the components in an app.
   */
  id: PropTypes.string,

  /**
   * The height and the width (in px) of the container
   * in which the molecules will be displayed.
   * Default: width:1000px / height:500px
   * It should be in JSON format.
  */
  viewportStyle: PropTypes.exact({
    width: PropTypes.string,
    height: PropTypes.string
  }),

  /**
   * Parameters (in JSON format) for the stage object of ngl.
   * Currently implemented are the quality of the visualisation
   * and the background colorFor a full list see:
   * http://nglviewer.org/ngl/api/file/src/stage/stage.js.html
   */
  stageParameters: PropTypes.exact({
    quality: PropTypes.string,
    backgroundColor: PropTypes.string,
    cameraType: PropTypes.string,
  }),

   /**
   * Parameters (in JSON format) for exporting the image
   */
  imageParameters: PropTypes.exact({
    antialias: PropTypes.bool,
    transparent: PropTypes.bool,
    trim: PropTypes.bool,
    defaultFilename: PropTypes.string
  }),

  /**
   * flag if download image was pressed
   */
  downloadImage: PropTypes.bool,

  /**
   * Variable which defines how many molecules should be shown and/or which chain
   * The following format needs to be used:
   * pdbID1.chain:start-end@atom1,atom2_pdbID2.chain:start-end
   * . indicates that only one chain should be shown
   * : indicates that a specific amino acid range should be shown (e.g. 1-50)
   * @ indicates that chosen atoms should be highlighted (e.g. @50,100,150)
   *  _ indicates that more than one protein should be shown
   */
  pdbString: PropTypes.string,

  /**
   * The data (in JSON format) that will be used to display the molecule
   * filename: name of the used pdb/cif file
   * ext: file extensions (pdb or cif)
   * selectedValue: pdbString
   * molStyles: selected molecule representation (cartoon, stick, sphere)
   * chain: selected chain
   * color: chain color
   * aaRange: selected residues range
   * config.input: content of the pdb file
   * config.type: format of config.input
   * resetView: bool if the view should be resettet
   * uploaded: bool if the file was uploaded
   */
  data: PropTypes.arrayOf(
    PropTypes.exact({
      filename: PropTypes.string.isRequired,
      ext: PropTypes.string.isRequired,
      selectedValue: PropTypes.string.isRequired,
      chain: PropTypes.string.isRequired,
      aaRange: PropTypes.string.isRequired,
      color: PropTypes.string.isRequired,
      chosen:PropTypes.exact({
        residues:PropTypes.string.isRequired,
        atoms:PropTypes.string.isRequired,
      }),
      config: PropTypes.exact({
        type: PropTypes.string.isRequired,
        input: PropTypes.string.isRequired
      }),
      resetView: PropTypes.bool.isRequired, 
      uploaded: PropTypes.bool.isRequired,
    })
  ),
  /**
   * The data (in JSON format) that will be used to style the displayed molecule
   * representations: one or multiple selected molecule representation
   *  - Possible molecule styles:
   *    'backbone,'ball+stick','cartoon', 'hyperball','licorice','line',
   *    'ribbon',''rope','spacefill','surface','trace','tube'
   *  - Possible additional representations:
   *    'axes','axes+box','helixorient','unitcell'
   * chosenAtomsColor: color of the 'ball+stick' representation of the chosen atoms
   */
  molStyles:
    PropTypes.exact({ 
      representations: PropTypes.arrayOf(PropTypes.string),
      chosenAtomsColor: PropTypes.string.isRequired,
      chosenAtomsRadius: PropTypes.number.isRequired,
      molSpacingXaxis: PropTypes.number.isRequired
    })
}

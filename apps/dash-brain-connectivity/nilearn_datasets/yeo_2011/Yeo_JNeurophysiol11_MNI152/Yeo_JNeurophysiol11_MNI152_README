Background
==========
Resting state fMRI data from 1000 subjects were registered using surface-based alignment. A clustering approach was employed to identify and replicate 7 and 17 networks of functionally coupled regions across the cerebral cortex. The results revealed local networks confined to sensory and motor cortices as well as distributed networks of association regions that form interdigitated circuits. Within the sensory and motor cortices, functional connectivity followed topographic representations across adjacent areas. In association cortex, the connectivity patterns often showed abrupt transitions between network boundaries, forming largely parallel circuits.


Files in Folder (Excluding README)
==================================
1) FSL_MNI152_FreeSurferConformed_1mm.nii.gz
2) Yeo2011_7Networks_MNI152_FreeSurferConformed1mm.nii.gz
3) Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz
4) Yeo2011_7Networks_ColorLUT.txt
5) Yeo2011_17Networks_MNI152_FreeSurferConformed1mm.nii.gz
6) Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz
7) Yeo2011_17Networks_ColorLUT.txt


Descriptions of files
=====================
1) FSL_MNI152_FreeSurferConformed_1mm.nii.gz is the FSL MNI152 1mm template interpolated and intensity normalized into a 256 x 256 x 256 1mm-isotropic volume (obtained by putting the FSL MNI152 1mm template through recon-all using FreeSurfer 4.5.0)

2) Yeo2011_7Networks_MNI152_FreeSurferConformed1mm.nii.gz is a volume consisting of 7 cortical networks projected into MNI152 space.The cortical ribbon is defined by putting the FSL MNI152 1mm template through recon-all using FreeSurfer 4.5.0. The slices of this volume is shown in Yeo et al., 2011. 

3) Yeo2011_7Networks_MNI152_FreeSurferConformed1mm.nii.gz is a volume consisting of 7 cortical networks projected into MNI152 space. The cortical ribbon is defined in a more liberal fashion than in (2). More specifically, the cortical mask is obtained by nonlinear warping 1000 subjects (from Yeo et al. 2011, Buckner et al. 2011) into MNI152 space via the FreeSurfer recon-all pipeline. An initial mask is first obtained where a voxel is decided to be a cortical voxel if the cortex of more than 150 subjects (out of 1000 subjects) were mapped to the voxel or if the voxel is labeled as part of the cortical ribbon from (2). This cortical mask is then smoothed, and small holes and islands in the masks are removed by simple morphological operations. 

4) Yeo2011_7Networks_ColorLUT.txt is a FreeSurfer readable text file specifying how the 7 networks are named, numbered and colored in Yeo et al. 2011:

  0            NONE   0   0   0   0
  1     7Networks_1 120  18 134   0
  2     7Networks_2  70 130 180   0
  3     7Networks_3   0 118  14   0
  4     7Networks_4 196  58 250   0
  5     7Networks_5 220 248 164   0
  6     7Networks_6 230 148  34   0
  7     7Networks_7 205  62  78   0

In particular the networks are numbered from 7Networks_1 to 7Networks_7. The first column of the text file specifies the value of voxels in the nifty values corresponding to the particular network. The second column of the text file specifies the named of the networks. For example, from the text file, voxels whose values = 3 corresponds to the network 7Networks_3.Columns 3 to 5 corresponds to the R, G, B values (ranges from 0 to 255) of the networks. Last column is all zeros (FreeSurfer's default).

5) Yeo2011_17Networks_MNI152_FreeSurferConformed1mm.nii.gz is a volume consisting of 17 cortical networks projected into MNI152 space. The cortical ribbon is defined in the same fashion as (2), i.e., by putting the FSL MNI152 1mm template through recon-all using FreeSurfer 4.5.0. 

6) Yeo2011_17Networks_MNI152_FreeSurferConformed1mm.nii.gz is a volume consisting of 17 cortical networks projected into MNI152 space. The cortical ribbon is defined in a more liberal fashion than in (2) and in the same way as (3). 

7) Yeo2011_17Networks_ColorLUT.txt is a FreeSurfer readable text file specifying how the 17 networks are named, numbered and colored in Yeo et al. 2011:

  0            NONE   0   0   0   0
  1    17Networks_1 120  18 134   0
  2    17Networks_2 255   0   0   0
  3    17Networks_3  70 130 180   0
  4    17Networks_4  42 204 164   0
  5    17Networks_5  74 155  60   0
  6    17Networks_6   0 118  14   0
  7    17Networks_7 196  58 250   0
  8    17Networks_8 255 152 213   0
  9    17Networks_9 220 248 164   0
 10   17Networks_10 122 135  50   0
 11   17Networks_11 119 140 176   0
 12   17Networks_12 230 148  34   0
 13   17Networks_13 135  50  74   0
 14   17Networks_14  12  48 255   0
 15   17Networks_15   0   0 130   0
 16   17Networks_16 255 255   0   0
 17   17Networks_17 205  62  78   0


Example Usage
=============
1) Except for the colortables, all the volumes are nifty volumes which can be read using any software like freeview (FreeSurfer), fslview (FSL), etc. 

2) To overlay the 7-network volume on the FSL MNI152 1mm template in freeview (assuming FreeSurfer is already installed) with the appropriate colors, type the following assuming the working directory is in the same directory as this README:

freeview -v FSL_MNI152_FreeSurferConformed_1mm.nii.gz Yeo2011_7Networks_MNI152_FreeSurferConformed1mm.nii.gz:colormap=lut:lut=Yeo2011_7Networks_ColorLUT.txt


Other Downloads
===============
1) The surface parcellations in Caret space can be viewed using Webcaret here: http://sumsdb.wustl.edu:8081/sums/directory.do?id=8286317

2) The surface parcellations in FreeSurfer space can also be downloaded here: http://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation_Yeo2011

3) Seed-based fcMRI Movies can be downloaded here: http://www.youtube.com/YeoKrienen


References
==========
Yeo BT, Krienen FM, Sepulcre J, Sabuncu MR, Lashkari D, Hollinshead M, Roffman JL, Smoller JW, Zollei L., Polimeni JR, Fischl B, Liu H, Buckner RL (2011) The organization of the human cerebral cortex estimated by intrinsic functional connectivity. J. Neurophysiol. In Press.

The paper can be downloaded from PubMed here: http://www.ncbi.nlm.nih.gov/pubmed/21653723


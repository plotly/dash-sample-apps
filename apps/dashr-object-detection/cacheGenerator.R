library(R.cache)
library(data.table)
library(compiler)
library(glue)
DEBUG <- T
FRAMERATE <- 24.0
f_read <- cmpfun(fread)
data_table <- cmpfun(data.table)
gum <- cmpfun(glue)
load_data <- cmpfun(function(path){
  #  # Load data about a specific footage (given by the path). It returns a listionary of useful variables such as
  #  #  the dataframe containing all the detection and bounds localization, the number of classes inside that footage,
  #  #  the matrix of all the classes in string, the given class with padding, and the root of the number of classes,
  #  #  rounded.
  #  #  Load DT containing all the processed object detections inside the video
  Info.DT <- f_read(path)
  
  #cls.vec <- the list of detected object classes for the .csv file
  
  cls.vec <- Info.DT[, class_str]
  n.cls <- length(cls.vec)
  title.vec <- unique(cls.vec)
  
  # Gets the smallest value needed to add to the end of the classes list to get a square matrix
  #   # n := argmin_{n : int} (n \geq n.c) \land (\sqrt{n} : int)
  #   # to find n we do
  
  r.rnd <- ceiling(sqrt(n.cls))
  
  pad.ent<- as.integer(r.rnd^2 - n.cls)
  pad.seq <- numeric(pad.ent)
  #Pad the class.vec
  cls.wpad <- c(cls.vec, pad.seq)
  cls.mtx <- apply(matrix(cls.wpad, nrow=r.rnd, ncol=r.rnd, byrow=T),2,rev)
  dat_list = list(INFO_DT = Info.DT, N_CLS=n.cls, MTX_CLS = cls.mtx, WPAD_CLS = cls.wpad, RT_RND = r.rnd, TITL=title.vec)
  if(DEBUG)
  {
    print(gum('{path} loaded.'))
  }
  return(dat_list)
})

dat.list <- list( james_bond = load_data("data/james_bond_object_data.csv"), 
                  
                  zebra =  load_data("data/Zebra_object_data.csv"), 
                  
                  car_show_drone = load_data("data/CarShowDrone_object_data.csv"),  
                  
                  car_footage = load_data("data/CarFootage_object_data.csv"),
                  
                  DroneCanalFestival = load_data("data/DroneCanalFestivalDetectionData.csv"),
                  
                  DroneCarFestival2 =  load_data("data/DroneCarFestival2DetectionData.csv"),
                  
                  FarmDrone = load_data("data/FarmDroneDetectionData.csv"),
                  
                  ManCCTV = load_data("data/ManCCTVDetectionData.csv"),
                  
                  RestaurantHoldup =load_data("data/RestaurantHoldupDetectionData.csv"))

url.list <-  list(regular=data_table(james_bond = 'https://www.youtube.com/watch?v=g9S5GndUhko', 
                                     
                                     zebra =  'https://www.youtube.com/watch?v=TVvtD3AVt10', 
                                     
                                     car_show_drone = 'https://www.youtube.com/watch?v=gPtn6hD7o8g',  
                                     
                                     car_footage = 'https://www.youtube.com/watch?v=qX3bDxHuq6I',
                                     
                                     DroneCanalFestival = 'https://youtu.be/0oucTt2OW7M',
                                     
                                     DroneCarFestival2 =  'https://youtu.be/vhJ7MHsJvwY',
                                     
                                     FarmDrone = 'https://youtu.be/aXfKuaP8v_A',
                                     
                                     ManCCTV =  'https://youtu.be/BYZORBIxgbc',
                                     
                                     RestaurantHoldup = 'https://youtu.be/WDin4qqgpac'), 
                  
                  bounding_box=data_table( james_bond = 'https://www.youtube.com/watch?v=g9S5GndUhko', 
                                           
                                           zebra =  'https://www.youtube.com/watch?v=G2pbZgyWQ5E', 
                                           
                                           car_show_drone = 'https://www.youtube.com/watch?v=9F5FdcVmLOY',  
                                           
                                           car_footage = 'https://www.youtube.com/watch?v=EhnNosq1Lrc',
                                           
                                           DroneCanalFestival = 'https://youtu.be/6ZZmsnwk2HQ',
                                           
                                           DroneCarFestival2 =  'https://youtu.be/2Gr4RQ-JHIs',
                                           
                                           FarmDrone = 'https://youtu.be/pvvW5yZlpyc',
                                           
                                           ManCCTV = 'https://youtu.be/1oMrHLrtOZw',
                                           
                                           RestaurantHoldup ='https://youtu.be/HOIKOwixYEY')
) 

saveCache(dat.list, pathname ="assets/datList.Rcache")
saveCache(url.list, pathname ="assets/urlList.Rcache")


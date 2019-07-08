
  COUNTIES1 = c('Albany','Allegany','bronx', 'Broome','Cattaraugus','Cayuga',
    'Chautauqua','Chemung','Chenango','Clinton','Columbia','Cortland','Delaware','Dutchess',
    'Erie','Essex','Franklin','Fulton','Genesee','Greene','Hamilton','Herkimer','Jefferson','Kings','Lewis',
    'Livingston','Madison','Monroe','Montgomery','Nassau', 'New York','Niagara','Oneida','Onondaga','Ontario',
    'Orange','Orleans','Oswego','Otsego','Putnam','Queens','Rensselaer','Richmond',
    'Rockland','St. Lawrence','Saratoga','Schenectady','Schoharie','Schuyler',
    'Seneca','Steuben','Suffolk','Sullivan','Tioga',
    'Tompkins','Ulster','Warren','Washington','Wayne','Westchester','Wyoming','Yates')
  
  
  WELL_STATUSES1 <- c('AC','AR', 'CA', 'DC', 'DD', 'DG', 'EX', 'IN', 'NR', 'PA', 'PI', 'PB', 'PM', 'RE', 'RW', 'SI', 'TA',
                     'TR', 'UN', 'UL', 'UM', 'VP')
  WELL_STATUSES2 <- c(
    'Active','Application Received to Drill/Plug/Convert','Cancelled','Drilling Completed','Drilled Deeper','Drilling in Progress',
    'Expired Permit', 'Inactive', 'Not Reported on AWR', 'Plugged and Abandoned', 'Permit Issued', 'Plugged Back', 
    'Plugged Back Multilateral', 'Refunded Fee','Released - Water Well','Shut-In','Temporarily Abandoned','Transferred Permit','Unknown','Unknown Located',
    'Unknown Not Found','Voided Permit')
  
  
  WELL_TYPES1 <- c('BR', 'Confidential', 'DH', 'DS', 'DW', 'GD', 'GE', 'GW', 
                   'IG', 'IW', 'LP', 'MB', 'MM', 'MS', 'NL', 'OB', 'OD', 'OE', 'OW', 'SG', 'ST', 'TH', 'UN')
  
  WELL_TYPES2 = c('Brine','Confidential','Dry Hole','Disposal','Dry Wildcat','Gas Development','Gas Extension','Gas Wildcat',
                    'Gas Injection','Oil Injection','Liquefied Petroleum Gas Storage','Monitoring Brine','Monitoring Miscellaneous',
                    'Monitoring Storage','Not Listed','Observation Well','Oil Development','Oil Extension','Oil Wildcat',
                    'Stratigraphic','Storage','Geothermal','Unknown')
  
  WELL_COLORS1 <- c('GD', 'GE',  'GW', 'IG', 'OD', 'OE', 'OW', 'ST', 'BR', 'MB', 'IW', 'LP', 'MS', 'Confidential', 
                    'DH', 'DS', 'DW', 'MM', 'NL', 'OB', 'SG', 'TH', 'UN')
   
  
  WELL_COLORS2 = c('#FFEDA0','#FA9FB5','#A1D99B','#67BD65','#BFD3E6','#B3DE69','#FDBF6F','#FC9272','#D0D1E6','#ABD9E9',
  '#3690C0','#F87A72','#CA6BCC','#DD3497','#4EB3D3','#FFFF33','#FB9A99','#A6D853','#D4B9DA','#AEB0B8','#CCCCCC','#EAE5D9','#C29A84')
  
  
  COUNTIES <- vector(mode="list", length= 62)
  names(COUNTIES) <- seq(1,123,2)    
  for (i in 1:62) {
    COUNTIES[[i]] <- COUNTIES1[i]
  }
  
  #Print(COUNTIES)
  
  WELL_STATUSES <- vector(mode="list", length= length(WELL_STATUSES1))
  names(WELL_STATUSES) <- WELL_STATUSES1    
  for (i in 1:length(WELL_STATUSES1)) {
    WELL_STATUSES[[i]] <- WELL_STATUSES2[i]
  }
  
  #Print(WELL_STATUSES)
  
  WELL_TYPES <- vector(mode="list", length= length(WELL_TYPES1))
  names(WELL_TYPES) <- WELL_TYPES1    
  for (i in 1:length(WELL_TYPES1)) {
    WELL_TYPES[[i]] <- WELL_TYPES2[i]
  }
  
  #Print(WELL_STATUSES)
  
  WELL_TYPES <- vector(mode="list", length= length(WELL_TYPES1))
  names(WELL_TYPES) <- WELL_TYPES1    
  for (i in 1:length(WELL_TYPES1)) {
    WELL_TYPES[[i]] <- WELL_TYPES2[i]
  }
  
  WELL_COLORS <- vector(mode="list", length= length(WELL_TYPES1))
  names(WELL_COLORS) <- WELL_COLORS1    
  for (i in 1:length(WELL_COLORS1)) {
    WELL_COLORS[[i]] <- WELL_COLORS2[i]
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

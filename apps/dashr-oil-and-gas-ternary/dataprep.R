df <- read.csv("data/Well_ID-Table_1.csv")
df[c(7,8,9,12,13,14)] <- lapply(df[c(7,8,9,12,13,14)], function(x) trimws(as.character(x)))
df <- df[df['result'] == 'GAS',]
df <- df[df['tdfm'] != '000',] 
form_id <- unique(df['tdfm'])
form_id <- form_id$tdfm
form_id <- lapply(form_id, function(x) trimws(as.character(x)))


df_fm <- read.csv("data/FmCodes-Table_1.csv")
df_fm <- as.data.frame(lapply(df_fm, trimws))

df['fm_name'] <- df['tdfm']
df_fm <- df_fm[-c(1),] 

formation_name <- list()
i <- 1
for (formation in form_id) {
  formation_name[[i]] <- df_fm[df_fm['fm_code'] == formation,'fm_name']
  for (name in as.character(df$fm_name)) {
    if(name == formation){
      df[df$fm_name == name,'fm_name'] = trimws(as.character(formation_name[[i]]))
    }else{
    }
  }
  i <- i+1
}

ternary_ax_title <- list('Quartz', 'Carbonate', 'Clay')
rand_composition <- matrix(sample(1:100,1137,replace = TRUE),ncol = 379)

df['Quartz'] <- rand_composition[1,]
df['Carbonate'] <- rand_composition[2,]
df['Clay'] <- rand_composition[3,]

write.csv(df, 'data/test_compositionR.csv')

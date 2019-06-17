library(sas7bdat)

setwd('~/Dropbox (HMS)/RagGroup Team Folder/TEDDY/22417/TEDDY_DATA/')

files=list.files()
files=files[grep('sas7bdat',files)]
files=files[-grep('diet_record_dataset',files)]

for(file in files){
  d=read.sas7bdat(file)
  print(file)
  file=sub('\\.sas7bdat$', '', file)
  write.csv(d,paste(file,'.csv',sep=''))
}


files=list.files()

for(file in files){
  d=readRDS(file)
  print(file)
  file=sub('\\.Rda$', '', file)
  write.csv(d,paste(file,'.csv',sep=''))
}
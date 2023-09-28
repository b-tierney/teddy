#clean TEDDY output

library(tidyverse)

files = list.files()
metadata=files[grep('metadata',files)]
files=files[-grep('metadata',files)]

metadata_load = readRDS(metadata)
metadata_load = metadata_load %>% select(-X)
saveRDS(metadata_load,metadata)

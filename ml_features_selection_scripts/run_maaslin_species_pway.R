library(Maaslin2)
library(data.table)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

abundance_file=args[1]
metadata_file=args[2]
HLA_list=args[3]
condition_list = args[4]
output_dir = args[5]
metadata_suffix = args[6]

if(metadata_suffix == "species") {
  d = fread(abundance_file,sep="\t",header=TRUE,data.table=FALSE)
  rownames(d) = d[,1]
  d = d[,-1]
} else {
  d = fread(abundance_file,sep="\t",header=TRUE,data.table=FALSE)
  rownames(d) = d[,1]
  d = d[,-1]
  colnames(d) = gsub("_human_free_Abundance","",colnames(d))
}


# read in metadata
metadata = read.csv(metadata_file)
metadata = metadata[,-1]
# get data in both microbiome and metadata

run_maaslin = function(d,metadata,condition,HLA,output_dir) {
  if(condition == "MIAA") {
    controls = metadata %>% dplyr::filter(is.na(age_first_MIAA)) %>% mutate(condition = 0)
    cases = metadata %>% dplyr::filter(age_at_collection<age_first_MIAA,!is.na(age_first_MIAA)) %>% mutate(condition = 1)
  } else if (condition == "GAD") {
    controls = metadata %>% dplyr::filter(is.na(age_first_GAD)) %>% mutate(condition = 0)
    cases = metadata %>% dplyr::filter(age_at_collection<age_first_GAD,!is.na(age_first_GAD)) %>% mutate(condition = 1)
  } else if (condition == "IA2A"){
    controls = metadata %>% dplyr::filter(is.na(age_first_IA2A)) %>% mutate(condition = 0)
    cases = metadata %>% dplyr::filter(age_at_collection<age_first_IA2A,!is.na(age_first_IA2A)) %>% mutate(condition = 1)
  } else if (condition == "seroconverters"){
    controls = metadata %>% dplyr::filter(t1d_sero_control == 'control') %>% mutate(condition = 0)
    cases = metadata %>% dplyr::filter(age_at_collection<age_mult_persist,!is.na(age_mult_persist)) %>% mutate(condition = 1)
  } else if (condition == "triple_converters_vs_T1D") {
    controls = metadata %>% dplyr::filter(t1d_sero_control == 'seroconverted') %>% filter(three_persist_conf == TRUE) %>% mutate(condition = 0)
    cases = metadata %>% dplyr::filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
  } else if (condition == "T1D") {
    controls = metadata %>% dplyr::filter(t1d == FALSE) %>% mutate(condition = 0)
    cases = metadata %>% dplyr::filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
  } else if (condition == "serconverters_or_T1D") {
    controls = metadata %>% dplyr::filter(t1d_sero_control == 'control') %>% mutate(condition = 0)
    cases_T1D = metadata %>% dplyr::filter(t1d == TRUE, T1D_Outcome=='Before') %>% mutate(condition = 1)
    cases_sero = metadata %>% dplyr::filter(age_at_collection<age_mult_persist,!is.na(age_mult_persist)) %>% mutate(condition = 1)
    cases_run = unique(c(cases_T1D$Run,cases_sero$Run))
    cases = metadata[match(cases_run,metadata$Run),]  %>% mutate(condition = 1)
  }
  metadata = rbind(controls,cases)
  # get metadata specific to HLA

  if(HLA == "DR3_DR4_only") {
    # remove subjects that are not eligable. only 68 samples so shouldn't be a huge deal to remove
    metadata = metadata %>% filter(HLA_Category.x!= "Not*Eligible")
    metadata = metadata %>% filter(HLA_Category.x== "DR4*030X/0302*DR3*0501/0201")
  } else if (HLA == "DR4_DR4_only") {
    metadata = metadata %>% filter(HLA_Category.x!= "Not*Eligible")
    metadata = metadata %>% filter(HLA_Category.x== "DR4*030X/0302*DR4*030X/0302")
  } else if (HLA == "DR4_DR8_only") {
    metadata = metadata %>% filter(HLA_Category.x!= "Not*Eligible")
    metadata = metadata %>% filter(HLA_Category.x== "DR4*030X/0302*DR8*0401/0402")
  } else if (HLA == "DR3_DR3_only") {
    metadata = metadata %>% filter(HLA_Category.x!= "Not*Eligible")
    metadata = metadata %>% filter(HLA_Category.x== "DR3*0501/0201*DR3*0501/0201")
  } else if (HLA == "DR4_DR1_only") {
    metadata = metadata %>% filter(HLA_Category.x!= "Not*Eligible")
    metadata = metadata %>% filter(HLA_Category.x== "DR4*030X/0302*DR1*0101/0501")
  } else if (HLA == "DR4_DR13") {
    metadata = metadata %>% filter(HLA_Category.x!= "Not*Eligible")
    metadata = metadata %>% filter(HLA_Category.x== "DR4*030X/0302*DR13*0102/0604")
  }

  # get samples in both metadata and abundance
  metadata_abundance_samples = intersect(metadata$Run,colnames(d))

  d = d[,match(metadata_abundance_samples,colnames(d))]
  d = t(d)
  metadata = metadata[match(metadata_abundance_samples,metadata$Run),]
  rownames(metadata) = metadata$Run

  abundance_prefix = gsub(".tsv","",basename(abundance_file))

  output_folder = paste(output_dir,"/",condition,"_",HLA,"_",abundance_prefix,"_maalin_output",sep="")

  maaslin_output = Maaslin2(input_data=d,input_metadata=metadata,output=output_folder,min_abundance = 0.0,min_prevalence=0.3,max_significance=0.05,normalization="NONE",transform="LOG",analysis_method = "LM",fixed_effects  = c("condition", "age_at_collection"),random_effects=c("maskid"),correction = "BH",standardize = TRUE,cores = 1,plot_scatter=FALSE)
  return(maaslin_output)
}

HLA_list = strsplit(HLA_list,split=",")[[1]]
condition_list = strsplit(condition_list,split=",")[[1]]

for(HLA in HLA_list) {
  print(HLA)
  for(condition in condition_list) {
    print(condition)
    maaslin_output = run_maaslin(d,metadata,condition,HLA,output_dir)
  }
}
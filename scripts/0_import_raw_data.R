################################################################################
#############                                                      #############
#############     Assessing the causal effect of plumage traits    #############
#############     on reproduction using Mendelian randomization    #############        
#############                                                      #############
#############                    0. Data Import                    #############
#############                                                      #############
#############                   By: Zach Laubach                   #############
#############                 created: 9 Sept 2025                 #############
#############              last updated: 10 June 2026              #############
################################################################################


  ### PURPOSE: Load barn swallow trait, reproductive success, and genome meta
              # data from raw data files
  
  
  # Code Blocks
    # 1. Configure work space
    # 2. Import pheno. and repro. data
    # 3. Import genomic data
    # 4. Import relatedness data
    # 5. Export data 
  



###############################################################################
##############             1.  Configure work space              ##############
###############################################################################

  ### 1.1 Global options
    ## a) clear global environment
      rm(list = ls())

    ## b) prevent R from automatically reading character strings as factors
      options(stringsAsFactors = FALSE)

      
  ### 1.2 Install and load CRAN packages   
    ## a) Data Manipulation and Descriptive Stats Packages
      # load tidyverse package
        library ('tidyverse')
      
      # load vcfR package
        library ('vcfR')
 
      # load here package
        library ('here')

        
  ### 1.3 Get Version and Session Info
    R.Version()
    sessionInfo()
    
    # Developed in:   
    # R version 4.6.0 (2026-04-24)
    # Platform: x86_64-apple-darwin20
    # Running under: macOS Tahoe 26.5.1
    
  
  ### 1.4 Set working directory 
    setwd(here())
    
    
      
###############################################################################
##############        2. Import pheno. and repro. data           ##############
###############################################################################    
      
  ### 2.1 Import phenotype data file
    ## a)  Import AdultData2008-2024.csv as df
      pheno_reprod_df <- read_csv(here('raw_data/AdultData2008-2024.csv'))
    
  
  ### 2.2 Import paternity data
    ## a) Import ferts_alive_d12_2021-2023.csv as df
      paternity_df <- read_csv(here('raw_data/ferts_alive_d12_2021-2023.csv'))
   
    ## b) Import no_nestlings_analyzed.csv as df
      # a list of adult birds for which no nestlings were sampled and analyzed 
      no_nestlings_df <- read_csv(here('raw_data/no_nestlings_analyzed.csv'))
      
    ## c) Import unkwn_IDobs.csv as df
      # a list of adult birds whose social nest was not confirmed via ID obs
      unkwn_IDobs_df <- read_csv(here('raw_data/unkwn_IDobs.csv'))  
      
    ## d) Import true_zero_fert_males.csv as df
      # a list of males who sired no nestlings
      zero_fert_males_df <- read_csv(here('raw_data/true_zero_fert_males.csv'))  
      
     
    
###############################################################################
##############              3. Import genomic data               ##############
###############################################################################    
    
  ### 3.1 Import imputed snps for candidate loci data
    ## a) imputed snps for ventral color candidate loci data 
      snps_plum_df <- read_csv(here('raw_data/',
                                    'candidate_snps_vent_color.csv'))   
      
    ## b) imputed snps for ventral color candidate loci data 
      snps_TS_df <- read_csv(here('raw_data/',
                               'candidate_snps_ts.csv'))   
      
    ## c) left join snps_plum_df to snps_TS_df
      snps_df <- snps_plum_df  %>%
        left_join(snps_TS_df,
                  by = c('sample.id' = 'sample.id'), 
                  copy = F) 
  
        
  ### 3.2 Import genome meta data
    ## a) create a list of txt files that match the pattern
      coverage_list <- list.files(path = here('cov/'),
                                  pattern = '^coverage.*\\.txt$', 
                                  recursive = T, full.names = T)
      
    ## b) load files from list into a single data frame
      coverage_df <- do.call(rbind, lapply(coverage_list, 
                                           read.table))  
    

         
###############################################################################
##############            4. Import relatedness data             ##############
###############################################################################    
      
  ### 4.1 Import relatedness data
    ## a) 2021-2023 relatedness data
      relate_df <- read.table(here('raw_data/',
                                'CO-2021-2023.thin10k.GL-recalc.relate'), 
                               sep = '\t', header = TRUE)
 
      
  
###############################################################################
##############                  6. Export data                   ##############
###############################################################################
      
  ### 6.1 Export data to an RData file     
    ## a) Save and export raw repro. success data 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      save(file = here('data/1_pheno_reprod_df.RData'), 
           list = c('pheno_reprod_df', 'paternity_df', 'no_nestlings_df',
                    'unkwn_IDobs_df', 'zero_fert_males_df'))
      
    ## b) Save and export imputed snps for candidate loci
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      save(file = here('data/1_snps_df.RData'), 
           list = c('snps_df'))
      
    ## d) Save and export genomic meta data
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      save(file = here('data/1_genome_meta_df.RData'), 
           list = c('coverage_df'))
      
    ## e) Save and export raw relatedness data
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      save(file = here('data/2_relate_df.RData'), 
           list = c('relate_df'))
      
    
      
   
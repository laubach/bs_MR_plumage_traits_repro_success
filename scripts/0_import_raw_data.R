################################################################################
#############                                                      #############
#############     Assessing the causal effect of plumage traits    #############
#############     on reproduction using Mendelian randomization    #############        
#############                                                      #############
#############                    0. Data Import                    #############
#############                                                      #############
#############                   By: Zach Laubach                   #############
#############                 created: 9 Sept 2025                 #############
#############             last updated: 16 April 2026              #############
################################################################################


  ### PURPOSE: Load barn swallow trait, reproductive success, and genome meta
              # data from raw data files
  
  
  # Code Blocks
    # 1. Configure work space
    # 2. Import pheno. and repro. data
    # 3. Import genome meta data
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
      # load tidyverse packages
        library ('tidyverse')
 
      # load here packages
        library ('here')

        
  ### 1.3 Get Version and Session Info
    R.Version()
    sessionInfo()
    
    # Developed in:   
    # R version 4.4.2 (2024-10-31)
    # Platform: x86_64-apple-darwin20
    # Running under: macOS Sequoia 15.1.1
    
  
  ### 1.4 Set working directory 
    setwd(here())
    
    
      
###############################################################################
##############        2. Import pheno. and repro. data           ##############
###############################################################################    
      
  ### 2.1 Import sample data file
    ## a)  Import AdultData2008-2024.csv as df
      pheno_reprod_df <- read_csv(here('raw_data/AdultData2008-2024.csv'))
    
    
    
###############################################################################
##############            3. Import genome meta data             ##############
###############################################################################    
      
  ### 3.1  Import genome meta data
    ## a) Import CO_genomes_2021-2023_metadata.csv as df
      genome_meta_df <- read_csv(here('raw_data/',
                                      'CO_genomes_2021-2023_metadata.csv'))    

     
         
###############################################################################
##############            4. Import relatedness data             ##############
###############################################################################    
      
  ### 4.1 Import relatedness data
    ## a) 2021-2023 relatedness data
      relate_df <- read.table(here('raw_data/',
                                'CO-2021-2023.thin10k.GL-recalc.relate'), 
                               sep = '\t', header = TRUE)
    
  
    
###############################################################################
##############              5. Import genomic data               ##############
###############################################################################    
    
  ### 5.1 Import genomic data (vcfgz files)
    ## a) 2021-2023 relatedness data
    
    vcf_temp_df <- read.vcfR(
      here('raw_data/hirundo_rustica.recalc.batch2.mask.CO-2023.maf05.thin10k.vcf.gz'))
  
    

      
  
###############################################################################
##############                  6. Export data                   ##############
###############################################################################
      
  ### 6.1 Export data to an RData file     
    ## a) Save and export raw repro. success data 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      save(file = here('data/1_pheno_reprod_df.RData'), 
           list = c('pheno_reprod_df'))
      

    ## b) Save and export genomic meta data
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      save(file = here('data/1_genome_meta_df.RData'), 
           list = c('genome_meta_df'))
      
    ## c) Save and export raw relatedness data
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      save(file = here('data/2_relate_df.RData'), 
           list = c('relate_df'))
      
    ## d) Save and export raw genomics data
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      save(file = here('data/3_vcf_temp_df.RData'), 
           list = c('vcf_temp_df'))
      
   
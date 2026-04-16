################################################################################
#############                                                      #############
#############     Assessing the causal effect of plumage traits    #############
#############     on reproduction using Mendelian randomization    #############   
#############                                                      #############
#############                     1. Tidy data                     #############
#############                                                      #############
#############                   By: Zach Laubach                   #############
#############                 created: 9 Sept 2025                 #############
#############               last updated: 9 Sept 2025              #############
################################################################################


  ### PURPOSE: Tidy the 
  
  # Code Blocks
    # 1. Configure work space
    # 2. Import data 
    # 3. Tidy attribute data frame
    # 4. Export data
  


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
    
    # load lubridate package
      library ('lubridate')

    
    # load here package
      library ('here')
    
    
  ### 1.3 Get Version and Session Info
    R.Version()
    sessionInfo()
    
    # Developed in:   
    # R version 4.4.2 (2024-10-31)
    # Platform: x86_64-apple-darwin20
    # Running under: macOS Sequoia 15.1.
    
  
  ### 1.4 Set working directory 
    setwd(here())
  
  
  ### 1.5 Source functions
    ## a) Source scripts path
      source_path <- paste("~/WD/Git/source_code/")
    
    ## b) all_char_to_lower function
      source(file = paste0(source_path, "all_char_to_lower.R"))
    
    ## c) format_var_names function
      source(file = paste0(source_path, "format_var_names.R"))
      
      

###############################################################################
##############                  2. Import data                   ##############
###############################################################################    

  ### 2.1 Import sample data files
    ## a) Load phenotype and reproductive data 
      load(here('data/1_pheno_reprod_df.RData'))
      
    ## b) Load genome meta data 
      load(here('data/1_genome_meta_df.RData'))
      
    ## c) Load relatedness data
      load(here('data/1_relate_df.RData'))

###############################################################################
##############          3. Tidy attribute data frame             ##############
############################################################################### 
  
  ### 3.1 Basic tidying genome_meta_df
    ## a) Format data to all lower case
      genome_meta_df <- AllCharactersToLower(genome_meta_df)
      
    ## b) Format variable names
      genome_meta_df <- FormatVarNames(genome_meta_df)
      
    ## c) Filter adults from genome_meta_df
      genome_meta_df <- genome_meta_df  %>%
        filter(!is.na(type)) %>%
        filter(type != 'chick' & type != 'egg' & type != 'kid')
      
  ### 3.2 Basic tidying pheno_reprod_df
    ## a) Format data to all lower case
      pheno_reprod_df <- AllCharactersToLower(pheno_reprod_df)
      
    ## b) Format variable names
      pheno_reprod_df <- FormatVarNames(pheno_reprod_df)
      
    ## c) Filter adults from pheno_reprod_df
      pheno_reprod_df <- pheno_reprod_df  %>%
        filter(!is.na(year)) %>%
        filter(year >= 2021 & year <= 2023) 
      
    ## c) left join pheno_reprod_df to genome_meta_df
      genome_pheno_repro_df <- genome_meta_df  %>%
        left_join(pheno_reprod_df,
                  by = c('band' = 'band', 'year' = 'year'), 
                  copy = F) 
    
    ## d) filter samples with no reproductive data 
      genome_pheno_repro_df <- genome_pheno_repro_df %>%
        filter(!is.na(eggs.1) | !is.na(eggs.2) | !is.na(eggs.3) | 
                 !is.na(eggs.4))
      
    ## e) sample size
      n_tot = length(unique(genome_pheno_repro_df$band))
    
    ## f) repeats
      repeated_ids <- genome_pheno_repro_df %>%
        pull(band) %>%
        duplicated() %>%
        {genome_pheno_repro_df$band[.]}
      
      
    

###############################################################################
##############                  4. Export data                   ##############
###############################################################################
      
  ### 4.1 Export data to an RData file 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      
    ## a) Save and export data for CHR 2015 pre-manipulation >= 20 RSSI
      save(file = here('data/2_genome_pheno_repro_df.RData'), 
           list = c('genome_pheno_repro_df'))
      
      

    

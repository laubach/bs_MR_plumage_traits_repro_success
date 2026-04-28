################################################################################
#############                                                      #############
#############     Assessing the causal effect of plumage traits    #############
#############     on reproduction using Mendelian randomization    #############   
#############                                                      #############
#############            3. Descriptive stats summaries            #############
#############                                                      #############
#############                   By: Zach Laubach                   #############
#############                 created: 9 Sept 2025                 #############
#############             last updated: 27 April 2026              #############
################################################################################


  ### PURPOSE: Tidy the 
  
  # Code Blocks
    # 1. Configure work space
    # 2. Import data 
    # 3. Descriptive stats relatedness and snps
    # 4. Descriptive stats traits
    # 5. Descriptive stats reproduction
    # 6. Export data
  


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
    
    # # load lubridate package
    #   library ('lubridate')

    
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
    ## a) load snps, phenotype, and reproduction data
      load(here('data/2-4_snps_pheno_repro_df.RData'))
      
    ## b) Load relatedness data
      load(here('data/3-4_relate_df.RData'))

      
      
###############################################################################
##############     3. Descriptive stats relatedness and snps     ##############
############################################################################### 
  
  ### 3.1  
    ## a) 
      
      
      
###############################################################################
##############            4. Descriptive stats traits            ##############
############################################################################### 
      
      ### 4.1  
      ## a) 
      
      
      
###############################################################################
##############        5. Descriptive stats reproduction          ##############
############################################################################### 
      
      ### 5.1  
      ## a)  
      
      

###############################################################################
##############                  6. Export data                   ##############
###############################################################################
      
  ### 6.1 Export data to an RData file 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      
    ## a) Save and export relatedness data
      save(file = here('data/4_MR_data.RData'), 
           list = c('snps_pheno_repro_df', 'relate_mat'))
      
      

    

################################################################################
#############                                                      #############
#############     Assessing the causal effect of plumage traits    #############
#############     on reproduction using Mendelian randomization    #############   
#############                                                      #############
#############           2. Tidy genetic relatedness data           #############
#############                                                      #############
#############                   By: Zach Laubach                   #############
#############                 created: 9 Sept 2025                 #############
#############              last updated: 10 June 2026              #############
################################################################################


  ### PURPOSE: Tidy the genetic relatedness data 
  
  # Code Blocks
    # 1. Configure work space
    # 2. Import data 
    # 3. Tidy relatedness data frame 
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
    
    # # load lubridate package
    #   library ('lubridate')

    
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
      load(here('data/2_relate_df.RData'))
      
      

###############################################################################
##############          3. Tidy relatedness data frame           ##############
############################################################################### 
  
  ### 3.1  Basic tidying relate_d
    ## a) Format variable names
      relate_df <- FormatVarNames(relate_df)
      
    ## b) intersect relate_df with individuals from snps_pheno_repro_df
      relate_coef_df <- relate_df %>%
        filter (ind1 %in% snps_pheno_repro_df$sample.id & 
                ind2 %in% snps_pheno_repro_df$sample.id)
    
    ## f) tidy relate_coef_df; (pi_hat is the relatedness coef)  
      relate_coef_df <- relate_coef_df %>%
        rename('coef.relate' = 'pi.hat') %>%
        select(c(ind1, ind2, coef.relate))
      
    
    ## g) Check if there are 270 individuals in both columns
      length(unique(relate_coef_df$ind1))
      length(unique(relate_coef_df$ind2))
           
    ## h) Check if all individuals are in both ind1 and ind2 columns
      # to generate the matrix from the pairwise lists, each ID must occur at 
      # least once in each ind. ID column
        ind1.unique <- setdiff(relate_coef_df$ind1, relate_coef_df$ind2)
        ind2.unique <- setdiff(relate_coef_df$ind2, relate_coef_df$ind1)
    
      # identify which row has both unique ind ID
        row <- which(relate_coef_df$ind1 == ind1.unique & 
                       relate_coef_df$ind2 == ind2.unique)
          
      # Swap ind1 with ind2 at the shared row
        relate_coef_df[row, c('ind1', 'ind2')] <- 
          relate_coef_df[row, c('ind2', 'ind1')]
      
      
    # ## g) concatenate ind1 and ind2
    #   relate_coef_df$pair.id <- paste0(relate_coef_df$ind1,'_',relate_coef_df$ind2)
    # 
    # ## h) check for repeat pair.id
    #   repeated_pair_ids <- relate_coef_df %>%
    #     pull(pair.id) %>%
    #     duplicated() %>%
    #     {relate_coef_df$pair.id[.]}
      # 0 repeats
      
    # ## i) average coef.relate when estimated in multiple years and tidy
    #   relate_coef_df <- relate_coef_df %>% 
    #     group_by(pair.id) %>% 
    #     mutate(coef.relate = mean(coef.relate)) %>%
    #     distinct(pair.id, .keep_all=TRUE) %>%
    #     ungroup() %>%
    #     #summarise_at(vars(matches('coef.relate')), mean)
    #     select(-c(pair.id))
    
    # ## i) remove pair ID
    #   relate_coef_df <- relate_coef_df %>%
    #     select(-c(pair.id))
      
    ## i) convert to a simple data frame
      relate_coef_df <- as.data.frame(relate_coef_df)
      class(relate_coef_df)

    ## j) transform long (pairwise) df to wide using xtabs - makes an asymmetric matrix
      relate_mat_asym <- xtabs(coef.relate ~ ind1 + ind2, data = relate_coef_df)

    ## k) convert to a matrix instead of xtab object
      relate_mat_asym <- unclass(relate_mat_asym)
    
    ## l) add square matrix by its transpose to get a symmetrical matrix
      relate_mat <- relate_mat_asym + t(relate_mat_asym)
      
    ## m) # Add diagonal values (setting self-relatedness to 1.0)
      diag(relate_mat) <- 1.0

  
  ### 3.2 Subset relatedness matrices to match phenotype data frames
    ## a) Create a subset relatedness matrix matching snps_plum_repro_df
      relate_mat_plum_snps <- relate_mat[rownames(relate_mat) %in% 
                                            snps_plum_repro_df$sample.id, 
                                          colnames(relate_mat) %in% 
                                           snps_plum_repro_df$sample.id]
      
    ## b) Create a subset relatedness matrix matching snps_TS_repro_df
      relate_mat_TS_snps <- relate_mat[rownames(relate_mat) %in% 
                                            snps_TS_repro_df$sample.id, 
                                          colnames(relate_mat) %in% 
                                            snps_TS_repro_df$sample.id]   
      
    ## c) Create a subset relatedness matrix matching snps_plum_repro_f_df
      relate_mat_plum_repro_f <- relate_mat[rownames(relate_mat) %in%
                                              snps_plum_repro_f_df$sample.id, 
                                            colnames(relate_mat) %in%
                                              snps_plum_repro_f_df$sample.id]
      
    ## d) Create a subset relatedness matrix matching snps_plum_repro_m_df
      relate_mat_plum_repro_m <- relate_mat[rownames(relate_mat) %in%
                                            snps_plum_repro_m_df$sample.id, 
                                          colnames(relate_mat) %in%
                                            snps_plum_repro_m_df$sample.id]
      
    ## e) Create a subset relatedness matrix matching snps_TS_repro_f_df
      relate_mat_TS_repro_f <- relate_mat[rownames(relate_mat) %in%
                                             snps_TS_repro_f_df$sample.id, 
                                           colnames(relate_mat) %in%
                                             snps_TS_repro_f_df$sample.id]
      
    ## f) Create a subset relatedness matrix matching snps_TS_repro_m_df
      relate_mat_TS_repro_m <- relate_mat[rownames(relate_mat) %in%
                                          snps_TS_repro_m_df$sample.id, 
                                        colnames(relate_mat) %in%
                                          snps_TS_repro_m_df$sample.id]
      
      

  # ### 3.5 toy data example tidying relate_df
  #   ## convert relate_coef_df (pi_hat is the relatedness coef)
  #   
  #   pairwise_data <- data.frame(
  #     ID1 = c("A1", "C1", "A1", "B1", "B1", "D1"),
  #     ID2 = c("B1", "A1", "D1", "C1", "D1", "C1"),
  #     Relatedness = c(0.85, 0.45, 0.96, 0.85, 0.56, 0.45)
  #   )
  #   
  #   # pairwise_data <- data.frame(
  #   #   ID1 = c("A1", "A1", "A1", "B1", "B1", "C1"),
  #   #   ID2 = c("B1", "C1", "D1", "C1", "D1", "D1"),
  #   #   Relatedness = c(0.85, 0.45, 0.96, 0.85, 0.56, 0.45)
  #   # )
  #   
  #   #dist_matrix2 <- list2dist(pairwise_data)
  #   r_matrix_asym <- xtabs(Relatedness ~ ID1 + ID2, data = pairwise_data)
  #   
  #   # add square matrix by its transpose to get a symmetrical matrix
  #   r_matrix <- r_matrix_asym + t(r_matrix_asym)
  # 
  #   # # DOESN'T WORK
  #   # # convert relate_mat_asym to a square matrix using 'Matrix'
  #   # # NOTE forceSymmetric fills in the missing triangle to make it a complete square matrix
  #   # r_matrix <- forceSymmetric(as.matrix(r_matrix_asym), uplo = "U")
  #   
  #   #r_matrix <- as.matrix(r_matrix)
  #   r_matrix <- unclass(r_matrix)
  #   class(r_matrix)
  #   
  #   # Add diagonal values (setting self-relatedness to 1.0)
  #   diag(r_matrix) <- 1.0
  #   
    

###############################################################################
##############                  4. Export data                   ##############
###############################################################################
      
  ### 4.1 Export data to an RData file 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      
    ## a) Save and export relatedness data
      save(file = here('data/3-5_relate_df.RData'), 
           list = c('relate_mat', 'relate_coef_df', 
                    'relate_mat_plum_snps', 'relate_mat_TS_snps'
                    , 'relate_mat_plum_repro_f', 'relate_mat_plum_repro_m',
                    'relate_mat_TS_repro_f', 'relate_mat_TS_repro_m'
                    ))
      
      

    

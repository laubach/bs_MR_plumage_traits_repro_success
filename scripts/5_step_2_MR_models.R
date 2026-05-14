################################################################################
#############                                                      #############
#############     Assessing the causal effect of plumage traits    #############
#############     on reproduction using Mendelian randomization    #############   
#############                                                      #############
#############       5. Step 2: Mendelian randomization models      #############
#############                                                      #############
#############                   By: Zach Laubach                   #############
#############                 created: 9 Sept 2025                 #############
#############              last updated: 12 May 2026               #############
################################################################################


  ### PURPOSE: Tidy the 
  
  # Code Blocks
    # 1. Configure work space
    # 2. Import data 
    # 3. Model fit/diagnostic checks 
    # 4. Step 2: Female fit values on reproduction
    # 5. Step 2: Male fit values on reproduction
    # 7. Export data
  
# Update the import/export file numbers
# What should be used for weights (cohen's D or effect sizes from Drew's paper
# Include relatedness in both steps
# is fecundity the best/available measure of reproductive success
# shoud we run models stratified by sex?
# should we exclude or just run sensitivity analyses for experiment?
# Is the average brightness for all ventral measures appropriate?
# What is the best measure of TS (avg, longest, right?)
# Are color measures standardized across years?

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
    
      # load here package
        library ('here')
      
    ## b) Modelling Packages    
      # load glmmTMB package
        library ('glmmTMB')
      
      # load bbmle package
        library('bbmle') # for AICtab
      
      # load DHARMa package
        library('DHARMa')
      
      # load performance package
        library('performance')
      
    ## c) Graphing/visualization packages        
      # load ggplot2 package
        library('ggplot2')
      
    
      # load gridExtra package 
        library ('gridExtra')  
      
      # load rempsyc to make publication ready tables
      # see https://cran.r-project.org/web/packages/rempsyc/vignettes/table.html
        library('rempsyc')
      
      # pkgs <- c('flextable', 'broom', 'report', 'effectsize')
      # install_if_not_installed(pkgs)
      
      # load flextable package 
        library ('flextable')  
      
      # load broom package
        #library ('broom')
      
      # load broom.mixed package 
        library ('broom.mixed')
      
      # load report package 
        library ('report')
      
      # load effectsize package 
        #library ('effectsize')
      
    
    
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
    ## a) load snps, phenotype, reproduction data, and relatedness data
      load(here('data/5_pred_vals_repro_data.RData'))
      

      
###############################################################################
##############          3. Model fit/diagnostic checks           ##############
###############################################################################         


      
###############################################################################
##############   5. Step 2: Female fit values on reproduction    ##############
###############################################################################              
    
  # For Step 2: Female models, we will regress measures of reproductive
      # success on the predicted values of each plumage trait
      # (e.g., reproductive success ~ predicted trait values) 
      # to get the unconfounded effect of traits on reproductive success.
      

  ### 5.1 SNPs predicted breast brightness on reproduction - female models
    ## a) glmmTMB model for effect of predicted breast brightness on the 
        # number of nestlings fledged from nest attempts 1 and 2
      snps_pred_plum_repro_zinf_nb1_f_mod <-  glmmTMB(sum.fecundity ~ 
                            scale(snps.plum.pred.val) +
                          # proportional covariance matrix based on relatedness
                            propto(0 + sample.id | 1, relate_mat_plum_repro_f),
                            data = snps_plum_repro_f_df,
                            ziformula=~1,
                            family = nbinom1)
      
    ## b) model diagnostics
      performance::check_overdispersion(snps_pred_plum_repro_zinf_nb1_f_mod) 
      performance::check_zeroinflation(snps_pred_plum_repro_zinf_nb1_f_mod)
    
    ## c) model summary
      summary(snps_pred_plum_repro_zinf_nb1_f_mod)
      
    ## d) tidy effect estimates
      snps_pred_plum_repro_zinf_nb1_f_tidy <- 
        tidy(snps_pred_plum_repro_zinf_nb1_f_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = T) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value'))  
      
    ## e) label model
      snps_pred_plum_repro_zinf_nb1_f_tidy$model <-
        c('SNPs instrument - female')
      
    ## f) tidy snps_pred_plum_repro_zinf_nb1_f_tidy by making a new 
        # column '95% CI'
      snps_pred_plum_repro_zinf_nb1_f_tidy$'95% CI' = 
        paste(round(snps_pred_plum_repro_zinf_nb1_f_tidy$conf.low, 3),
              round(snps_pred_plum_repro_zinf_nb1_f_tidy$conf.high, 3), 
              sep=", ")
      
    ## g) tidy snps_pred_plum_repro_zinf_nb1_f_tidy by removing extra 
        # columns and rows
      snps_pred_plum_repro_zinf_nb1_f_tidy <- 
        snps_pred_plum_repro_zinf_nb1_f_tidy %>%
        rename(c(`Rate Ratio (exponentiated)` = `estimate`, 
                 `P-value` = `p.value`)) %>%
        select(c(`Rate Ratio (exponentiated)`, `95% CI`, `P-value`, `model`))
      
    ## h) move 'model' to first column
      snps_pred_plum_repro_zinf_nb1_f_tidy <- 
        snps_pred_plum_repro_zinf_nb1_f_tidy %>%
        relocate(model, .before = `Rate Ratio (exponentiated)`) 
    
      
  ### 5.2 SNPs predicted breast brightness on reproduction - male models
    ## a) glmmTMB model for effect of predicted breast brightness on the 
      # number of nestlings fledged from nest attempts 1 and 2
      snps_pred_plum_repro_zinf_nb1_m_mod <-  glmmTMB(ep.repro ~ 
                          scale(snps.plum.pred.val) +
                          # proportional covariance matrix based on relatedness
                          propto(0 + sample.id | 1, relate_mat_plum_repro_m),
                          data = snps_plum_repro_m_df,
                          ziformula=~1,
                          family = nbinom1)
      
    ## b) model diagnostics
      performance::check_overdispersion(snps_pred_plum_repro_zinf_nb1_m_mod) 
      performance::check_zeroinflation(snps_pred_plum_repro_zinf_nb1_m_mod)
      
    ## c) model summary
      summary(snps_pred_plum_repro_zinf_nb1_m_mod)
      
    ## d) tidy effect estimates
      snps_pred_plum_repro_zinf_nb1_m_tidy <- 
        tidy(snps_pred_plum_repro_zinf_nb1_m_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = T) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value'))  
      
    ## e) label model
      snps_pred_plum_repro_zinf_nb1_m_tidy$model <-
        c('SNPs instrument - male')
    
    ## f) tidy snps_pred_plum_repro_zinf_nb1_m_tidy by making a new 
        # column '95% CI'
      snps_pred_plum_repro_zinf_nb1_m_tidy$'95% CI' = 
        paste(round(snps_pred_plum_repro_zinf_nb1_m_tidy$conf.low, 3),
              round(snps_pred_plum_repro_zinf_nb1_m_tidy$conf.high, 3), 
              sep=", ")
      
    ## g) tidy snps_pred_plum_repro_zinf_nb1_m_tidy by removing extra 
        # columns and rows
      snps_pred_plum_repro_zinf_nb1_m_tidy <- 
        snps_pred_plum_repro_zinf_nb1_m_tidy %>%
        rename(c(`Rate Ratio (exponentiated)` = `estimate`, 
                 `P-value` = `p.value`)) %>%
        select(c(`Rate Ratio (exponentiated)`, `95% CI`, `P-value`, `model`))
      
    ## h) move 'model' to first column
      snps_pred_plum_repro_zinf_nb1_m_tidy <- 
        snps_pred_plum_repro_zinf_nb1_m_tidy %>%
        relocate(model, .before = `Rate Ratio (exponentiated)`) 
    
      
  ### 5.3 PGS predicted breast brightness on reproduction - female models
    ## a) glmmTMB model for effect of predicted breast brightness on the 
      # number of nestlings fledged from nest attempts 1 and 2
      pgs_pred_plum_repro_zinf_nb1_f_mod <-  glmmTMB(sum.fecundity ~ 
                            scale(pgs.plum.pred.val) +
                        # proportional covariance matrix based on relatedness
                            propto(0 + sample.id | 1, relate_mat_plum_repro_f),
                            data = snps_plum_repro_f_df,
                            ziformula=~1,
                              family = nbinom1)
      
    ## b) model diagnostics
      performance::check_overdispersion(pgs_pred_plum_repro_zinf_nb1_f_mod) 
      performance::check_zeroinflation(pgs_pred_plum_repro_zinf_nb1_f_mod)
      
    ## c) model summary
      summary(pgs_pred_plum_repro_zinf_nb1_f_mod)
      
    ## d) tidy effect estimates
      pgs_pred_plum_repro_zinf_nb1_f_tidy <- 
        tidy(pgs_pred_plum_repro_zinf_nb1_f_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = T) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value'))  
      
    ## e) label model
      pgs_pred_plum_repro_zinf_nb1_f_tidy$model <-
        c('PGS instrument - female')
      
    ## f) tidy pgs_pred_plum_repro_zinf_nb1_f_tidy by making a new 
      # column '95% CI'
      pgs_pred_plum_repro_zinf_nb1_f_tidy$'95% CI' = 
        paste(round(pgs_pred_plum_repro_zinf_nb1_f_tidy$conf.low, 3),
              round(pgs_pred_plum_repro_zinf_nb1_f_tidy$conf.high, 3), 
              sep=", ")
      
    ## g) tidy pgs_pred_plum_repro_zinf_nb1_f_tidy by removing extra 
      # columns and rows
      pgs_pred_plum_repro_zinf_nb1_f_tidy <- 
        pgs_pred_plum_repro_zinf_nb1_f_tidy %>%
        rename(c(`Rate Ratio (exponentiated)` = `estimate`, 
                 `P-value` = `p.value`)) %>%
        select(c(`Rate Ratio (exponentiated)`, `95% CI`, `P-value`, `model`))
      
    ## h) move 'model' to first column
      pgs_pred_plum_repro_zinf_nb1_f_tidy <- 
        pgs_pred_plum_repro_zinf_nb1_f_tidy %>%
        relocate(model, .before = `Rate Ratio (exponentiated)`) 
      
      
  ### 5.4 PGS predicted breast brightness on reproduction - male models
    ## a) glmmTMB model for effect of predicted breast brightness on the 
      # number of nestlings fledged from nest attempts 1 and 2
      pgs_pred_plum_repro_zinf_nb1_m_mod <-  glmmTMB(ep.repro ~ 
                            scale(pgs.plum.pred.val) +
                          # proportional covariance matrix based on relatedness
                            propto(0 + sample.id | 1, relate_mat_plum_repro_m),
                            data = snps_plum_repro_m_df,
                            ziformula=~1,
                              family = nbinom1)
      
    ## b) model diagnostics
      performance::check_overdispersion(pgs_pred_plum_repro_zinf_nb1_m_mod) 
      performance::check_zeroinflation(pgs_pred_plum_repro_zinf_nb1_m_mod)
      
    ## c) model summary
      summary(pgs_pred_plum_repro_zinf_nb1_m_mod)
      
    ## d) tidy effect estimates
      pgs_pred_plum_repro_zinf_nb1_m_tidy <- 
        tidy(pgs_pred_plum_repro_zinf_nb1_m_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = T) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value'))  
      
    ## e) label model
      pgs_pred_plum_repro_zinf_nb1_m_tidy$model <-
        c('PGS instrument - male')
      
    ## f) tidy pgs_pred_plum_repro_zinf_nb1_m_tidy by making a new 
      # column '95% CI'
      pgs_pred_plum_repro_zinf_nb1_m_tidy$'95% CI' = 
        paste(round(pgs_pred_plum_repro_zinf_nb1_m_tidy$conf.low, 3),
              round(pgs_pred_plum_repro_zinf_nb1_m_tidy$conf.high, 3), 
              sep=", ")
      
    ## g) tidy pgs_pred_plum_repro_zinf_nb1_m_tidy by removing extra 
      # columns and rows
      pgs_pred_plum_repro_zinf_nb1_m_tidy <- 
        pgs_pred_plum_repro_zinf_nb1_m_tidy %>%
        rename(c(`Rate Ratio (exponentiated)` = `estimate`, 
                 `P-value` = `p.value`)) %>%
        select(c(`Rate Ratio (exponentiated)`, `95% CI`, `P-value`, `model`))
      
    ## h) move 'model' to first column
      pgs_pred_plum_repro_zinf_nb1_m_tidy <- 
        pgs_pred_plum_repro_zinf_nb1_m_tidy %>%
        relocate(model, .before = `Rate Ratio (exponentiated)`)    
      

  ### 5.5 Make formatted tables of the associations between the instrument
        # breast brightness
    # a) Combine Step 2: Female regression estimates into a tidy table
      trait_tot_repro_smmry <- 
        bind_rows(snps_pred_plum_repro_zinf_nb1_f_tidy,
                  snps_pred_plum_repro_zinf_nb1_m_tidy,
                  pgs_pred_plum_repro_zinf_nb1_f_tidy,
                  pgs_pred_plum_repro_zinf_nb1_m_tidy)
      
    ## b) format table creating flextable object 
      trait_tot_repro_smmry_frmt <- nice_table(trait_tot_repro_smmry, 
                                            separate.header = T)
      
    ## c) preview the table as a word doc    
      #print(trait_tot_repro_smmry_frmt, preview = 'docx')
      
    ## d) save trait_tot_repro_smmry_frmt as a word doc
      flextable::save_as_docx(trait_tot_repro_smmry_frmt, 
                      path = here('output/trait_tot_repro_smmry_frmt.docx')) 
      
      
      

###############################################################################
##############    6. Results visualisation    ##############
###############################################################################         

  ### 6.1 graphs
    

###############################################################################
##############                  4. Export data                   ##############
###############################################################################
      
  ### 4.1 Export data to an RData file 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      
    ## a) Save and export relatedness data
      save(file = here('data/3_relatedness_data.RData'), 
           list = c('relate_mat', 'relate_coef_df'))
      
      

    

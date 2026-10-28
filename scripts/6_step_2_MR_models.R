################################################################################
#############                                                      #############
#############     Assessing the causal effect of plumage traits    #############
#############     on reproduction using Mendelian randomization    #############   
#############                                                      #############
#############       5. Step 2: Mendelian randomization models      #############
#############                                                      #############
#############                   By: Zach Laubach                   #############
#############                 created: 9 Sept 2025                 #############
#############              last updated: 10 June 2026              #############
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
    
      # load patchwork packages
        library('patchwork')  
    
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
    ## a) load snps, phenotype, reproduction data, and relatedness data
      load(here('data/5_pred_vals_repro_data.RData'))
      
    ## b) load natural plumage reproduction models for graphing
      load(here('data/5_nat_plum_repro_mods.RData'))
      

      
###############################################################################
##############          3. Model fit/diagnostic checks           ##############
###############################################################################         


      
###############################################################################
##############   3. Step 2: Female fit values on reproduction    ##############
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
      
      
      
      ###############################################################################
      ##############   4. Step 2: Mmale fit values on reproduction    ##############
      ###############################################################################      
     
      
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
      

      
      
      
      
      
      
      ###############################################################################
      ##############   5. Step 2: Tidy    ##############
      ###############################################################################        
      
      
  ### 5.5 Make formatted tables of the associations between the instrument
        # breast brightness
    ## a) Combine Step 2: Female regression estimates into a tidy table
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

  ### 6.1 Format models summaries for graphing
    # a) Tidy female_plum_repro_sum for graphing
      # rename model 
      female_plum_repro_sum[3, 1] <- 'natural color'
      # filter to retain on specific mode
      female_plum_repro_sum <- female_plum_repro_sum %>%
        filter(model == 'natural color')
      
    # b) Tidy male_plum_repro_sum for graphing
      # rename model 
      male_plum_repro_sum[2, 1] <- 'natural color'
      # filter to retain on specific mode
      male_plum_repro_sum <- male_plum_repro_sum %>%
        filter(model == 'natural color')
      
    ## c) Combine models tidy table
      trait_tot_repro_graphing_tbl <- 
        bind_rows(female_plum_repro_sum,
                  snps_pred_plum_repro_zinf_nb1_f_tidy,
                  pgs_pred_plum_repro_zinf_nb1_f_tidy,
                  male_plum_repro_sum,
                  snps_pred_plum_repro_zinf_nb1_m_tidy,
                  pgs_pred_plum_repro_zinf_nb1_m_tidy)
      
    ## d) add a variable 'sex' for subsestting graphs
      trait_tot_repro_graphing_tbl$sex <- c('female', 'female', 'female', 
                                            'male', 'male', 'male')
      
    ## e) separate 95% CI into lower and upper
      trait_tot_repro_graphing_tbl <-  trait_tot_repro_graphing_tbl %>%
      separate_wider_delim(
        cols = `95% CI`, 
        delim = ', ', 
        names = c('95% CI Lower', '95% CI Upper')) %>%
        mutate(across(c(`95% CI Lower`, `95% CI Upper`), as.numeric))
      
    # b) Tidy trait_tot_repro_graphing_tbl for graphing
      # rename model 
      trait_tot_repro_graphing_tbl[2, 1] <- 'SNPs instrument'
      trait_tot_repro_graphing_tbl[3, 1] <- 'PGS instrument'
      trait_tot_repro_graphing_tbl[5, 1] <- 'SNPs instrument'
      trait_tot_repro_graphing_tbl[6, 1] <- 'PGS instrument'
      
    # c) set the order of the models
      trait_tot_repro_graphing_tbl$model <- 
        factor(trait_tot_repro_graphing_tbl$model, 
               levels = c('natural color', 'SNPs instrument', 
                          'PGS instrument'))

      
      levels(trait_tot_repro_graphing_tbl$model)
      
      
  ### 6.2 Graph model summaries     
    ## a) Graph results of female natural vs genetic instrument models of 
      # breast brightness and total fecundity
      f_natural_vs_genetic_inst_plot <- 
        ggplot(data = (trait_tot_repro_graphing_tbl %>%
                         filter(`sex` == 'female')), 
               aes(x = `model`, 
                   y = `Rate Ratio (exponentiated)`, 
                   color = `model`)) +
        ylim(0.25, 2) +
        #facet_wrap(~ sex) +
        
        geom_hline(yintercept = 1, color = 'red',
                   linetype = 2, linewidth = 2) + # line at null behind coefs
        geom_point(size = 6) +
        geom_errorbar(aes(ymin=(`95% CI Lower`), 
                          ymax=(`95% CI Upper`)), linewidth = 2,  linetype = 1,
                      width = 0.25, position = position_dodge(0.5)) +
        scale_color_manual(values=c('slateblue4', 'slateblue3', 
                                    'slateblue1'
        )) +
        #coord_flip() + # flip x and y axes
        labs(title = 'Natural plumage vs genetic instrument') +
        theme(plot.title = element_text(hjust = 0.5, size = 20)) + # center title
        theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size title and axes labels
        theme(text = element_text(size=22, face = 'bold')) +
        theme(legend.position = 'none') +
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = 'white')) +
        # add major axes
        theme(axis.line = element_line(colour = 'black',
                                       linewidth = 0.5, linetype = 'solid')) +
        # change axes font style, color, size, angle, margin, and legend
        theme(axis.text.x = element_text(face='bold', color='black', 
                                         size=18, angle=0,
                                         margin = margin(t = 10, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face='bold', color='black', 
                                         size=22, angle=0, 
                                         margin = margin(t = 0, r = 10, 
                                                         b = 0, l = 10)),
              legend.title=element_blank(),
              legend.text=element_text(size=14),
              legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(italic('Female models'))) +
        ylab(expression
             (bold('Reproduction Rate Ratio (95% CI)')))  
      
    ## b) view plot  
      print(f_natural_vs_genetic_inst_plot)
      
    ## c) Graph results of male natural vs genetic instrument models of 
      # breast brightness and total fecundity
      m_natural_vs_genetic_inst_plot <- 
        ggplot(data = (trait_tot_repro_graphing_tbl %>%
                         filter(`sex` == 'male')), 
               aes(x = `model`, 
                   y = `Rate Ratio (exponentiated)`, 
                   color = `model`)) +
        ylim(0.25, 2) +
        #facet_wrap(~ sex) +
        
        geom_hline(yintercept = 1, color = 'red',
                   linetype = 2, linewidth = 2) + # line at null behind coefs
        geom_point(size = 6) +
        geom_errorbar(aes(ymin=(`95% CI Lower`), 
                          ymax=(`95% CI Upper`)), linewidth = 2,  linetype = 1,
                      width = 0.25, position = position_dodge(0.5)) +
        scale_color_manual(values=c('springgreen4', 'springgreen3', 
                                    'springgreen1'
        )) +
        #coord_flip() + # flip x and y axes
        labs(title = 'Natural plumage vs genetic instrument') +
        theme(plot.title = element_text(hjust = 0.5, size = 20)) + # center title
        theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size title and axes labels
        theme(text = element_text(size=22, face = 'bold')) +
        theme(legend.position = 'none') +
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = 'white')) +
        # add major axes
        theme(axis.line = element_line(colour = 'black',
                                       linewidth = 0.5, linetype = 'solid')) +
        # change axes font style, color, size, angle, margin, and legend
        theme(axis.text.x = element_text(face='bold', color='black', 
                                         size=18, angle=0,
                                         margin = margin(t = 10, r = 0, 
                                                         b = 10, l = 0)),
              # axis.text.y = element_text(face='bold', color='black', 
              #                            size=22, angle=0, 
              #                            margin = margin(t = 0, r = 0, 
              #                                            b = 0, l = 10)),
              axis.title.y = element_blank(), # Removes the y-axis title (label)
              axis.text.y  = element_blank(), # Removes the tick labels
              legend.title=element_blank(),
              legend.text=element_text(size=14),
              legend.position = 'none', #c(0.91, 0.94),
              legend.key = element_blank()) +
        xlab(expression(italic('Male models'))) +
        ylab(expression
             (bold('Reproduction Rate Ratio (95% CI)'))) 
      
    ## d) view plot  
      print(m_natural_vs_genetic_inst_plot)
      
    ## e) make panel plot using patchwork
      panel_plot <- (f_natural_vs_genetic_inst_plot | 
                       m_natural_vs_genetic_inst_plot)
      
    ## f) view panel plot
      panel_plot
      
    ## g) save panel plot
      # use ggsave to save the plot
      ggsave('panel_plot.pdf', 
             plot = panel_plot, 
             device = NULL,
             path = paste0(here(),'/output'), 
             scale = 1, width = 16,
             height = 7,
             units = c('in'), dpi = 300, limitsize = TRUE)
      
      
     

      

    

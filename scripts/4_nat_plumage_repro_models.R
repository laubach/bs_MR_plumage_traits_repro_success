################################################################################
#############                                                      #############
#############     Assessing the causal effect of plumage traits    #############
#############     on reproduction using Mendelian randomization    #############   
#############                                                      #############
#############              4. Natural plumage models               #############
#############                                                      #############
#############                   By: Zach Laubach                   #############
#############                 created: 9 Sept 2025                 #############
#############               last updated: 20 Oct 2026              #############
################################################################################


  ### PURPOSE: Tidy the 
  
  # Code Blocks
    # 1. Configure work space
    # 2. Import data 
    # 3. Model natural plumage traits and reproduction
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
      
      # load MuMIn
        library('MuMIn')
      
    ## c) Graphing/visualization packages        
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
    ## a) load snps, phenotype, and reproduction data
      load(here('data/2-4_snps_pheno_repro_df.RData'))
      
    ## b) Load relatedness data
      load(here('data/3-5_relate_df.RData'))


      
###############################################################################
############## 3. Model natural plumage traits and reproduction ##############
###############################################################################      
    
  ### 3.1 Model fitting to identify the appropriate distribution for count
      # data
    ## a) Fit glmmTMB model (tot.repro ~ breast.avg.bright) 
      # with Poisson distribution
      snps_pred_plum_repro_poiss_mod <-  glmmTMB(tot.repro ~ 
                               breast.avg.bright
                               + sex  +
                         # proportional covariance matrix based on relatedness
                               propto(0 + sample.id | 1, relate_mat_plum_snps),
                                                 data = snps_plum_repro_df,
                                                 family = poisson)
      
      
    ## b) udpate the Poisson model to include zero-inflation
      snps_pred_plum_repro_zinf_poiss_mod <- 
        update(snps_pred_plum_repro_poiss_mod, ziformula=~1)
      
    ## c) udpate the Poisson model to fit a negative binomial 2 distribution
      # where variance = µ(1 + µ/k)
      snps_pred_plum_repro_nb2_mod <- 
        update(snps_pred_plum_repro_poiss_mod, family = nbinom2)
      
    ## d) udpate the negative binomial 2 model to include zero-inflation
      snps_pred_plum_repro_zinf_nb2_mod <- 
        update(snps_pred_plum_repro_nb2_mod, ziformula=~1)
      
    ## e) udpate the Poisson model to fit a negative binomial 1 distribution
      # where variance = ϕµ
      snps_pred_plum_repro_nb1_mod <- 
        update(snps_pred_plum_repro_poiss_mod, family = nbinom1)
      
    ## f) udpate the negative binomial 1 model to include zero-inflation
      snps_pred_plum_repro_zinf_nb1_mod <- 
        update(snps_pred_plum_repro_nb1_mod,  ziformula=~1)
      
    ## g) assess delta AIC using bbmle to determine best fitting model   
      AICtab(snps_pred_plum_repro_poiss_mod, 
             snps_pred_plum_repro_zinf_poiss_mod,
             snps_pred_plum_repro_nb2_mod, 
             snps_pred_plum_repro_zinf_nb2_mod, 
             snps_pred_plum_repro_nb1_mod, 
             snps_pred_plum_repro_zinf_nb1_mod)
      
      
  ### 3.2 check diagnostics of best fitting model using 'DHARMa' and 
      # 'performance' packages
    ## a)  simulate residuals
      resid_best_fit_mod <- 
        simulateResiduals(fittedModel = snps_pred_plum_repro_zinf_nb1_mod, 
                          n = 250)
      
    ## b) plot residuals
      plot(resid_best_fit_mod)
      
    ## c) assess over dispersion
      testDispersion(resid_best_fit_mod)
      performance::check_overdispersion(snps_pred_plum_repro_zinf_nb1_mod) 
      
    ## d) assess zero inflation
      testZeroInflation(resid_best_fit_mod)
      performance::check_zeroinflation(snps_pred_plum_repro_zinf_nb1_mod)
      
      
  ### 3.3 Models of natural breast brightness on reproduction - female models
    ## a) glmmTMB model for effect of natural breast brightness on the 
      # number of nestlings fledged from nest attempts 1 and 2
      breast_repro_zinf_nb1_f_mod <- glmmTMB(tot.repro ~
                      scale(breast.avg.bright) +
                      # proportional covariance matrix based on relatedness
                          propto(0 + sample.id | 1, relate_mat_plum_repro_f),
                          data = snps_plum_repro_f_df,
                          ziformula=~1,
                              family = nbinom1)
      
    ## b) model diagnostics
      performance::check_overdispersion(breast_repro_zinf_nb1_f_mod) 
      performance::check_zeroinflation(breast_repro_zinf_nb1_f_mod)
      
    ## c) model summary
      summary(breast_repro_zinf_nb1_f_mod)
      
    ## d) tidy effect estimates
      breast_repro_zinf_nb1_f_tidy <- 
        tidy(breast_repro_zinf_nb1_f_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = T) %>%
        filter(!(term == '(Intercept)'| term == 'sexm')) %>%
        select(c('term', 'estimate', #'std.error', 
                 'conf.low', 'conf.high', 
                 #'statistic', 
                 'p.value'))  
      
    ## e) label model
      breast_repro_zinf_nb1_f_tidy$model <-
        c('total reproduction ~ breast color')
      
    ## f) tidy breast_repro_zinf_nb1_f_tidy by making a new column '95% CI'
      breast_repro_zinf_nb1_f_tidy$'95% CI' = 
        paste(round(breast_repro_zinf_nb1_f_tidy$conf.low, 3),
              round(breast_repro_zinf_nb1_f_tidy$conf.high, 3), 
              sep=", ")
      
    ## g) tidy breast_repro_zinf_nb1_f_tidy by removing extra columns and rows
      breast_repro_zinf_nb1_f_tidy <- breast_repro_zinf_nb1_f_tidy %>%
        rename(c(`Rate Ratio (exponentiated)` = `estimate`, 
                 `P-value` = `p.value`)) %>%
        select(c(`Rate Ratio (exponentiated)`, `95% CI`, `P-value`, `model`))
      
    ## h) move 'model' to first column
      breast_repro_zinf_nb1_f_tidy <- breast_repro_zinf_nb1_f_tidy %>%
        relocate(model, .before = `Rate Ratio (exponentiated)`) 
      
    ## i) glmmTMB model for effect of natural breast brightness on the 
      # number of EP nestlings fledged from nest attempts 1 and 2
      breast_ep_repro_zinf_nb1_f_mod <-  glmmTMB(ep.repro ~
                  scale(breast.avg.bright) +
                  # proportional covariance matrix based on relatedness
                  propto(0 + sample.id | 1, relate_mat_plum_repro_f),
                  data = snps_plum_repro_f_df,
                  ziformula=~1,
                    family = nbinom1)
      
    ## j) model diagnostics
      performance::check_overdispersion(breast_ep_repro_zinf_nb1_f_mod) 
      performance::check_zeroinflation(breast_ep_repro_zinf_nb1_f_mod)
      
    ## k) model summary
      summary(breast_ep_repro_zinf_nb1_f_mod)
      
    ## l) tidy effect estimates
      breast_ep_repro_zinf_nb1_f_tidy <- 
        tidy(breast_ep_repro_zinf_nb1_f_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = T) %>%
        filter(!(term == '(Intercept)'| term == 'sexm')) %>%
        select(c('term', 'estimate', #'std.error', 
                 'conf.low', 'conf.high', 
                 #'statistic', 
                 'p.value'))  
      
    ## m) label model
      breast_ep_repro_zinf_nb1_f_tidy$model <-
        c('EP ~ breast color')
      
    ## n) tidy breast_ep_repro_zinf_nb1_f_tidy by making a new column '95% CI'
      breast_ep_repro_zinf_nb1_f_tidy$'95% CI' = 
        paste(round(breast_ep_repro_zinf_nb1_f_tidy$conf.low, 3),
              round(breast_ep_repro_zinf_nb1_f_tidy$conf.high, 3), 
              sep=", ")
      
    ## o) tidy breast_ep_repro_zinf_nb1_f_tidy by removing extra columns and rows
      breast_ep_repro_zinf_nb1_f_tidy <- breast_ep_repro_zinf_nb1_f_tidy %>%
        rename(c(`Rate Ratio (exponentiated)` = `estimate`, 
                 `P-value` = `p.value`)) %>%
        select(c(`Rate Ratio (exponentiated)`, `95% CI`, `P-value`, `model`))
      
    ## p) move 'model' to first column
      breast_ep_repro_zinf_nb1_f_tidy <- breast_ep_repro_zinf_nb1_f_tidy %>%
        relocate(model, .before = `Rate Ratio (exponentiated)`) 
      
    ## q) glmmTMB model for effect of natural breast brightness on the 
      # number of eggs (fecundity) from nest attempts 1 and 2
      breast_fecund_zinf_nb1_f_mod <- glmmTMB(sum.fecundity ~
                  scale(breast.avg.bright) +
                  # proportional covariance matrix based on relatedness
                  propto(0 + sample.id | 1, relate_mat_plum_repro_f),
                  data = snps_plum_repro_f_df,
                  ziformula=~1,
                    family = nbinom1)
      
    ## r) model diagnostics
      performance::check_overdispersion(breast_fecund_zinf_nb1_f_mod) 
      performance::check_zeroinflation(breast_fecund_zinf_nb1_f_mod)
      
    ## s) model summary
      summary(breast_fecund_zinf_nb1_f_mod)
      
    ## t) tidy effect estimates
      breast_fecund_zinf_nb1_f_tidy <- 
        tidy(breast_fecund_zinf_nb1_f_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = T) %>%
        filter(!(term == '(Intercept)'| term == 'sexm')) %>%
        select(c('term', 'estimate', #'std.error', 
                 'conf.low', 'conf.high', 
                 #'statistic', 
                 'p.value'))  
      
    ## u label model
      breast_fecund_zinf_nb1_f_tidy$model <-
        c('fecundity ~ breast color')
      
    ## v) tidy breast_fecund_zinf_nb1_f_tidy by making a new column '95% CI'
      breast_fecund_zinf_nb1_f_tidy$'95% CI' = 
        paste(round(breast_fecund_zinf_nb1_f_tidy$conf.low, 3),
              round(breast_fecund_zinf_nb1_f_tidy$conf.high, 3), 
              sep=", ")
      
    ## w) tidy breast_fecund_zinf_nb1_f_tidy by removing extra columns and rows
      breast_fecund_zinf_nb1_f_tidy <- breast_fecund_zinf_nb1_f_tidy %>%
        rename(c(`Rate Ratio (exponentiated)` = `estimate`, 
                 `P-value` = `p.value`)) %>%
        select(c(`Rate Ratio (exponentiated)`, `95% CI`, `P-value`, `model`))
      
    ## x) move 'model' to first column
      breast_fecund_zinf_nb1_f_tidy <- breast_fecund_zinf_nb1_f_tidy %>%
        relocate(model, .before = `Rate Ratio (exponentiated)`) 
      
   
  ### 3.4 Models of natural tail streamer on reproduction - female models
    ## a) glmmTMB model for effect of natural tail streamer on the 
      # number of nestlings fledged from nest attempts 1 and 2
      TS_repro_zinf_nb1_f_mod <- glmmTMB(tot.repro ~
                            scale(rt.intact.avg.ts) +   
                        # proportional covariance matrix based on relatedness
                            propto(0 + sample.id | 1, relate_mat_plum_repro_f),
                            data = snps_plum_repro_f_df,
                            ziformula=~1,
                              family = nbinom1)
      
    ## b) model diagnostics
      performance::check_overdispersion(TS_repro_zinf_nb1_f_mod) 
      performance::check_zeroinflation(TS_repro_zinf_nb1_f_mod)
      
    ## c) model summary
      summary(TS_repro_zinf_nb1_f_mod)
      
    ## d) tidy effect estimates
      TS_repro_zinf_nb1_f_tidy <- 
        tidy(TS_repro_zinf_nb1_f_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = T) %>%
        filter(!(term == '(Intercept)'| term == 'sexm')) %>%
        select(c('term', 'estimate', #'std.error', 
                 'conf.low', 'conf.high', 
                 #'statistic', 
                 'p.value'))  
      
    ## e) label model
      TS_repro_zinf_nb1_f_tidy$model <-
        c('total reproduction ~ TS')
      
    ## f) tidy TS_repro_zinf_nb1_f_tidy by making a new column '95% CI'
      TS_repro_zinf_nb1_f_tidy$'95% CI' = 
        paste(round(TS_repro_zinf_nb1_f_tidy$conf.low, 3),
              round(TS_repro_zinf_nb1_f_tidy$conf.high, 3), 
              sep=", ")
      
    ## g) tidy TS_repro_zinf_nb1_f_tidy by removing extra columns and rows
      TS_repro_zinf_nb1_f_tidy <- TS_repro_zinf_nb1_f_tidy %>%
        rename(c(`Rate Ratio (exponentiated)` = `estimate`, 
                 `P-value` = `p.value`)) %>%
        select(c(`Rate Ratio (exponentiated)`, `95% CI`, `P-value`, `model`))
      
    ## h) move 'model' to first column
      TS_repro_zinf_nb1_f_tidy <- TS_repro_zinf_nb1_f_tidy %>%
        relocate(model, .before = `Rate Ratio (exponentiated)`) 
      
    ## i) glmmTMB model for effect of natural tail streamer on the 
      # number of EP nestlings fledged from nest attempts 1 and 2
      TS_ep_repro_zinf_nb1_f_mod <-  glmmTMB(ep.repro ~
                          scale(rt.intact.avg.ts) +
                       # proportional covariance matrix based on relatedness
                          propto(0 + sample.id | 1, relate_mat_plum_repro_f),
                          data = snps_plum_repro_f_df,
                          ziformula=~1,
                             family = nbinom1)
      
    ## j) model diagnostics
      performance::check_overdispersion(TS_ep_repro_zinf_nb1_f_mod) 
      performance::check_zeroinflation(TS_ep_repro_zinf_nb1_f_mod)
      
    ## k) model summary
      summary(TS_ep_repro_zinf_nb1_f_mod)
      
    ## l) tidy effect estimates
      TS_ep_repro_zinf_nb1_f_tidy <- 
        tidy(TS_ep_repro_zinf_nb1_f_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = T) %>%
        filter(!(term == '(Intercept)'| term == 'sexm')) %>%
        select(c('term', 'estimate', #'std.error', 
                 'conf.low', 'conf.high', 
                 #'statistic', 
                 'p.value'))  
      
    ## m) label model
      TS_ep_repro_zinf_nb1_f_tidy$model <-
        c('EP ~ TS')
      
    ## n) tidy TS_ep_repro_zinf_nb1_f_tidy by making a new column '95% CI'
      TS_ep_repro_zinf_nb1_f_tidy$'95% CI' = 
        paste(round(TS_ep_repro_zinf_nb1_f_tidy$conf.low, 3),
              round(TS_ep_repro_zinf_nb1_f_tidy$conf.high, 3), 
              sep=", ")
      
    ## o) tidy TS_ep_repro_zinf_nb1_f_tidy by removing extra columns and rows
      TS_ep_repro_zinf_nb1_f_tidy <- TS_ep_repro_zinf_nb1_f_tidy %>%
        rename(c(`Rate Ratio (exponentiated)` = `estimate`, 
                 `P-value` = `p.value`)) %>%
        select(c(`Rate Ratio (exponentiated)`, `95% CI`, `P-value`, `model`))
      
    ## p) move 'model' to first column
      TS_ep_repro_zinf_nb1_f_tidy <- TS_ep_repro_zinf_nb1_f_tidy %>%
        relocate(model, .before = `Rate Ratio (exponentiated)`) 
      
    ## q) glmmTMB model for effect of natural tail streamer on the 
      # number of eggs (fecundity) from nest attempts 1 and 2
      TS_fecund_zinf_nb1_f_mod <- glmmTMB(sum.fecundity ~
                            scale(rt.intact.avg.ts) +
                        # proportional covariance matrix based on relatedness
                            propto(0 + sample.id | 1, relate_mat_plum_repro_f),
                            data = snps_plum_repro_f_df,
                            ziformula=~1,
                              family = nbinom1)
      
    ## r) model diagnostics
      performance::check_overdispersion(TS_fecund_zinf_nb1_f_mod) 
      performance::check_zeroinflation(TS_fecund_zinf_nb1_f_mod)
      
    ## s) model summary
      summary(TS_fecund_zinf_nb1_f_mod)
      
    ## t) tidy effect estimates
      TS_fecund_zinf_nb1_f_tidy <- 
        tidy(TS_fecund_zinf_nb1_f_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = T) %>%
        filter(!(term == '(Intercept)'| term == 'sexm')) %>%
        select(c('term', 'estimate', #'std.error', 
                 'conf.low', 'conf.high', 
                 #'statistic', 
                 'p.value'))  
      
    ## u) label model
      TS_fecund_zinf_nb1_f_tidy$model <-
        c('fecundity ~ TS')
      
    ## v) tidy TS_fecund_zinf_nb1_f_tidy by making a new column '95% CI'
      TS_fecund_zinf_nb1_f_tidy$'95% CI' = 
        paste(round(TS_fecund_zinf_nb1_f_tidy$conf.low, 3),
              round(TS_fecund_zinf_nb1_f_tidy$conf.high, 3), 
              sep=", ")
      
    ## w) tidy TS_fecund_zinf_nb1_f_tidy by removing extra columns and rows
      TS_fecund_zinf_nb1_f_tidy <- TS_fecund_zinf_nb1_f_tidy %>%
        rename(c(`Rate Ratio (exponentiated)` = `estimate`, 
                 `P-value` = `p.value`)) %>%
        select(c(`Rate Ratio (exponentiated)`, `95% CI`, `P-value`, `model`))
      
    ## x) move 'model' to first column
      TS_fecund_zinf_nb1_f_tidy <- TS_fecund_zinf_nb1_f_tidy %>%
        relocate(model, .before = `Rate Ratio (exponentiated)`)  
      

  ### 3.5 Make formatted tables of the associations between natural
      # breast brightness and tail streamer with reproduction - females
    ## a) row bind breast_repro_zinf_nb1_f_tidy, 
      # breast_ep_repro_zinf_nb1_f_tidy, TS_repro_zinf_nb1_f_tidy, 
      # breast_fecund_zinf_nb1_f_tidy, TS_repro_zinf_nb1_f_tidy
      # TS_ep_repro_zinf_nb1_f_tidy, and TS_fecund_zinf_nb1_f_tidy
      female_plum_repro_sum <- bind_rows(breast_repro_zinf_nb1_f_tidy,
                                       breast_ep_repro_zinf_nb1_f_tidy,
                                       breast_fecund_zinf_nb1_f_tidy,
                                       TS_repro_zinf_nb1_f_tidy, 
                                       TS_ep_repro_zinf_nb1_f_tidy,
                                       TS_fecund_zinf_nb1_f_tidy) 
      
    ## b) format table creating flextable object 
      female_plum_repro_sum_frmt <- nice_table(female_plum_repro_sum, 
                                             separate.header = T)
      
    ## c) preview the table as a word doc    
      #print(female_plum_repro_sum_frmt, preview = 'docx')
      
    ## d) save female_plum_repro_sum_frmt as a word doc
      flextable::save_as_docx(female_plum_repro_sum_frmt, 
                        path = here('output/female_plum_repro_sum_frmt.docx'))  
      
      
      
  ### 3.6 Models of natural breast brightness on reproduction - male models
    ## a) glmmTMB model for effect of natural breast brightness on the 
      # number of nestlings fledged from nest attempts 1 and 2
      breast_repro_zinf_nb1_m_mod <- glmmTMB(tot.repro ~
                            scale(breast.avg.bright) +
                         # proportional covariance matrix based on relatedness
                            propto(0 + sample.id | 1, relate_mat_plum_repro_m),
                            data = snps_plum_repro_m_df,
                            ziformula=~1,
                              family = nbinom1)
      
    ## b) model diagnostics
      performance::check_overdispersion(breast_repro_zinf_nb1_m_mod) 
      performance::check_zeroinflation(breast_repro_zinf_nb1_m_mod)
      
    ## c) model summary
      summary(breast_repro_zinf_nb1_m_mod)
      
    ## d) tidy effect estimates
      breast_repro_zinf_nb1_m_tidy <- 
        tidy(breast_repro_zinf_nb1_m_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = T) %>%
        filter(!(term == '(Intercept)'| term == 'sexm')) %>%
        select(c('term', 'estimate', #'std.error', 
                 'conf.low', 'conf.high', 
                 #'statistic', 
                 'p.value'))  
      
    ## e) label model
      breast_repro_zinf_nb1_m_tidy$model <-
        c('total reproduction ~ breast color')
      
    ## f) tidy breast_repro_zinf_nb1_m_tidy by making a new column '95% CI'
      breast_repro_zinf_nb1_m_tidy$'95% CI' = 
        paste(round(breast_repro_zinf_nb1_m_tidy$conf.low, 3),
              round(breast_repro_zinf_nb1_m_tidy$conf.high, 3), 
              sep=", ")
      
    ## g) tidy breast_repro_zinf_nb1_m_tidy by removing extra columns and rows
      breast_repro_zinf_nb1_m_tidy <- breast_repro_zinf_nb1_m_tidy %>%
        rename(c(`Rate Ratio (exponentiated)` = `estimate`, 
                 `P-value` = `p.value`)) %>%
        select(c(`Rate Ratio (exponentiated)`, `95% CI`, `P-value`, `model`))
      
    ## h) move 'model' to first column
      breast_repro_zinf_nb1_m_tidy <- breast_repro_zinf_nb1_m_tidy %>%
        relocate(model, .before = `Rate Ratio (exponentiated)`) 
      
    ## i) glmmTMB model for effect of natural breast brightness on the 
      # number of EP nestlings fledged from nest attempts 1 and 2
      breast_ep_repro_zinf_nb1_m_mod <-  glmmTMB(ep.repro ~
                            scale(breast.avg.bright) +
                        # proportional covariance matrix based on relatedness
                            propto(0 + sample.id | 1, relate_mat_plum_repro_m),
                            data = snps_plum_repro_m_df,
                            ziformula=~1,
                              family = nbinom1)
      
    ## j) model diagnostics
      performance::check_overdispersion(breast_ep_repro_zinf_nb1_m_mod) 
      performance::check_zeroinflation(breast_ep_repro_zinf_nb1_m_mod)
      
    ## k) model summary
      summary(breast_ep_repro_zinf_nb1_m_mod)
      
    ## l) tidy effect estimates
      breast_ep_repro_zinf_nb1_m_tidy <- 
        tidy(breast_ep_repro_zinf_nb1_m_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = T) %>%
        filter(!(term == '(Intercept)'| term == 'sexm')) %>%
        select(c('term', 'estimate', #'std.error', 
                 'conf.low', 'conf.high', 
                 #'statistic', 
                 'p.value'))  
      
    ## m) label model
      breast_ep_repro_zinf_nb1_m_tidy$model <-
        c('EP ~ breast color')
      
    ## n) tidy breast_ep_repro_zinf_nb1_m_tidy by making a new column '95% CI'
      breast_ep_repro_zinf_nb1_m_tidy$'95% CI' = 
        paste(round(breast_ep_repro_zinf_nb1_m_tidy$conf.low, 3),
              round(breast_ep_repro_zinf_nb1_m_tidy$conf.high, 3), 
              sep=", ")
      
    ## o) tidy breast_ep_repro_zinf_nb1_m_tidy by removing extra columns and rows
      breast_ep_repro_zinf_nb1_m_tidy <- breast_ep_repro_zinf_nb1_m_tidy %>%
        rename(c(`Rate Ratio (exponentiated)` = `estimate`, 
                 `P-value` = `p.value`)) %>%
        select(c(`Rate Ratio (exponentiated)`, `95% CI`, `P-value`, `model`))
      
    ## p) move 'model' to first column
      breast_ep_repro_zinf_nb1_m_tidy <- breast_ep_repro_zinf_nb1_m_tidy %>%
        relocate(model, .before = `Rate Ratio (exponentiated)`) 
      
      
  ### 3.7 Models of natural tail streamer on reproduction - male models
    ## a) glmmTMB model for effect of natural tail streamer on the 
      # number of nestlings fledged from nest attempts 1 and 2
      TS_repro_zinf_nb1_m_mod <- glmmTMB(tot.repro ~
                            scale(rt.intact.avg.ts) +   
                        # proportional covariance matrix based on relatedness
                            propto(0 + sample.id | 1, relate_mat_plum_repro_m),
                            data = snps_plum_repro_m_df,
                            ziformula=~1,
                              family = nbinom1)
      
    ## b) model diagnostics
      performance::check_overdispersion(TS_repro_zinf_nb1_m_mod) 
      performance::check_zeroinflation(TS_repro_zinf_nb1_m_mod)
      
    ## c) model summary
      summary(TS_repro_zinf_nb1_m_mod)
      
    ## d) tidy effect estimates
      TS_repro_zinf_nb1_m_tidy <- 
        tidy(TS_repro_zinf_nb1_m_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = T) %>%
        filter(!(term == '(Intercept)'| term == 'sexm')) %>%
        select(c('term', 'estimate', #'std.error', 
                 'conf.low', 'conf.high', 
                 #'statistic', 
                 'p.value'))  
      
    ## e) label model
      TS_repro_zinf_nb1_m_tidy$model <-
        c('total reproduction ~ TS')
      
    ## f) tidy TS_repro_zinf_nb1_m_tidy by making a new column '95% CI'
      TS_repro_zinf_nb1_m_tidy$'95% CI' = 
        paste(round(TS_repro_zinf_nb1_m_tidy$conf.low, 3),
              round(TS_repro_zinf_nb1_m_tidy$conf.high, 3), 
              sep=", ")
      
    ## g) tidy TS_repro_zinf_nb1_m_tidy by removing extra columns and rows
      TS_repro_zinf_nb1_m_tidy <- TS_repro_zinf_nb1_m_tidy %>%
        rename(c(`Rate Ratio (exponentiated)` = `estimate`, 
                 `P-value` = `p.value`)) %>%
        select(c(`Rate Ratio (exponentiated)`, `95% CI`, `P-value`, `model`))
      
    ## h) move 'model' to first column
      TS_repro_zinf_nb1_m_tidy <- TS_repro_zinf_nb1_m_tidy %>%
        relocate(model, .before = `Rate Ratio (exponentiated)`) 
      
    ## i) glmmTMB model for effect of natural tail streamer on the 
      # number of EP nestlings fledged from nest attempts 1 and 2
      TS_ep_repro_zinf_nb1_m_mod <-  glmmTMB(ep.repro ~
                            scale(rt.intact.avg.ts) +
                        # proportional covariance matrix based on relatedness
                            propto(0 + sample.id | 1, relate_mat_plum_repro_m),
                            data = snps_plum_repro_m_df,
                            ziformula=~1,
                              family = nbinom1)
      
    ## j) model diagnostics
      performance::check_overdispersion(TS_ep_repro_zinf_nb1_m_mod) 
      performance::check_zeroinflation(TS_ep_repro_zinf_nb1_m_mod)
      
    ## k) model summary
      summary(TS_ep_repro_zinf_nb1_m_mod)
      
    ## l) tidy effect estimates
      TS_ep_repro_zinf_nb1_m_tidy <- 
        tidy(TS_ep_repro_zinf_nb1_m_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = T) %>%
        filter(!(term == '(Intercept)'| term == 'sexm')) %>%
        select(c('term', 'estimate', #'std.error', 
                 'conf.low', 'conf.high', 
                 #'statistic', 
                 'p.value'))  
      
    ## m) label model
      TS_ep_repro_zinf_nb1_m_tidy$model <-
        c('EP ~ TS')
      
    ## n) tidy TS_ep_repro_zinf_nb1_m_tidy by making a new column '95% CI'
      TS_ep_repro_zinf_nb1_m_tidy$'95% CI' = 
        paste(round(TS_ep_repro_zinf_nb1_m_tidy$conf.low, 3),
              round(TS_ep_repro_zinf_nb1_m_tidy$conf.high, 3), 
              sep=", ")
      
    ## o) tidy TS_ep_repro_zinf_nb1_m_tidy by removing extra columns and rows
      TS_ep_repro_zinf_nb1_m_tidy <- TS_ep_repro_zinf_nb1_m_tidy %>%
        rename(c(`Rate Ratio (exponentiated)` = `estimate`, 
                 `P-value` = `p.value`)) %>%
        select(c(`Rate Ratio (exponentiated)`, `95% CI`, `P-value`, `model`))
      
    ## p) move 'model' to first column
      TS_ep_repro_zinf_nb1_m_tidy <- TS_ep_repro_zinf_nb1_m_tidy %>%
        relocate(model, .before = `Rate Ratio (exponentiated)`) 
      
      
  ### 3.8 Make formatted tables of the associations between natural
      # breast brightness and tail streamer with reproduction - males
    ## a) row bind breast_repro_zinf_nb1_m_tidy, 
      # breast_ep_repro_zinf_nb1_m_tidy, TS_repro_zinf_nb1_m_tidy, and 
      # TS_ep_repro_zinf_nb1_m_tidy,
      male_plum_repro_sum <- bind_rows(breast_repro_zinf_nb1_m_tidy,
                                      breast_ep_repro_zinf_nb1_m_tidy,
                                      TS_repro_zinf_nb1_m_tidy, 
                                      TS_ep_repro_zinf_nb1_m_tidy) 
      
    ## b) format table creating flextable object 
      male_plum_repro_sum_frmt <- nice_table(male_plum_repro_sum, 
                                            separate.header = T)
      
    ## c) preview the table as a word doc    
      #print(male_plum_repro_sum_frmt, preview = 'docx')
      
    ## d) save male_plum_repro_sum_frmt as a word doc
      flextable::save_as_docx(male_plum_repro_sum_frmt, 
                              
                              path = here('output/male_plum_repro_sum_frmt.docx'))  
      




###############################################################################
##############                  4. Export data                   ##############
###############################################################################

  ### 4.1 Export data to an RData file
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.

    ## a) Save and export natural plumage reproduction models for graphing
      save(file = here('data/5_nat_plum_repro_mods.RData'),
           list = c('female_plum_repro_sum', 'male_plum_repro_sum'))

  

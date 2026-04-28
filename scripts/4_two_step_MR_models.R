################################################################################
#############                                                      #############
#############     Assessing the causal effect of plumage traits    #############
#############     on reproduction using Mendelian randomization    #############   
#############                                                      #############
#############      4. Two step Mendelian randomization models      #############
#############                                                      #############
#############                   By: Zach Laubach                   #############
#############                 created: 9 Sept 2025                 #############
#############             last updated: 27 April 2026              #############
################################################################################


  ### PURPOSE: Tidy the 
  
  # Code Blocks
    # 1. Configure work space
    # 2. Import data 
    # 3. Step 1: Both sexes snp effects on traits
    # 4. Model fit/diagnostic checks 
    # 5. Step 2: Female fit values on reproduction
    # 6. Step 2: Male fit values on reproduction
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
      # load GMMAT package
        library ('GMMAT')
      
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
    ## a) load snps, phenotype, and reproduction data
      load(here('data/2-4_snps_pheno_repro_df.RData'))
      
    ## b) Load relatedness data
      load(here('data/3-4_relate_df.RData'))
      
    

###############################################################################
##############    3. Step 1: Both sexes snp effects on traits    ##############
###############################################################################           
      
  # For Step 1, we will regress plumage traits on each candidate snp - the 
      # genetic instrument (e.g., plumage ~ lead SNPs for each plumage locus) 
      # to get predicted plumage values for use in Step 2.
    
  
  ### 3.1 KITLG effect on breast brightness - both sexes model
    ## a) glmmkin model for effect of kitlg on breast brightness
      kitlg_plum_mod <- glmmkin(breast.avg.bright ~ kitlg, 
                                  data = snps_plum_repro_df, 
                                  kins = relate_mat_plum_snps, 
                                  id = 'sample.id',
                                  family = gaussian(link = 'identity'))

    ## b) model summary  
      summary(kitlg_plum_mod)
      
    ## c) extract coefficients
      kitlg_plum_tidy <- data.frame(
        term = names(kitlg_plum_mod$coefficients),
        estimate = kitlg_plum_mod$coefficients,
        std.error = sqrt(diag(kitlg_plum_mod$cov)),
        stringsAsFactors = FALSE
      )
      
    ## d) compute Wald stats and P-values
      kitlg_plum_tidy <- kitlg_plum_tidy %>%
        mutate(
          statistic = estimate / std.error,
          p.value = 2 * (1 - pnorm(abs(statistic)))
        )
     
    ## e) calculate Wald type 95% CI and tidy data framer 
      kitlg_plum_tidy <- kitlg_plum_tidy %>%
        mutate(conf.low = estimate - 1.96 * std.error) %>%
        mutate(conf.high = estimate + 1.96 * std.error) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'conf.low', 'conf.high'))
      
    ## f) label model
      kitlg_plum_tidy$model <- c('KITLG')    
      
    ## g) Get fitted values (response scale)
      #kitlg_plum_pred_val <- fitted.values(kitlg_plum_mod)
      kitlg_plum_pred_val <- predict(kitlg_plum_mod, type = "response")
      
    ## h) Combine fitted values with original data
      snps_plum_repro_df$kitlg.plum.pred.val <- kitlg_plum_pred_val
      
      
  ### 3.2 SLC45A2 effect on breast brightness - both sexes model
    ## a) glmmkin model for effect of slc45a2 on breast brightness
      slc45a2_plum_mod <- glmmkin(breast.avg.bright ~ slc45a2, 
                                    data = snps_plum_repro_df, 
                                  kins = relate_mat_plum_snps, 
                                  id = 'sample.id',
                                  family = gaussian(link = 'identity'))
      
    ## b) summarize model  
      summary(slc45a2_plum_mod)
      
    ## c) extract coefficients
      slc45a2_plum_tidy <- data.frame(
        term = names(slc45a2_plum_mod$coefficients),
        estimate = slc45a2_plum_mod$coefficients,
        std.error = sqrt(diag(slc45a2_plum_mod$cov)),
        stringsAsFactors = FALSE
      )
      
    ## d) compute Wald stats and P-values
      slc45a2_plum_tidy <- slc45a2_plum_tidy %>%
        mutate(
          statistic = estimate / std.error,
          p.value = 2 * (1 - pnorm(abs(statistic)))
        )
      
    ## e) calculate Wald type 95% CI and tidy data framer 
      slc45a2_plum_tidy <- slc45a2_plum_tidy %>%
        mutate(conf.low = estimate - 1.96 * std.error) %>%
        mutate(conf.high = estimate + 1.96 * std.error) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'conf.low', 'conf.high'))
      
    ## f) label model
      slc45a2_plum_tidy$model <- c('SLC45A2')  
      
    ## g) Get fitted values (response scale)Get fitted values (response scale)
      slc45a2_plum_pred_val <- predict(slc45a2_plum_mod, type = "response")
      
    ## h) Combine with original data
      snps_plum_repro_df$slc45a2.plum.pred.val <- slc45a2_plum_pred_val
      
      
  ### 3.3 ROR2 effect on breast brightness - both sexes model
    ## a) glmmkin model for effect of ror2 on breast brightness
      ror2_plum_mod <- glmmkin(breast.avg.bright ~ ror2, 
                                  data = snps_plum_repro_df, 
                               kins = relate_mat_plum_snps, 
                               id = 'sample.id',
                               family = gaussian(link = 'identity'))
      
    ## b) summarize model  
      summary(ror2_plum_mod)
      
    ## c) extract coefficients
      ror2_plum_tidy <- data.frame(
        term = names(ror2_plum_mod$coefficients),
        estimate = ror2_plum_mod$coefficients,
        std.error = sqrt(diag(ror2_plum_mod$cov)),
        stringsAsFactors = FALSE
      )
      
    ## d) compute Wald stats and P-values
      ror2_plum_tidy <- ror2_plum_tidy %>%
        mutate(
          statistic = estimate / std.error,
          p.value = 2 * (1 - pnorm(abs(statistic)))
        )
      
    ## e) calculate Wald type 95% CI and tidy data framer 
      ror2_plum_tidy <- ror2_plum_tidy %>%
        mutate(conf.low = estimate - 1.96 * std.error) %>%
        mutate(conf.high = estimate + 1.96 * std.error) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'conf.low', 'conf.high'))
      
    ## f) label model
      ror2_plum_tidy$model <- c('ROR2')  
      
    ## g) Get fitted values (response scale)
     ror2_plum_pred_val <- predict(ror2_plum_mod, type = "response")
      
    ## h) Combine with original data
      snps_plum_repro_df$ror2.plum.pred.val <- ror2_plum_pred_val
      
      
  ### 3.4 PDE1C effect on tail streamer length  - both sexes model
    ## a) glmmkin model for effect of pde1c on tail streamer length
      pde1c_ts_mod <- glmmkin(max.avg.ts ~ pde1c, data = snps_TS_repro_df, 
                              kins = relate_mat_TS_snps, 
                              id = 'sample.id',
                              family = gaussian(link = 'identity'))
      
    ## b) summarize model 
      summary(pde1c_ts_mod)
      
    ## c) extract coefficients
      pde1c_ts_tidy <- data.frame(
        term = names(pde1c_ts_mod$coefficients),
        estimate = pde1c_ts_mod$coefficients,
        std.error = sqrt(diag(pde1c_ts_mod$cov)),
        stringsAsFactors = FALSE
      )
      
    ## d) compute Wald stats and P-values
      pde1c_ts_tidy <- pde1c_ts_tidy %>%
        mutate(
          statistic = estimate / std.error,
          p.value = 2 * (1 - pnorm(abs(statistic)))
        )
      
    ## e) calculate Wald type 95% CI and tidy data framer 
      pde1c_ts_tidy <- pde1c_ts_tidy %>%
        mutate(conf.low = estimate - 1.96 * std.error) %>%
        mutate(conf.high = estimate + 1.96 * std.error) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'conf.low', 'conf.high'))
      
      ## f) label model
      pde1c_ts_tidy$model <- c('PDE1C')  
      
    ## c) Get fitted values (response scale)
      pde1c_ts_pred_val <- predict(pde1c_ts_mod, type = "response")
      
    ## d) Combine with original data
      snps_TS_repro_df$pde1c.ts.pred.val <- pde1c_ts_pred_val
      
       
  ### 3.5 Are candidate snps associated with  plumage traits? 
    # a) Combine snp instrument regression estimates into a tidy table
      snps_trait_smmry <- 
        bind_rows(kitlg_plum_tidy,
                  slc45a2_plum_tidy,
                  ror2_plum_tidy,
                  pde1c_ts_tidy)
      
    ## b) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of 'model' variable 
      snps_trait_smmry <- 
        transform(snps_trait_smmry, 
                  model = factor(model,
                                 levels = c('KITLG',
                                            'SLC45A2',
                                            'ROR2',
                                            'PDE1C')))  
      
      # change order of columns
      snps_trait_smmry <- snps_trait_smmry %>%
        relocate(model, .before = term)
      
    ## c) format table creating flextable object 
      snps_trait_frmt <- 
        nice_table(snps_trait_smmry, separate.header = T)
      
    ## d) preview the table as a word doc    
      #print(snps_trait_frmt, preview = 'docx')
      
    ## e) save snps_trait_frmt as a word doc
      flextable::save_as_docx(snps_trait_frmt, 
                              path = here('output/snps_trait_frmt.docx')) 
      
      
  ### 3.6 Tidy predicted breast brightness values from Step 1 for use in Step 2
    ## a) select the predicted values
      snps_plum_repro_df <- snps_plum_repro_df %>%
        select(c(sample.id, kitlg.plum.pred.val, slc45a2.plum.pred.val,
                 ror2.plum.pred.val))
    
    ## b) Left join the predicted values for breast plumage brightness
        # to snps_plum_repro_f_df
      snps_plum_repro_f_df <- snps_plum_repro_f_df  %>%
        left_join(snps_plum_repro_df,
                  by = c('sample.id' = 'sample.id'), 
                  copy = F) 
    
    ## c) Left join the predicted values for breast plumage brightness
      # to snps_plum_repro_m_df
      snps_plum_repro_m_df <- snps_plum_repro_m_df  %>%
        left_join(snps_plum_repro_df,
                  by = c('sample.id' = 'sample.id'), 
                  copy = F) 

      
  ### 3.7 Tidy predicted tail streamer values from Step 1 for use in Step 2
    ## a) select the predicted values
      snps_TS_repro_df <- snps_TS_repro_df %>%
        select(c(sample.id, pde1c.ts.pred.val))
      
    ## b) Left join the predicted values for breast plumage brightness
      # to snps_TS_repro_f_df
      snps_TS_repro_f_df <- snps_TS_repro_f_df  %>%
        left_join(snps_TS_repro_df,
                  by = c('sample.id' = 'sample.id'), 
                  copy = F) 
      
    ## c) Left join the predicted values for breast plumage brightness
      # to snps_TS_repro_m_df
      snps_TS_repro_m_df <- snps_TS_repro_m_df  %>%
        left_join(snps_TS_repro_df,
                  by = c('sample.id' = 'sample.id'), 
                  copy = F) 
      
      
      
###############################################################################
##############          4. Model fit/diagnostic checks           ##############
###############################################################################         
      
      
      # kitlg_fit_repro_mod <- glmmTMB(NCalls~(FT+ArrivalTime)*SexParent+ 
      #                                  offset(log(BroodSize))+(1|Nest), 
      #                                data=Owls, ziformula=~1, family=poisson)
      
  ### 4.1 Model fitting to identify the appropriate distribution for count
        # data
    ## a) Fit model with Poisson distribution
      kitlg_fit_repro_f_poiss_mod <- glmmTMB(tot.repro ~ kitlg.plum.pred.val, 
                                    data = snps_plum_repro_f_df, 
                                    family = poisson)
      
    ## b) udpate the Poisson model to include zero-inflation
      kitlg_fit_repro_f_zinf_poiss_mod <- update(kitlg_fit_repro_f_poiss_mod,
                                                 ziformula=~1)

    ## c) udpate the Poisson model to fit a negative binomial 2 distribution
        # where variance = µ(1 + µ/k)
      kitlg_fit_repro_f_nb2_mod <- update(kitlg_fit_repro_f_poiss_mod,
                                                 family = nbinom2)
      
    ## d) udpate the negative binomial 2 model to include zero-inflation
      kitlg_fit_repro_f_zinf_nb2_mod <- update(kitlg_fit_repro_f_nb2_mod,
                                               ziformula=~1)
      
    ## e) udpate the Poisson model to fit a negative binomial 1 distribution
      # where variance = ϕµ
      kitlg_fit_repro_f_nb1_mod <- update(kitlg_fit_repro_f_poiss_mod,
                                          family = nbinom1)
      
    ## f) udpate the negative binomial 1 model to include zero-inflation
      kitlg_fit_repro_f_zinf_nb1_mod <- update(kitlg_fit_repro_f_nb1_mod,
                                               ziformula=~1)

        
    ## g) assess delta AIC using bbmle to determine best fitting model   
      AICtab(kitlg_fit_repro_f_poiss_mod, kitlg_fit_repro_f_zinf_poiss_mod,
             kitlg_fit_repro_f_nb2_mod, kitlg_fit_repro_f_zinf_nb2_mod, 
             kitlg_fit_repro_f_nb1_mod, kitlg_fit_repro_f_zinf_nb1_mod)
    
        
  ### 4.2 check diagnostics of best fitting model using 'DHARMa' and 
      # 'performance' packages
    ## a)  simulate residuals
      resid_best_fit_mod <- simulateResiduals(fittedModel = kitlg_fit_repro_f_zinf_poiss_mod, 
                               n = 250)
   
    ## b) plot residuals
      plot(resid_best_fit_mod)
      
    ## c) assess over dispersion
      testDispersion(resid_best_fit_mod)
      performance::check_overdispersion(kitlg_fit_repro_f_zinf_poiss_mod) 
      
    ## d) assess zero inflation
      testZeroInflation(resid_best_fit_mod)
      performance::check_zeroinflation(kitlg_fit_repro_f_zinf_poiss_mod)

      
      
###############################################################################
##############   5. Step 2: Female fit values on reproduction    ##############
###############################################################################              
    
  # For Step 2: Female models, we will regress measures of reproductive
      # on the predicted values of each plumage trait
      # (e.g., reproductive success ~ predicted trait values) 
      # to get the unconfounded effect of traits on reproductive success.
      

  ### 5.1 KITLG predicted breast brightness on reproduction - female models
    ## a) glmmTMB model for effect of predicted breast brightness on the 
        # number of nestlings fledged from nest attempts 1 and 2
      kitlg_fit_repro_f_zinf_poiss_mod <- glmmTMB(tot.repro ~ 
                                              kitlg.plum.pred.val, 
                                              data = snps_plum_repro_f_df,
                                              ziformula=~1,
                                              family = poisson)
      
    ## b) model diagnostics
      performance::check_overdispersion(kitlg_fit_repro_f_zinf_poiss_mod) 
      performance::check_zeroinflation(kitlg_fit_repro_f_zinf_poiss_mod)
    
    ## c) model summary
      summary(kitlg_fit_repro_f_zinf_poiss_mod)
      
    ## d) tidy effect estimates
      kitlg_fit_repro_f_zinf_poiss_tidy <- 
        tidy(kitlg_fit_repro_f_zinf_poiss_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = T) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'conf.low', 'conf.high'))
      
    ## e) label model
      kitlg_fit_repro_f_zinf_poiss_tidy$model <- 
        c('KITLG')  
      
  ### 5.2 SLC45A2 predicted breast brightness on reproduction - female models
    ## a) glmmTMB model for effect of predicted breast brightness on the 
      # number of nestlings fledged from nest attempts 1 and 2
      slc45a2_fit_repro_f_zinf_poiss_mod <- glmmTMB(tot.repro ~ 
                                                slc45a2.plum.pred.val, 
                                                data = snps_plum_repro_f_df,
                                                ziformula=~1,
                                                family = poisson)
      
    ## b) model diagnostics
      performance::check_overdispersion(slc45a2_fit_repro_f_zinf_poiss_mod) 
      performance::check_zeroinflation(slc45a2_fit_repro_f_zinf_poiss_mod)
      
    ## c) model summary
      summary(slc45a2_fit_repro_f_zinf_poiss_mod)
      
    ## d) tidy effect estimates
      slc45a2_fit_repro_f_zinf_poiss_tidy <- 
        tidy(slc45a2_fit_repro_f_zinf_poiss_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = T) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'conf.low', 'conf.high'))
      
    ## e) label model
      slc45a2_fit_repro_f_zinf_poiss_tidy$model <- 
        c('SLC45A2')  
      
      
  ### 5.3 ROR2 predicted breast brightness on reproduction - female models
    ## a) glmmTMB model for effect of predicted breast brightness on the 
      # number of nestlings fledged from nest attempts 1 and 2
      ror2_fit_repro_f_zinf_poiss_mod <- glmmTMB(tot.repro ~ 
                                                  ror2.plum.pred.val, 
                                                  data = snps_plum_repro_f_df,
                                                  ziformula=~1,
                                                  family = poisson)
      
    ## b) model diagnostics
      performance::check_overdispersion(ror2_fit_repro_f_zinf_poiss_mod) 
      performance::check_zeroinflation(ror2_fit_repro_f_zinf_poiss_mod)
      
    ## c) model summary
      summary(ror2_fit_repro_f_zinf_poiss_mod)
      
    ## d) tidy effect estimates
      ror2_fit_repro_f_zinf_poiss_tidy <- 
        tidy(ror2_fit_repro_f_zinf_poiss_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = T) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'conf.low', 'conf.high'))
      
    ## e) label model
      ror2_fit_repro_f_zinf_poiss_tidy$model <- 
        c('ROR2') 
      
      
  ### 5.4 PDE1C predicted tail streamer length on reproduction - female models
    ## a) glmmTMB model for effect of predicted TS length on the 
      # number of nestlings fledged from nest attempts 1 and 2
      pde1c_fit_repro_f_zinf_poiss_mod <- glmmTMB(tot.repro ~ 
                                              pde1c.ts.pred.val, 
                                              data = snps_TS_repro_f_df,
                                              ziformula=~1,
                                              family = poisson)
      
    ## b) model diagnostics
      performance::check_overdispersion(pde1c_fit_repro_f_zinf_poiss_mod) 
      performance::check_zeroinflation(pde1c_fit_repro_f_zinf_poiss_mod)
      
    ## c) model summary
      summary(pde1c_fit_repro_f_zinf_poiss_mod)
      
    ## d) tidy effect estimates
      pde1c_fit_repro_f_zinf_poiss_tidy <- 
        tidy(pde1c_fit_repro_f_zinf_poiss_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = T) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'conf.low', 'conf.high'))
      
    ## e) label model
      pde1c_fit_repro_f_zinf_poiss_tidy$model <- 
        c('PDE1C') 
      
      
  ### 5.5 Do female plumage traits cause differences in reproduction? 
    # a) Combine Step 2: Female regression estimates into a tidy table
      trait_tot_repro_f_smmry <- 
        bind_rows(kitlg_fit_repro_f_zinf_poiss_tidy,
                  slc45a2_fit_repro_f_zinf_poiss_tidy,
                  ror2_fit_repro_f_zinf_poiss_tidy,
                  pde1c_fit_repro_f_zinf_poiss_tidy)
      
    ## b) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of 'model' variable 
      trait_tot_repro_f_smmry <- 
        transform(trait_tot_repro_f_smmry, 
                  model = factor(model,
                                 levels = c('KITLG',
                                            'SLC45A2',
                                            'ROR2',
                                            'PDE1C')))  
      
      # change order of columns
      trait_tot_repro_f_smmry <- trait_tot_repro_f_smmry %>%
        relocate(model, .before = term)
      
    ## c) format table creating flextable object 
      trait_tot_repro_f_frmt <- 
        nice_table(trait_tot_repro_f_smmry, separate.header = T)
      
    ## d) preview the table as a word doc    
      #print(trait_tot_repro_f_frmt, preview = 'docx')
      
    ## e) save trait_tot_repro_f_frmt as a word doc
      flextable::save_as_docx(trait_tot_repro_f_frmt, 
                              path = here('output/trait_tot_repro_f_frmt.docx')) 
      
     

###############################################################################
##############    6. Step 2: Male fit values on reproduction     ##############
###############################################################################         

  ### 6.1 KITLG predicted breast brightness on reproduction - male models
    ## a) glmmTMB model for effect of predicted breast brightness on male's
        # total paternity (WP + EP) from nest attempts 1 and 2
      kitlg_fit_repro_m_zinf_poiss_mod <- glmmTMB(tot.repro ~ 
                                                kitlg.plum.pred.val, 
                                                data = snps_plum_repro_m_df,
                                                ziformula=~1,
                                                family = poisson)
      
    ## b) model diagnostics
      performance::check_overdispersion(kitlg_fit_repro_m_zinf_poiss_mod) 
      performance::check_zeroinflation(kitlg_fit_repro_m_zinf_poiss_mod)
      
    ## c) model summary
      summary(kitlg_fit_repro_m_zinf_poiss_mod)
      
    ## d) tidy effect estimates
      kitlg_fit_repro_m_zinf_poiss_tidy <- 
        tidy(kitlg_fit_repro_m_zinf_poiss_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = T) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'conf.low', 'conf.high'))
      
    ## e) label model
      kitlg_fit_repro_m_zinf_poiss_tidy$model <- 
        c('KITLG') 
      
      
  ### 6.2 SLC45A2 predicted breast brightness on reproduction - male models
    ## a) glmmTMB model for effect of predicted breast brightness on male's
      # total paternity (WP + EP) from nest attempts 1 and 2 
      slc45a2_fit_repro_m_zinf_poiss_mod <- glmmTMB(tot.repro ~ 
                                                slc45a2.plum.pred.val, 
                                                data = snps_plum_repro_m_df,
                                                ziformula=~1,
                                                family = poisson)
      
    ## b) model diagnostics
      performance::check_overdispersion(slc45a2_fit_repro_m_zinf_poiss_mod) 
      performance::check_zeroinflation(slc45a2_fit_repro_m_zinf_poiss_mod)
      
    ## c) model summary
      summary(slc45a2_fit_repro_m_zinf_poiss_mod)
      
    ## d) tidy effect estimates
      slc45a2_fit_repro_m_zinf_poiss_tidy <- 
        tidy(slc45a2_fit_repro_m_zinf_poiss_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = T) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'conf.low', 'conf.high'))
      
    ## e) label model
      slc45a2_fit_repro_m_zinf_poiss_tidy$model <- 
        c('SLC45A2') 
      
      
  ### 6.3 ROR2 predicted breast brightness on reproduction - male models
    ## a) glmmTMB model for effect of predicted breast brightness on male's
      # total paternity (WP + EP) from nest attempts 1 and 2 
      ror2_fit_repro_m_zinf_poiss_mod <- glmmTMB(tot.repro ~ 
                                              ror2.plum.pred.val, 
                                              data = snps_plum_repro_m_df,
                                              ziformula=~1,
                                              family = poisson)
      
    ## b) model diagnostics
      performance::check_overdispersion(ror2_fit_repro_m_zinf_poiss_mod) 
      performance::check_zeroinflation(ror2_fit_repro_m_zinf_poiss_mod)
      
    ## c) model summary
      summary(ror2_fit_repro_m_zinf_poiss_mod)
      
    ## d) tidy effect estimates
      ror2_fit_repro_m_zinf_poiss_tidy <- 
        tidy(ror2_fit_repro_m_zinf_poiss_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = T) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'conf.low', 'conf.high'))
      
    ## e) label model
      ror2_fit_repro_m_zinf_poiss_tidy$model <- 
        c('ROR2') 
      
      
  ### 6.4 PDE1C predicted tail streamer length on reproduction - male models
    ## a) glmmTMB model for effect of predicted TS on male's
      # total paternity (WP + EP) from nest attempts 1 and 2 
      pde1c_fit_repro_m_zinf_poiss_mod <- glmmTMB(tot.repro ~ 
                                              pde1c.ts.pred.val, 
                                              data = snps_TS_repro_m_df,
                                              ziformula=~1,
                                              family = poisson)
      
    ## b) model diagnostics
      performance::check_overdispersion(pde1c_fit_repro_m_zinf_poiss_mod) 
      performance::check_zeroinflation(pde1c_fit_repro_m_zinf_poiss_mod)
     
    ## c) model summary
      summary(pde1c_fit_repro_m_zinf_poiss_mod)
     
    ## d) tidy effect estimates
      pde1c_fit_repro_m_zinf_poiss_tidy <-
        tidy(pde1c_fit_repro_m_zinf_poiss_mod, 
            effects = 'fixed',
            conf.int = T,
            conf.method = 'Wald',
            conf.level = 0.95, 
            exponentiate = T) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'conf.low', 'conf.high'))
     
    ## e) label model
      pde1c_fit_repro_m_zinf_poiss_tidy$model <- 
        c('PDE1C')  
      
      
  ### 6.5 Do male plumage traits cause differences in reproduction? 
    # a) Combine Step 2: Male regression estimates into a tidy table
      trait_tot_repro_m_smmry <- 
        bind_rows(kitlg_fit_repro_m_zinf_poiss_tidy,
                  slc45a2_fit_repro_m_zinf_poiss_tidy,
                  ror2_fit_repro_m_zinf_poiss_tidy,
                  pde1c_fit_repro_m_zinf_poiss_tidy)
      
    ## b) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of 'model' variable 
      trait_tot_repro_m_smmry <- 
        transform(trait_tot_repro_m_smmry, 
                  model = factor(model,
                                 levels = c('KITLG',
                                            'SLC45A2',
                                            'ROR2',
                                            'PDE1C')))  
      
      # change order of columns
      trait_tot_repro_m_smmry <- trait_tot_repro_m_smmry %>%
        relocate(model, .before = term)
      
    ## c) format table creating flextable object 
      trait_tot_repro_m_frmt <- 
        nice_table(trait_tot_repro_m_smmry, separate.header = T)
      
    ## d) preview the table as a word doc    
      #print(trait_tot_repro_m_frmt, preview = 'docx')
      
    ## e) save trait_tot_repro_m_frmt as a word doc
      flextable::save_as_docx(trait_tot_repro_m_frmt, 
                        path = here('output/trait_tot_repro_m_frmt.docx')) 
      
      
      
    
      

      
      
      ################################################################################      
      #       
      #       
      #       
      #   ### 3.1 KITLG breast brightness association - Female model
      #       kitlg_plum_f_mod <- glmmkin(breast.avg.bright ~ kitlg, 
      #                                 data = snps_plum_fecund_f_df, 
      #                         kins = relate_mat_plum_fecund_f, id = 'sample.id',
      #                         family = gaussian(link = 'identity'))
      #  
      #       summary(kitlg_plum_f_mod)
      # 
      #       
      #       #Get fitted values (response scale)
      #       #kitlg_plum_pred_val <- fitted.values(kitlg_plum_f_mod)
      #       kitlg_plum_f_pred_val <- predict(kitlg_plum_f_mod, type = "response")
      #       
      #       # Combine with original data
      #       snps_plum_fecund_f_df$kitlg.plum.f.pred.val <- kitlg_plum_f_pred_val
      #       
      #       
      #   ### 3.2 SLC45A2 breast brightness association - Female model
      #       slc45a2_plum_f_mod <- glmmkin(breast.avg.bright ~ slc45a2, 
      #                                   data = snps_plum_fecund_f_df, 
      #                                   kins = relate_mat_plum_fecund_f, id = 'sample.id',
      #                                   family = gaussian(link = 'identity'))
      #       
      #       
      #       summary(slc45a2_plum_f_mod)
      #       
      #       
      #       #Get fitted values (response scale)
      #       slc45a2_plum_f_pred_val <- predict(slc45a2_plum_f_mod, type = "response")
      #       
      #       # Combine with original data
      #       snps_plum_fecund_f_df$slc45a2.plum.f.pred.val <- slc45a2_plum_f_pred_val
      #       
      #       
      #     ### 3.3 ROR2 breast brightness association - Female model
      #       ror2_plum_f_mod <- glmmkin(breast.avg.bright ~ ror2, 
      #                                     data = snps_plum_fecund_f_df, 
      #                                     kins = relate_mat_plum_fecund_f, id = 'sample.id',
      #                                     family = gaussian(link = 'identity'))
      #       
      #       summary(ror2_plum_f_mod)
      #       
      #       #Get fitted values (response scale)
      #       ror2_plum_f_pred_val <- predict(ror2_plum_f_mod, type = "response")
      #       
      #       # Combine with original data
      #       snps_plum_fecund_f_df$ror2.plum.f.pred.val <- ror2_plum_f_pred_val
      #       
      #       
      # ################################################################################  
      #       
      #       
      #   ### 3.1 KITLG breast brightness association - male model
      #       kitlg_plum_m_mod <- glmmkin(breast.avg.bright ~ kitlg, 
      #                                   data = snps_plum_pat_m_df, 
      #                                   kins = relate_mat_plum_pat_m, id = 'sample.id',
      #                                   family = gaussian(link = 'identity'))
      #       
      #       summary(kitlg_plum_m_mod)
      #       
      #       
      #       #Get fitted values (response scale)
      #       kitlg_plum_m_pred_val <- predict(kitlg_plum_m_mod, type = "response")
      #       
      #       # Combine with original data
      #       snps_plum_pat_m_df$kitlg.plum.m.pred.val <- kitlg_plum_m_pred_val
      #       
      #       
      #       ### 3.2 SLC45A2 breast brightness association - male model
      #       slc45a2_plum_m_mod <- glmmkin(breast.avg.bright ~ slc45a2, 
      #                                     data = snps_plum_pat_m_df, 
      #                                     kins = relate_mat_plum_pat_m, id = 'sample.id',
      #                                     family = gaussian(link = 'identity'))
      #       
      #       
      #       summary(slc45a2_plum_m_mod)
      #       
      #       
      #       #Get fitted values (response scale)
      #       slc45a2_plum_m_pred_val <- predict(slc45a2_plum_m_mod, type = "response")
      #       
      #       # Combine with original data
      #       snps_plum_pat_m_df$slc45a2.plum.m.pred.val <- slc45a2_plum_m_pred_val
      #       
      #       
      #       ### 3.3 ROR2 breast brightness association - male model
      #       ror2_plum_m_mod <- glmmkin(breast.avg.bright ~ ror2, 
      #                                  data = snps_plum_pat_m_df, 
      #                                  kins = relate_mat_plum_pat_m, id = 'sample.id',
      #                                  family = gaussian(link = 'identity'))
      #       
      #       summary(ror2_plum_m_mod)
      #       
      #       #Get fitted values (response scale)
      #       ror2_plum_m_pred_val <- predict(ror2_plum_m_mod, type = "response")
      #       
      #       # Combine with original data
      #       snps_plum_pat_m_df$ror2.plum.m.pred.val <- ror2_plum_m_pred_val          
      
      
              
  # TOY EXAMPLE    
      
      data("example")
      pheno <- example$pheno
      pheno2 <- example$pheno2
      GRM <- example$GRM
      
      # set up toy data
      
      
      
      pheno <- pheno %>%
        rename(c('repro' = 'age',
                 'gene' = 'sex')) %>%
        select(-disease)
      
      pheno$repro <- pheno$repro / 10
      
      model0 <- glmmkin(trait ~ gene, data = pheno, kins = GRM, id = "id",
                        family = gaussian(link = "identity"), )
      
      
      

###############################################################################
##############                  4. Export data                   ##############
###############################################################################
      
  ### 4.1 Export data to an RData file 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      
    ## a) Save and export relatedness data
      save(file = here('data/3_relatedness_data.RData'), 
           list = c('relate_mat', 'relate_coef_df'))
      
      

    

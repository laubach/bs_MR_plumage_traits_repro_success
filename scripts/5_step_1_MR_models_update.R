################################################################################
#############                                                      #############
#############     Assessing the causal effect of plumage traits    #############
#############     on reproduction using Mendelian randomization    #############   
#############                                                      #############
#############       5. Step 1: Mendelian randomization models      #############
#############                                                      #############
#############                   By: Zach Laubach                   #############
#############                 created: 9 Sept 2025                 #############
#############               last updated: 20 Oct 2026              #############
################################################################################


  ### PURPOSE: Determining the genetic instrument
  
  # Code Blocks
    # 1. Configure work space
    # 2. Import data 
    # 3. Step 1: Create female genetic instrument
    # 4. Step 1: Create male genetic instrument
    # 5. Step 1: Elastic net based female instrument
    # 6. Step 1: Elastic net based male instrument

    # 4. Sensitivity analysis
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
      
      # install.packages('rsample')
        library('rsample')
    
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
      
      # load partR2 package 
        # library('partR2')
       # library('r2glmm')
      
      # load car package
        library('car')
        library('carData')
      
      # load nlme package
        library('insight')
      
      # # load MixMod package
      #   library('MixMod')

      # install.packages('tidymodels')
        library('tidymodels')
      
      # install.packages('glmnet')
        library('glmnet')
      
      # install.packages('recipes')
        library('recipes')
      
      # install.packages('parsnip')
        library('parsnip')
  

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
        library ('broom')
      
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
      load(here('data/5_MR_data.RData'))
      
    ## b) Load relatedness data
      load(here('data/3-5_relate_df.RData'))
      
    

###############################################################################
##############    3. Step 1: Create female genetic instrument    ##############
###############################################################################           
      
  ### NOTE: Determine which snps within 10kb region on either side of 
          # top candidate snps from Schield et al. 2024 are associated
          # with the plumage traits in CO population. 
          # These snps have beta estimate that is at least as 0.5 as large
          # as the beta from the top snp from Schield et al. 2024.
          # Then assess independence of associated snps from same candidate
          # region.
      
      # The ideal model would actually include a SNPxSex interaction term so 
      # that you can calculate sex-specific genetically predicted plumage from 
      # step 1 for step 2 without losing power by stratifying. Kristen also 
      # suggested in her mouse models they would always run sex-stratified 
      # and sex-combined analyses, this would at least allow you to confirm 
      # that the associations with plumage are the same across sex (PDE1C 
      # would obviously already be different if it's monomorphic in males).
      # 1. When we condition on sex (e.g., phenotype ~ SNP + sex + relatedness) 
      # it is to control for any sex-specific differences that might impact the 
      # phenotype. Adjusting for sex helps reduce bias and improve precision of 
      # the SNP-phenotype association, especially if there are any baseline 
      # differences by sex in the phenotype. In contrast, the sex-by-SNP 
      # interaction addresses whether the effect of the SNP differs in males 
      # versus females (either in magnitude or direction).  
      # 
      # Proxy SNP: If there are other top loci that you can’t find the lead SNP 
      # in your imputed data for, then you can take a proxy SNP in LD with the 
      # first SNP (based on some other LD data) to represent that locus.
      # Independent signals in a single locus. You could conduct conditional 
      # analysis in your sample to identify these. After regressing out the 
      # lead SNP in the locus, are there any other SNPs in the region that meet 
      # some similarly high significance threshold with plumage? You could 
      # alternatively refer back to the original GWAS results to see if nearby 
      # SNPs have drastically different MAF/Effect estimates, which would 
      # suggest independence. Can you check then check their LD?
      # 3. If you are missing an original lead SNP but find a proxy (a SNP in 
      # LD with the original – we often use a r2>0.8 threshold) that is not 
      # monomorphic, then I think it is reasonable to use the proxy SNP instead. 
      # When I find proxy SNPs in my MR analyses, I replace the original SNP 
      # with the proxy SNP, too. 
      # Yes, the goal of the conditional analysis would be to test if SNPs in 
      # LD are independently associated with the outcome or if they are all 
      # tagging the same region. If they are independent, then you can include 
      # both in the PGS. I believe you do this in a stepwise manner 
      # (plumage ~ lead SNP + second SNP in LD; plumage ~ lead SNP + second SNP + third SNP in LD) 
      # until no signal remains.
     
      # discussion with Kayleigh 14/Oct - maybe stratify by sex for 
      # calculating the first stage. 
  
  # Developing instrument methodology
      # 1. Identify lead SNP from Schield et al and all SNPs +/- 10kb
          # since beyond this, SNPs not in LD and thus not a good proxy for
          # the lead SNP from the external population (hybrids)
      # 2. Assess in these windows if beta of proxy SNPs 
          # is >= 0.5 * beta of lead SNP
      # 3. Check if any of these proxy SNPs with large enough betas 
          # are polymorphic in the CO data
      # 4. Those SNPs that pass the above threshold are then assessed in 
          # regression models where the breast brightness is the outcome - 
          # identify those SNPs for which at least one allele is significant 
          # at relaxed alpha = 0.2. This is the list of CO candidate proxy SNPs 
          # for each locus
      # 5. if there is a single polymorphic SNP at a locus, carry it forward 
          # to test it in the multiple loci models - Described at the start
          # of section 4.
    
  ### 3.1 Tidy snps_plum_repro_df and snps_plum_repro_f_df
    ## a) transform snsp to factor 
      snps_plum_repro_df <- snps_plum_repro_df %>%
        mutate(across(c(kitlg.44277058:apc.46065597),
                      ~factor(., levels =  c('0', '1', '2'))))
      
    ## b) transform snsp to factor (female df)
      snps_plum_repro_f_df <- snps_plum_repro_f_df %>%
        mutate(across(c(kitlg.44277058:apc.46065597),
                      ~factor(., levels =  c('0', '1', '2'))))
      

  ### 3.2 SNPs effect on breast brightness - female
    ## a) glmmTMB model for effect of lead SNPs or closest candidate proxy 
        # SNPs based on Schiled et al - female
      snps_plum_mod_f <-  glmmTMB(breast.avg.bright ~ 
                                  kitlg.44285485 +
                                  plxnc1.46520400 +
                                  slc45a3.19332900 +
                                  bnc2.33382167 +
                                  gnaq.35901305 +
                                  ror2.44206228 +
                                  apc.46065597,
                                data = snps_plum_repro_f_df,
                                family = gaussian)
      
    ## b) model summary - female 
      summary(snps_plum_mod_f)
      
    ## c) tidy effect estimates - female
      snps_plum_mod_tidy_f <- 
        tidy(snps_plum_mod_f, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = F) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value')) 
      # NOTE: NO significant snps:
     
       
  ### 3.3 Calculate an F-statistic based on snps_plum_mod_f and the
      #  snps_plum_mod_reduced_f sums of squares
      
    ## a) run a reduced model with no snps - female
      snps_plum_mod_reduced_f <-  glmmTMB(breast.avg.bright ~
                                      #sex 
                                      + 1,
                       # proportional covariance matrix based on relatedness
                       # propto(0 + sample.id | 1, relate_mat_plum_snps),
                         data = snps_plum_repro_f_df,
                         family = gaussian)
      
    ## b) use stats::anova to compare full and reduced model
      snps_mod_comp_f <- tibble::tibble(anova(snps_plum_mod_f,
                                                  snps_plum_mod_reduced_f))
      
    ## c) add model variable
      snps_mod_comp_f$model <-  c('')
      
    ## d) move 'model' to first column
      snps_mod_comp_f <- snps_mod_comp_f %>%
        relocate(model, .before = `Df`)
      
    ## e) update values to indicate model
      snps_mod_comp_f[1,1] <- 'reduced SNPs model - female'
      snps_mod_comp_f[2,1] <- 'full SNPs model - female'
      
    ## f) Get fitted values (response scale)
      #snps_plum_pred_val <- fitted.values(snps_plum_mod_f)
      snps_plum_pred_val_f <- predict(snps_plum_mod_f,
                                    type = "response")
      
    ## g) Combine fitted values with original data
      snps_plum_repro_f_df$snps.plum.pred.val <- snps_plum_pred_val_f
    
    ## h) Determine residual and model degrees of freedom from full model
      df_resid_snps_f <- get_df(snps_plum_mod_f) # same as df.residual(model)
      df_mod_snps_f <- get_df(snps_plum_mod_f, 
                          type = "model") # same as attr(logLik(model), "df")
      
    ## i) Determine the model df from reduced model
      df_mod_snps_r_f <- get_df(snps_plum_mod_reduced_f, 
                          type = "model") # same as attr(logLik(model), "df")
      
    ## j) extract residual sums of squares from full model
      resid_snps_f <- residuals(snps_plum_mod_f, type = "response")
      resid_ss_snps_f <- sum(resid_snps_f^2)
      
    ## k) calculate the Mean Square error from full model
      ms_error_snps_f = resid_ss_snps_f/df_resid_snps_f
      
    ## l) extract full model predicted values
      y_hat_snps_f <- predict(snps_plum_mod_f,
                            type = "response")
      
    ## m) extract reduced model predicted values
      y_hat_snps_r_f <- predict(snps_plum_mod_reduced_f,
                              type = "response")
      
    ## n) extract actual response values - same for models
      y_snps_f <- model.response(model.frame(snps_plum_mod_f))
      
    ## o) Calculate Model Sum of Squares from full model
      mss_snps_f <- sum((y_hat_snps_f - mean(y_snps_f))^2)
      
    ## p) Calculate Model Sum of Squares from reduced model
      mss_snps_r_f <- sum((y_hat_snps_r_f - mean(y_snps_f))^2)
      
    ## q) calculate the extra-model mean square
      ms_model_snps_f = (mss_snps_f - mss_snps_r_f)/
        (df_mod_snps_f - df_mod_snps_r_f)
      
    ## r) calculate partial F-value that is attributable to instrument
      F_snps_f = ms_model_snps_f/ms_error_snps_f
      
    ## s) estimate the p-value for the F-value
      p_value_snps_f <- pf(F_snps_f, (df_mod_snps_f - df_mod_snps_r_f) , 
                           df_resid_snps_f,
                         lower.tail = FALSE)
      
    ## t) estimate the marginal and conditional R-square
      r.squaredGLMM(snps_plum_mod_f)
      
      
  ### 3.4 PGS effect on breast brightness - female
    ## a) transform snps to numeric - female df
      snps_plum_repro_f_df <- snps_plum_repro_f_df %>%
        mutate(across(c(kitlg.44277058:apc.46065597), 
                     as.numeric))
     
     # Lead and closest candidate proxy plumage snps
        # kitlg.44285485 = -4.5346084
        # plxnc1.46520400 = -2.44366682
        # slc45a3.19332900 = -3.227701
        # bnc2.33382167 = 0.06599935
        # gnaq.35901305 = -0.1531523
        # ror2.44206228 = -3.5401068
        # apc.46065597 = -4.7054829
     
    ## b) Beta values from Schield et al used to create PGS instrument - females
      betas_lead_snps <- c(-4.5346084, -2.44366682, -3.227701, 0.06599935, 
                  -0.1531523, -3.5401068, -4.7054829)
    
    ## c) create PGS instrument based on lead SNPs from Schield et al - females
      snps_plum_repro_f_df <- snps_plum_repro_f_df  %>%
        rowwise() %>%
        mutate(pgs.breast = sum(c_across(c(c(kitlg.44285485, plxnc1.46520400,
                                            slc45a3.19332900, bnc2.33382167,
                                            gnaq.35901305, ror2.44206228,
                                            apc.46065597))) 
                               * betas_lead_snps)) %>%
        ungroup()
     
    ## d) glmmTMB model for effect of lead SNPs or closest candidate proxy 
     # SNPs based on Schiled et al - female
      pgs_plum_mod_f <-  glmmTMB(breast.avg.bright ~ 
                                     pgs.breast,
                                     #sex +
                      # proportional covariance matrix based on relatedness
                      # propto(0 + sample.id | 1, relate_mat_plum_snps),
                                   data = snps_plum_repro_f_df,
                                   family = gaussian)

    ## e) model summary  
     summary(pgs_plum_mod_f)
     
    ## f) tidy PGS lead SNPs effect estimates - female
      pgs_plum_mod_tidy_f <- 
        tidy(pgs_plum_mod_f, 
              effects = 'fixed',
              conf.int = T,
              conf.method = 'Wald',
              conf.level = 0.95, 
              exponentiate = F) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value')) 
     
    
  ### 3.5 Calculate an F-statistic based on pgs_plum_mod_f and the
     #  pgs_plum_mod_reduced_f sums of squares
     
    ## a) run a reduced model with no pgs - female
      pgs_plum_mod_reduced_f <-  glmmTMB(breast.avg.bright ~
                                           #sex 
                                           + 1,
                      # proportional covariance matrix based on relatedness
                      # propto(0 + sample.id | 1, relate_mat_plum_pgs),
                                         data = snps_plum_repro_f_df,
                                         family = gaussian)
     
    ## b) use stats::anova to compare full and reduced model
      pgs_mod_comp_f <- tibble::tibble(anova(pgs_plum_mod_f,
                                             pgs_plum_mod_reduced_f))
     
    ## c) add model variable
      pgs_mod_comp_f$model <-  c('')
     
    ## d) move 'model' to first column
      pgs_mod_comp_f <- pgs_mod_comp_f %>%
        relocate(model, .before = `Df`)
     
    ## e) update values to indicate model
      pgs_mod_comp_f[1,1] <- 'reduced PGS model - female'
      pgs_mod_comp_f[2,1] <- 'full PGS model - female'
     
    ## f) Get fitted values (response scale)
      #pgs_plum_pred_val <- fitted.values(pgs_plum_mod_f)
      pgs_plum_pred_val_f <- predict(pgs_plum_mod_f,
                                     type = "response")
     
    ## g) Combine fitted values with original data
      snps_plum_repro_f_df$pgs.plum.pred.val <- pgs_plum_pred_val_f
     
    ## h) Determine residual and model degrees of freedom from full model
      df_resid_pgs_f <- get_df(pgs_plum_mod_f) # same as df.residual(model)
      df_mod_pgs_f <- get_df(pgs_plum_mod_f, 
                             type = "model") # same as attr(logLik(model), "df")
     
    ## i) Determine the model df from reduced model
      df_mod_pgs_r_f <- get_df(pgs_plum_mod_reduced_f, 
                               type = "model") # same as attr(logLik(model), "df")
     
    ## j) extract residual sums of squares from full model
      resid_pgs_f <- residuals(pgs_plum_mod_f, type = "response")
      resid_ss_pgs_f <- sum(resid_pgs_f^2)
     
    ## k) calculate the Mean Square error from full model
      ms_error_pgs_f = resid_ss_pgs_f/df_resid_pgs_f
     
    ## l) extract full model predicted values
      y_hat_pgs_f <- predict(pgs_plum_mod_f,
                             type = "response")
     
    ## m) extract reduced model predicted values
      y_hat_pgs_r_f <- predict(pgs_plum_mod_reduced_f,
                               type = "response")
     
    ## n) extract actual response values - same for models
      y_pgs_f <- model.response(model.frame(pgs_plum_mod_f))
     
    ## o) Calculate Model Sum of Squares from full model
      mss_pgs_f <- sum((y_hat_pgs_f - mean(y_pgs_f))^2)
     
    ## p) Calculate Model Sum of Squares from reduced model
      mss_pgs_r_f <- sum((y_hat_pgs_r_f - mean(y_pgs_f))^2)
     
    ## q) calculate the extra-model mean square
      ms_model_pgs_f = (mss_pgs_f - mss_pgs_r_f)/
       (df_mod_pgs_f - df_mod_pgs_r_f)
     
    ## r) calculate partial F-value that is attributable to instrument
      F_pgs_f = ms_model_pgs_f/ms_error_pgs_f
     
    ## s) estimate the p-value for the F-value
      p_value_pgs_f <- pf(F_pgs_f, (df_mod_pgs_f - df_mod_pgs_r_f) , 
                          df_resid_pgs_f,
                          lower.tail = FALSE)
     
    ## t) estimate the marginal and conditional R-square
      r.squaredGLMM(pgs_plum_mod_f)
  
      
      
  ### 4.7 Make formatted tables of the associations between the instrument
      # breast brightness
      ## a) row bind snps_mod_comparison and pgs_plum_mod_f,
      instrument_mod_sum_f <- bind_rows(snps_plum_mod_tidy_f, 
                                        pgs_plum_mod_tidy_f) 
      
      ## b) tidy instrument_mod_sum
      instrument_mod_sum <- instrument_mod_sum %>%
        rename(c('P' = 'Pr(>Chisq)')) %>%
        select(c(`model`, `Df`, `AIC`, `logLik`,`Chisq`, `Chi Df`, `P`))
      
      ## c) format table creating flextable object 
      instrument_mod_sum_frmt <- nice_table(instrument_mod_sum, 
                                            separate.header = T)
      
      ## c) preview the table as a word doc    
      #print(instrument_mod_sum_frmt, preview = 'docx')
      
      ## d) save instrument_mod_sum_frmt as a word doc
      flextable::save_as_docx(instrument_mod_sum_frmt, 
                              path = here('output/instrument_mod_sum_frmt.docx'))  
      
      
      ### 4.8 Tidy predicted breast brightness values from Step 1 for use in Step 2
      ## a) select the predicted values
      snps_pred_plum_repro_df <- snps_plum_repro_df %>%
        select(c(sample.id, snps.plum.pred.val,
                 pgs.plum.pred.val))
      
      ## b) Left join the predicted values for breast plumage brightness
      # to snps_plum_repro_f_df
      snps_plum_repro_f_df <- snps_plum_repro_f_df  %>%
        left_join(snps_pred_plum_repro_df,
                  by = c('sample.id' = 'sample.id'), 
                  copy = F) 
      
      ## c) Left join the predicted values for breast plumage brightness
      # to snps_plum_repro_m_df
      snps_plum_repro_m_df <- snps_plum_repro_m_df  %>%
        left_join(snps_pred_plum_repro_df,
                  by = c('sample.id' = 'sample.id'), 
                  copy = F) 
      
     
     
###############################################################################
##############     4. Step 1: Create male genetic instrument     ##############
###############################################################################           
      
  ### 4.1 Tidy snps_plum_repro_m_df
    ## a) transform snsp to factor (male df)
      snps_plum_repro_m_df <- snps_plum_repro_m_df %>%
        mutate(across(c(kitlg.44277058:apc.46065597),
                      ~factor(., levels =  c('0', '1', '2'))))
      
      
  ### 4.2 SNPs effect on breast brightness - male
    ## a) glmmTMB model for effect of lead SNPs or closest candidate proxy 
      # SNPs based on Schiled et al - male
      snps_plum_mod_m <-  glmmTMB(breast.avg.bright ~ 
                                    kitlg.44285485 +
                                    plxnc1.46520400 +
                                    slc45a3.19332900 +
                                    bnc2.33382167 +
                                    gnaq.35901305 +
                                    ror2.44206228 +
                                    apc.46065597,
                                  data = snps_plum_repro_m_df,
                                  family = gaussian)
      
    ## b) model summary - male 
      summary(snps_plum_mod_m)
      
    ## c) tidy effect estimates - male
      snps_plum_mod_tidy_m <- 
        tidy(snps_plum_mod_m, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = F) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value')) 
      # NOTE: NO significant snps:
      
      
  ### 4.3 Calculate an F-statistic based on snps_plum_mod_m and the
      #  snps_plum_mod_reduced_m sums of squares
      
    ## a) run a reduced model with no snps - male
      snps_plum_mod_reduced_m <-  glmmTMB(breast.avg.bright ~
                                            #sex 
                                            + 1,
                    # proportional covariance matrix based on relatedness
                    # propto(0 + sample.id | 1, relate_mat_plum_snps),
                                          data = snps_plum_repro_m_df,
                                          family = gaussian)
      
    ## b) use stats::anova to compare full and reduced model
      snps_mod_comp_m <- tibble::tibble(anova(snps_plum_mod_m,
                                              snps_plum_mod_reduced_m))
      
    ## c) add model variable
      snps_mod_comp_m$model <-  c('')
      
    ## d) move 'model' to first column
      snps_mod_comp_m <- snps_mod_comp_m %>%
        relocate(model, .before = `Df`)
      
    ## e) update values to indicate model
      snps_mod_comp_m[1,1] <- 'reduced SNPs model - male'
      snps_mod_comp_m[2,1] <- 'full SNPs model - male'
      
    ## f) Get fitted values (response scale)
      #snps_plum_pred_val <- fitted.values(snps_plum_mod_m)
      snps_plum_pred_val_m <- predict(snps_plum_mod_m,
                                      type = "response")
      
    ## g) Combine fitted values with original data
      snps_plum_repro_m_df$snps.plum.pred.val <- snps_plum_pred_val_m
      
    ## h) Determine residual and model degrees of freedom from full model
      df_resid_snps_m <- get_df(snps_plum_mod_m) # same as df.residual(model)
      df_mod_snps_m <- get_df(snps_plum_mod_m, 
                              type = "model") # same as attr(logLik(model), "df")
      
    ## i) Determine the model df from reduced model
      df_mod_snps_r_m <- get_df(snps_plum_mod_reduced_m, 
                                type = "model") # same as attr(logLik(model), "df")
      
    ## j) extract residual sums of squares from full model
      resid_snps_m <- residuals(snps_plum_mod_m, type = "response")
      resid_ss_snps_m <- sum(resid_snps_m^2)
      
    ## k) calculate the Mean Square error from full model
      ms_error_snps_m = resid_ss_snps_m/df_resid_snps_m
      
    ## l) extract full model predicted values
      y_hat_snps_m <- predict(snps_plum_mod_m,
                              type = "response")
      
    ## m) extract reduced model predicted values
      y_hat_snps_r_m <- predict(snps_plum_mod_reduced_m,
                                type = "response")
      
    ## n) extract actual response values - same for models
      y_snps_m <- model.response(model.frame(snps_plum_mod_m))
      
    ## o) Calculate Model Sum of Squares from full model
      mss_snps_m <- sum((y_hat_snps_m - mean(y_snps_m))^2)
      
    ## p) Calculate Model Sum of Squares from reduced model
      mss_snps_r_m <- sum((y_hat_snps_r_m - mean(y_snps_m))^2)
      
    ## q) calculate the extra-model mean square
      ms_model_snps_m = (mss_snps_m - mss_snps_r_m)/
        (df_mod_snps_m - df_mod_snps_r_m)
      
    ## r) calculate partial F-value that is attributable to instrument
      F_snps_m = ms_model_snps_m/ms_error_snps_m
      
    ## s) estimate the p-value for the F-value
      p_value_snps_m <- pf(F_snps_m, (df_mod_snps_m - df_mod_snps_r_m) , 
                           df_resid_snps_m,
                           lower.tail = FALSE)
      
    ## t) estimate the marginal and conditional R-square
      r.squaredGLMM(snps_plum_mod_m)
      
      
  ### 4.4 PGS effect on breast brightness - male
    ## a) transform snps to numeric - male df
      snps_plum_repro_m_df <- snps_plum_repro_m_df %>%
        mutate(across(c(kitlg.44277058:apc.46065597), 
                      as.numeric))
      
      # Lead and closest candidate proxy plumage snps
      # kitlg.44285485 = -4.5346084
      # plxnc1.46520400 = -2.44366682
      # slc45a3.19332900 = -3.227701
      # bnc2.33382167 = 0.06599935
      # gnaq.35901305 = -0.1531523
      # ror2.44206228 = -3.5401068
      # apc.46065597 = -4.7054829
      
    ## b) Beta values from Schield et al used to create PGS instrument - males
      betas_lead_snps <- c(-4.5346084, -2.44366682, -3.227701, 0.06599935, 
                           -0.1531523, -3.5401068, -4.7054829)
      
    ## c) create PGS instrument based on lead SNPs from Schield et al - males
      snps_plum_repro_m_df <- snps_plum_repro_m_df  %>%
        rowwise() %>%
        mutate(pgs.breast = sum(c_across(c(c(kitlg.44285485, plxnc1.46520400,
                                             slc45a3.19332900, bnc2.33382167,
                                             gnaq.35901305, ror2.44206228,
                                             apc.46065597))) 
                                * betas_lead_snps)) %>%
        ungroup()
      
      
    ## d) glmmTMB model for effect of lead SNPs or closest candidate proxy 
      # SNPs based on Schiled et al - male
      pgs_plum_mod_m <-  glmmTMB(breast.avg.bright ~ 
                                   pgs.breast,
                                 #sex +
                       # proportional covariance matrix based on relatedness
                       # propto(0 + sample.id | 1, relate_mat_plum_snps),
                                 data = snps_plum_repro_m_df,
                                 family = gaussian)
      
    ## e) model summary  
      summary(pgs_plum_mod_m)
      
    ## f) tidy PGS lead SNPs effect estimates - male
      pgs_plum_mod_tidy_m <- 
        tidy(pgs_plum_mod_m, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = F) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value')) 
      
      
  ### 4.5 Calculate an F-statistic based on pgs_plum_mod_m and the
      #  pgs_plum_mod_reduced_m sums of squares
      
    ## a) run a reduced model with no pgs - male
      pgs_plum_mod_reduced_m <-  glmmTMB(breast.avg.bright ~
                                           #sex 
                                           + 1,
                       # proportional covariance matrix based on relatedness
                       # propto(0 + sample.id | 1, relate_mat_plum_pgs),
                                         data = snps_plum_repro_m_df,
                                         family = gaussian)
      
    ## b) use stats::anova to compare full and reduced model
      pgs_mod_comp_m <- tibble::tibble(anova(pgs_plum_mod_m,
                                             pgs_plum_mod_reduced_m))
      
    ## c) add model variable
      pgs_mod_comp_m$model <-  c('')
      
    ## d) move 'model' to first column
      pgs_mod_comp_m <- pgs_mod_comp_m %>%
        relocate(model, .before = `Df`)
      
    ## e) update values to indicate model
      pgs_mod_comp_m[1,1] <- 'reduced PGS model - male'
      pgs_mod_comp_m[2,1] <- 'full PGS model - male'
      
    ## f) Get fitted values (response scale)
      #pgs_plum_pred_val <- fitted.values(pgs_plum_mod_m)
      pgs_plum_pred_val_m <- predict(pgs_plum_mod_m,
                                     type = "response")
      
    ## g) Combine fitted values with original data
      snps_plum_repro_m_df$pgs.plum.pred.val <- pgs_plum_pred_val_m
      
    ## h) Determine residual and model degrees of freedom from full model
      df_resid_pgs_m <- get_df(pgs_plum_mod_m) # same as df.residual(model)
      df_mod_pgs_m <- get_df(pgs_plum_mod_m, 
                         type = "model") # same as attr(logLik(model), "df")
      
    ## i) Determine the model df from reduced model
      df_mod_pgs_r_m <- get_df(pgs_plum_mod_reduced_m, 
                         type = "model") # same as attr(logLik(model), "df")
      
    ## j) extract residual sums of squares from full model
      resid_pgs_m <- residuals(pgs_plum_mod_m, type = "response")
      resid_ss_pgs_m <- sum(resid_pgs_m^2)
      
    ## k) calculate the Mean Square error from full model
      ms_error_pgs_m = resid_ss_pgs_m/df_resid_pgs_m
      
    ## l) extract full model predicted values
      y_hat_pgs_m <- predict(pgs_plum_mod_m,
                             type = "response")
      
    ## m) extract reduced model predicted values
      y_hat_pgs_r_m <- predict(pgs_plum_mod_reduced_m,
                               type = "response")
      
    ## n) extract actual response values - same for models
      y_pgs_m <- model.response(model.frame(pgs_plum_mod_m))
      
    ## o) Calculate Model Sum of Squares from full model
      mss_pgs_m <- sum((y_hat_pgs_m - mean(y_pgs_m))^2)
      
    ## p) Calculate Model Sum of Squares from reduced model
      mss_pgs_r_m <- sum((y_hat_pgs_r_m - mean(y_pgs_m))^2)
      
    ## q) calculate the extra-model mean square
      ms_model_pgs_m = (mss_pgs_m - mss_pgs_r_m)/
        (df_mod_pgs_m - df_mod_pgs_r_m)
      
    ## r) calculate partial F-value that is attributable to instrument
      F_pgs_m = ms_model_pgs_m/ms_error_pgs_m
      
    ## s) estimate the p-value for the F-value
      p_value_pgs_m <- pf(F_pgs_m, (df_mod_pgs_m - df_mod_pgs_r_m) , 
                          df_resid_pgs_m,
                          lower.tail = FALSE)
      
    ## t) estimate the marginal and conditional R-square
      r.squaredGLMM(pgs_plum_mod_m)      
      
  
     
###############################################################################
##############  #. Step 1: Create genetic instrument using elastic net  ##############
###############################################################################   
      
      
      # For Step 1, regress plumage traits on the candidate snps - the 
      # genetic instrument (e.g., plumage ~ SNPs for each plumage locus) 
      # to get predicted plumage values for use in Step 2.
      
    # Developing instrument methodology - Cnt'd
      # 6.  Tests of the best predictors of the plumage trait - for each  
          # locus where multiple candidate proxy SNPs identified in section 3,
          # use dredge to compare AICc for all pairwise combinations. Use
          # AICc summary to identify best fit models (delta AIC <2)
      # 7.  Tests of conditional associations - for each locus where multiple 
          # candidate proxy SNPs identified assess the 
          # conditional associations in each AICc determined best fit model
          # and make sure that multiple candidate proxy SNPs are independently
          # associated with the plumage trait (p < 0.2)
      # 8.  Additional decision criteria - if AIC and testing which SNPs are
          # are independently associated with plumage does not identify a
          # the best candidate proxy SNPs to carry forward to the instrument,
          # then compare F-statistics (full vs intercept only model)
      # 9.  F-test of multiple loci - Put all candidate proxy SNPs into a 
          # and determine what combination produces the best F-test.
      # 10. If there is a single polymorphic SNP at a locus, even if it is
          # not significant, add to the full model to obtain best F-test
    
    ## a) create a data frame that can be used in glmnet
      snps_mat_df <- snps_plum_repro_df %>%
        #filter(sex == 'm') %>%
        select(c(breast.avg.bright,kitlg.44277058:apc.46065597))
      
     
     ## b) use rsample to split the data frame into training and test data 
      snps_mat_df_split <- initial_split(snps_mat_df, prop = 0.75, 
                                         strata = breast.avg.bright)
      snps_mat_df_split_train <- training( snps_mat_df_split)
      snps_mat_df_split_test  <- testing( snps_mat_df_split)
      
    ## c) use recipes to specify the model parameterization formula and 
        # either standardize continuous predictors or create dummy variables 
        # for categorical  variables
     
      
      # list snps
      noquote(colnames(snps_mat_df))
      
      snps_mat_df_rec <- recipe(breast.avg.bright ~ 
                                  kitlg.44277058 + kitlg.44277464 + 
                                  kitlg.44277783 + kitlg.44278776 +  
                                  kitlg.44279385 + kitlg.44280329 + 
                                  kitlg.44284459 + kitlg.44284884 +
                                  kitlg.44285163 + kitlg.44285212 + 
                                  kitlg.44285279 + kitlg.44285485 +
                                  kitlg.44285821 + kitlg.44287574 +
                                  kitlg.44288056 + kitlg.44295622 +
                                  plcnc1.46510640 + plxnc1.46512949 +
                                  plxnc1.46513206 + plcnc1.46514224 +
                                  plxnc1.46518352 + plxnc1.46518900 +
                                  plxnc1.46520400 + plxnc1.46520543 + 
                                  plxnc1.46521340 + plxnc1.46521344 +
                                  plxnc1.46521577 + plxnc1.46521842 +
                                  plxnc1.46522706 + plxnc1.46524855 +
                                  plxnc1.46525798 + plxnc1.46526447 +
                                  plxnc1.46526602 + plcnc1.46527459 + 
                                  plxnc1.46528632 + plxnc1.46528737 +
                                  plxnc1.46528880 + plxnc1.46529619 +
                                  plxnc1.46529670 + slc45a3.19323094 +
                                  slc45a3.19324409 + slc45a3.19325078 +
                                  slc45a3.19328471 + slc45a3.19332900 +
                                  slc45a3.19333762 + bnc2.33382167 +
                                  bnc2.33385738 + bnc2.33387499 +
                                  gnaq.35898831 + gnaq.35901305 +
                                  gnaq.35909335 + ror2.44196487 +
                                  ror2.44201599 + ror2.44206228 +   
                                  ror2.44214826 + apc.46065597, 
                                data = snps_mat_df_split_train) %>%
                                step_dummy(all_nominal_predictors())
                              # step_normalize(all_predictors())
      
    ## d) Use tidymodels and parsnip to specify the model/estimator and 
        # to optimize tuning parameters (penalty and mixture)
      
      # For classification (use logistic_reg() for logistic regression)
      enet_spec <- linear_reg(
        penalty = tune(), # regularization (lambda)
        mixture = tune()  # balances alpha between 0 - ridge and 1 - lasso
      ) %>%
        set_engine('glmnet') %>%
        set_mode('regression')
        #set_mode("classification") # for logistic regression
      
    ## e) create a work flow that combines the model recipe/parameterization
        # and the specifications for the model/estimator
      enet_wf <- workflow() %>%
        add_recipe(snps_mat_df_rec) %>%
        add_model(enet_spec)
      
    ## f) specify cross validation folds (here 5) and set the tuning grid for 
        # different combinations of the penalty and mixture tunings
      snps_mat_df_folds <- vfold_cv(snps_mat_df_split_train, v = 5)
      
      enet_grid <- grid_regular( #grid_space_filling(
        penalty(), 
        mixture(), 
        levels = 5 # creates a tuning grid of 10 combinations of penalty and
                   # mixture ranges
        # penalty(range = c(-3, 2), trans = log10_trans()),
        # mixture(range = c(0, 1)),
        # size = 30
      )
      
    ## g) Search test combinations of the penalty and mixture values across
        # the cross-validation folds to tune the model by optimization of
        # of the root mean square error (RMSE)
      enet_tuned <- tune_grid(
        enet_wf,
        resamples = snps_mat_df_folds,
        grid = enet_grid,
        metrics = metric_set(rmse, mae, rsq)
        #metrics = metric_set(roc_auc, accuracy) # for logistic regression
      )
      
    ## h) Find the best parameters based on RMSE
      best_params <- select_best(enet_tuned, metric = 'rmse')
      best_params
      
      
    ## i) use the best tuning parameters and fit the model on the test data
      final_wf <- finalize_workflow(enet_wf, best_params)
      final_fit <- last_fit(final_wf, snps_mat_df_split)
      
    ## j) extract performance on the test set
      collect_metrics(final_fit)
      
    ## k) extract test set predictions
      test_predictions <- collect_predictions(final_fit)
      
    ## l) compare measures of model fit based combinations of 
        # penalty (regularization) and mixture (alpha)
      autoplot(enet_tuned)
      show_best(enet_tuned, metric = "rmse")
      
     ## m) extract the specific work flow with the optimized tuning parameters 
      fitted_wf <- extract_workflow(final_fit)
      
    ## n) extract the coefficients and save a data frame for the best predictors
        # of the outcome based on the best fit model
      fitted_wf_df <- fitted_wf %>% 
        extract_fit_parsnip() %>% 
        tidy() %>% 
        filter(term != "(Intercept)") %>% 
        mutate(absolute_estimate = abs(estimate)) %>% 
        arrange(desc(absolute_estimate))
      
      
      # install.packages('vip', repos = c("https://bgreenwell.r-universe.dev", "https://cloud.r-project.org"))
      # install.packages('vip')
      # library(vip)
      # 
      # fitted_wf %>% 
      #   extract_fit_parsnip() %>% 
      #   vip(num_features = 10, geom = "point") +
      #   labs(title = "Elastic Net Variable Importance",
      #        x = "Predictor Variables",
      #        y = "Importance (Absolute Coefficient)") +
      #   theme_minimal()
      
      
      
      
      # 
      # 
      # cv_fit <- cv.glmnet(
      #   x = X, 
      #   y = y, 
      #   alpha = 0.5,        # 0.5 balances Lasso (alpha=1) and Ridge (alpha=0)
      #   nfolds = 10,        # 10-fold cross-validation
      #   family = "gaussian" # use "binomial" for binary/case-control phenotypes
      # )
      # 
      # 
      # # Plot cross-validation error against log(lambda)
      # plot(cv_fit)
      # 
      # # Extract the best lambda (lambda.min or lambda.1se)
      # best_lambda <- cv_fit$lambda.min
      # 
      # # Extract non-zero SNP coefficients (the selected best predictors)
      # coefs <- coef(cv_fit, s = "lambda.min")
      # selected_snps <- coefs[coefs[, 1] != 0, ]
      # print(coefs)
      # # Print selected SNPs and their weights (excluding the intercept)
      # print(selected_snps[names(selected_snps) != "(Intercept)"])
      # 
      # 
      # ## a) transform snsp to factor (male df)
      # snps_plum_repro_m_df <- snps_plum_repro_df %>%
      #   filter(sex == 'm') %>%
      # 
      # class(snps_plum_repro_m_df$kitlg.44280329)
      # ## a) transform snsp to factor (male df)
      # snps_m_mat_df <- snps_plum_repro_m_df %>%
      #   select(c(breast.avg.bright,kitlg.44277058:apc.46065597))
      # 
      # X <- model.matrix(breast.avg.bright ~ ., data = snps_m_mat_df)[, -1]
      # y <- snps_m_mat_df$breast.avg.bright
      # alpha.target = 0.001
      # 
      # cv_fit <- cva.glmnet(breast.avg.bright ~ ., 
      #                      data = snps_m_mat_df, family = "gaussian", 
      #                      intercept = F)
      # 
      # get_best_alpha <- function(cva_model) {
      #   # Extract the minimum loss for each alpha tested
      #   losses <- sapply(cva_model$modlist, function(mod) min(mod$cvm))
      #   # Get the index of the absolute lowest loss
      #   best_idx <- which.min(losses)
      #   # Return the alpha value
      #   return(cva_model$alpha[best_idx])
      # }
      # 
      # best_alpha <- get_best_alpha(cv_fit)
      # print(paste("Optimal Alpha:", best_alpha))
      # 
      # 
      # # 4. Extract the model corresponding to the best alpha
      # best_model <- cv_fit$modlist[[which(cv_fit$alpha == best_alpha)]]
      # 
      # # 5. Extract selected variables (coefficients) at the optimal lambda
      # # Use "lambda.1se" for a more parsimonious/sparse model, or "lambda.min" 
      # best_coefficients <- coef(best_model, s = "lambda.1se")
      # 
      # # 6. Filter out variables shrunk to zero to see your selected features
      # coef_matrix <- as.matrix(best_coefficients)
      # selected_variables <- coef_matrix[coef_matrix[, 1] != 0, , drop = FALSE]
      # print(selected_variables)
      # 
      # library(glmnetUtils)
      # 
      # # Fits cv.glmnet across a sequence of alpha values automatically
      # cv_fit <- cva.glmnet(Y ~ ., data = X, intercept = F)
      # 
      # # Extract and view the optimal alpha value from the sequence
      # best_alpha <- cv_fit$alpha[which.min(sapply(cv_fit$modlist, function(x) min(x$cvm)))]
      # print(best_alpha)
      # 
      # # 5. Extract the optimal coefficients
      # best_coefs <- coef(cv_model, s = "lambda.min")
      # print(best_coefs)
      # summary(cv_model)
      # 
      # lasso_spec <- linear_reg(penalty = 0.1, mixture = 1) %>% # mixture = 1 means LASSO
      #   set_engine("glmnet")
      # 
      # 
      # 
      # 
      # 
      #       # Fit the model using a standard formula interface
      # lasso_fit <- lasso_spec %>%
      #   fit(breast.avg.bright ~ + plxnc1.46520543
      #       + plxnc1.46521340 # same
      #       + plxnc1.46521344 # same
      #       + plxnc1.46518352
      #       + plxnc1.46524855
      #       + plxnc1.46525798 # f-test
      #       + plxnc1.46526602
      #       + plxnc1.46512949
      #       + plxnc1.46528880, data = snps_plum_repro_m_df)
      # 
      # asso_coefs <- coef(lasso_fit, s = "lambda.1se")
      # print(lasso_fit)

    
        
 
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
  ### 4.7 Make formatted tables of the associations between the instrument
        # breast brightness
    ## a) row bind snps_mod_comparison and pgs_mod_comparison,
      instrument_mod_sum <- bind_rows(snps_mod_comparison, 
                                      pgs_mod_comparison) 
      
    ## b) tidy instrument_mod_sum
      instrument_mod_sum <- instrument_mod_sum %>%
        rename(c('P' = 'Pr(>Chisq)')) %>%
        select(c(`model`, `Df`, `AIC`, `logLik`,`Chisq`, `Chi Df`, `P`))

    ## c) format table creating flextable object 
      instrument_mod_sum_frmt <- nice_table(instrument_mod_sum, 
                                         separate.header = T)
      
    ## c) preview the table as a word doc    
      #print(instrument_mod_sum_frmt, preview = 'docx')
      
    ## d) save instrument_mod_sum_frmt as a word doc
      flextable::save_as_docx(instrument_mod_sum_frmt, 
                        path = here('output/instrument_mod_sum_frmt.docx'))  
      
      
  ### 4.8 Tidy predicted breast brightness values from Step 1 for use in Step 2
    ## a) select the predicted values
      snps_pred_plum_repro_df <- snps_plum_repro_df %>%
        select(c(sample.id, snps.plum.pred.val,
                 pgs.plum.pred.val))
    
    ## b) Left join the predicted values for breast plumage brightness
        # to snps_plum_repro_f_df
      snps_plum_repro_f_df <- snps_plum_repro_f_df  %>%
        left_join(snps_pred_plum_repro_df,
                  by = c('sample.id' = 'sample.id'), 
                  copy = F) 
    
    ## c) Left join the predicted values for breast plumage brightness
      # to snps_plum_repro_m_df
      snps_plum_repro_m_df <- snps_plum_repro_m_df  %>%
        left_join(snps_pred_plum_repro_df,
                  by = c('sample.id' = 'sample.id'), 
                  copy = F) 

      

###############################################################################
##############                  5. Export data                   ##############
###############################################################################
      
  ### 5.1 Export data to an RData file 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      
    ## a) Save and export relatedness data
      save(file = here('data/5_pred_vals_repro_data.RData'), 
           list = c('snps_plum_repro_df', 'snps_plum_repro_f_df', 
                    'snps_plum_repro_m_df', 'relate_mat_plum_snps',
                    'relate_mat_plum_repro_f', 'relate_mat_plum_repro_m'))
      
      

    

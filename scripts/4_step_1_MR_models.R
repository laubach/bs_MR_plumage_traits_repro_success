################################################################################
#############                                                      #############
#############     Assessing the causal effect of plumage traits    #############
#############     on reproduction using Mendelian randomization    #############   
#############                                                      #############
#############       4. Step 1: Mendelian randomization models      #############
#############                                                      #############
#############                   By: Zach Laubach                   #############
#############                 created: 9 Sept 2025                 #############
#############              last updated: 10 May 2026               #############
################################################################################


  ### PURPOSE: Tidy the 
  
  # Code Blocks
    # 1. Configure work space
    # 2. Import data 
    # 3. Assess strength of snp effects on traits
    # 4. Step 1: Create genetic instrument
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
    
      # load here package
        library ('here')
      
    ## b) Modelling Packages    
       # load glmmTMB package
        library ('glmmTMB')
      
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
      load(here('data/4_MR_data.RData'))
      
    ## b) Load relatedness data
      load(here('data/3-4_relate_df.RData'))
      
    

###############################################################################
##############    3. Assess strength of snp effects on traits    ##############
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
      
      
    
  ### 3.1 Tidy snps_plum_repro_df
    ## a) transform snsp to factor 
      snps_plum_repro_df <- snps_plum_repro_df %>%
        mutate(across(c(kitlg.44277058:apc.46065597),
                      ~factor(., levels =  c('0', '1', '2'))))
      
     
  ### 3.2 Assess KITLG snp associations with breast brightness
    ## a) make a list of snp names from KITLG
      kitlg_snps <- snps_plum_repro_df %>%
        select(kitlg.44277058:kitlg.44285485) %>%
        colnames()
      
      kitlg_snps <- as.list(kitlg_snps)
      
    ## b) Use map to run multiple glmmTMB models - one for each snp
      kitlg_snps_plum_mods <- map(kitlg_snps, ~ {
        # specify model formula including a proportional covariance matrix 
        # based on relatedness
        # test for sex x snp interactions - NO sig intx p < 0.1
        #form <- as.formula(paste('breast.avg.bright ~ sex *', .x,
        # if no intx, then control for sex  
        form <- as.formula(paste('breast.avg.bright ~ sex +', .x,
                        '+ propto(0 + sample.id | 1, relate_mat_plum_snps)'))
        glmmTMB(
          formula = form,
          data = snps_plum_repro_df,
          family = gaussian(),
        )
      })
      
    ## c) View Summary of first model
      summary(kitlg_snps_plum_mods[[1]])
      
    ## d) Test to make sure map working by comparing to a glmmTMB model for 
      # effect of kitlg.44277058 on breast brightness - shoud be the same
      # as the results from c
      snps_plum_mod <-  glmmTMB(breast.avg.bright ~ 
                                 sex + kitlg.44277058 +
                        # proportional covariance matrix based on relatedness
                          propto(0 + sample.id | 1, relate_mat_plum_snps),
                                data = snps_plum_repro_df,
                                family = gaussian)
      summary(snps_plum_mod)
      
    ## e) Use broom.mixed to tidy effect estimates
      kitlg_snps_plum_mods_tidy <- kitlg_snps_plum_mods %>%
        map_dfr(tidy, effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = F) %>%
        filter(!(term == '(Intercept)' | term == 'sexm')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value'))
      # NOTE: significant snps: 
        #kitlg.44277058 
      
      
      
  ### 3.3 Assess PLXNC1 snp associations with breast brightness
    ## a) make a list of snp names from PLXNC1
      plxnc1_snps <- snps_plum_repro_df %>%
        select(plxnc1.46512949:plxnc1.46529670) %>%
        colnames()
      
      plxnc1_snps <- as.list(plxnc1_snps)
      
    ## b) Use map to run multiple glmmTMB models - one for each snp
      plxnc1_snps_plum_mods <- map(plxnc1_snps, ~ {
        # specify model formula including a proportional covariance matrix 
        # based on relatedness
        # test for sex x snp interactions - NO sig intx p < 0.1
        #form <- as.formula(paste('breast.avg.bright ~ sex *', .x,
        # if no intx, then control for sex  
        form <- as.formula(paste('breast.avg.bright ~ sex +', .x,
                        '+ propto(0 + sample.id | 1, relate_mat_plum_snps)'))
        glmmTMB(
          formula = form,
          data = snps_plum_repro_df,
          family = gaussian(),
        )
      })
      
    ## c) Use broom.mixed to tidy effect estimates
      plxnc1_snps_plum_mods_tidy <- plxnc1_snps_plum_mods %>%
        map_dfr(tidy, effects = 'fixed',
                conf.int = T,
                conf.method = 'Wald',
                conf.level = 0.95, 
                exponentiate = F) %>%
        filter(!(term == '(Intercept)'| term == 'sexm')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value'))
      # NOTE: significant snps: 
        # sexm:plxnc1.46525798
        # sexm:plxnc1.46521842
        # sexm:plxnc1.46518352
        # sexm:plxnc1.46526602
        # plxnc1.46520543
        # plxnc1.46524855
        # plxnc1.46520400
       

      
  ### 3.4 SLC45A3 effect on breast brightness - both sexes model
    ## a) glmmTMB model for effect of SLC45A3 on breast brightness
      slc45a3_plum_mod <-  glmmTMB(breast.avg.bright ~ 
                              slc45a3.19332900
                              # test for intx - no intx
                              # * sex +
                              + sex  +
                          # proportional covariance matrix based on relatedness
                              propto(0 + sample.id | 1, relate_mat_plum_snps),
                              data = snps_plum_repro_df,
                                 family = gaussian)
      
    ## b) model summary  
      summary(slc45a3_plum_mod)
      
    ## c) tidy effect estimates
      slc45a3_plum_mod_tidy <- 
        tidy(slc45a3_plum_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = F) %>%
        filter(!(term == '(Intercept)'| term == 'sexm')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value'))
      # NOTE: NO significant snps: 
      
    
  ### 3.5 BNC2 effect on breast brightness - both sexes model
    ## a) glmmTMB model for effect of BNC2 on breast brightness
      bnc2_plum_mod <-  glmmTMB(breast.avg.bright ~ 
                               bnc2.33387499
                               # test for intx - no intx
                               #* sex +
                               + sex  +
                         # proportional covariance matrix based on relatedness
                               propto(0 + sample.id | 1, relate_mat_plum_snps),
                               data = snps_plum_repro_df,
                                   family = gaussian)
      
    ## b) model summary  
      summary(bnc2_plum_mod)
      
    ## c) tidy effect estimates
      bnc2_plum_mod_tidy <- 
        tidy(bnc2_plum_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = F) %>%
        filter(!(term == '(Intercept)'| term == 'sexm')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value')) 
      # NOTE: NO significant snps:
      
      
  ### 3.6 Assess ROR2 snp associations with breast brightness
    ## a) make a list of snp names from ROR2
      ror2_snps <- snps_plum_repro_df %>%
        select(ror2.44196487:ror2.44214826) %>%
        colnames()
     
      ror2_snps <- as.list(ror2_snps)
      
    ## b) Use map to run multiple glmmTMB models - one for each snp
      ror2_snps_plum_mods <- map(ror2_snps, ~ {
        # specify model formula including a proportional covariance matrix 
        # based on relatedness
        # test for sex x snp interactions - NO sig intx p < 0.1
        #form <- as.formula(paste('breast.avg.bright ~ sex *', .x,
        # if no intx, then control for sex  
        form <- as.formula(paste('breast.avg.bright ~ sex +', .x,
                          '+ propto(0 + sample.id | 1, relate_mat_plum_snps)'))
        glmmTMB(
          formula = form,
          data = snps_plum_repro_df,
          family = gaussian(),
        )
      })
      
    ## c) Use broom.mixed to tidy effect estimates
      ror2_snps_plum_mods_tidy <- ror2_snps_plum_mods %>%
        map_dfr(tidy, effects = 'fixed',
                conf.int = T,
                conf.method = 'Wald',
                conf.level = 0.95, 
                exponentiate = F) %>%
        filter(!(term == '(Intercept)'| term == 'sexm')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value'))
      # NOTE: significant snps: 
        # ror2.44206228
        # ror2.44196487
      
      
  ### 3.7 APC effect on breast brightness - both sexes model
    ## a) glmmTMB model for effect of APC on breast brightness
      apc_plum_mod <-  glmmTMB(breast.avg.bright ~ 
                                 apc.46065597
                               # test for intx - no intx
                               #* sex +
                               + sex  +
                          # proportional covariance matrix based on relatedness
                                propto(0 + sample.id | 1, relate_mat_plum_snps),
                                data = snps_plum_repro_df,
                                  family = gaussian)
     
    ## b) model summary  
      summary(apc_plum_mod)
      
    ## c) tidy effect estimates
      apc_plum_mod_tidy <- 
        tidy(apc_plum_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = F) %>%
        filter(!(term == '(Intercept)'| term == 'sexm')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value'))  
      # NOTE: NO significant snps:
      
      
  ### 3.8 Assess if multiple snps from the same candidate region have 
        # independent effect of plumage (where snp p < 0.1)
    ## a) glmmTMB model for effect of PLXNC1 snps on breast brightness
      # sexm:plxnc1.46525798
      # sexm:plxnc1.46521842
      # sexm:plxnc1.46518352
      # sexm:plxnc1.46526602
      # plxnc1.46520543
      # plxnc1.46524855
      # plxnc1.46520400
      plxnc1_all_plum_mods <-  glmmTMB(breast.avg.bright ~ 
                                plxnc1.46520400 
                              + plxnc1.46520543 
                              + plxnc1.46524855 
                              + plxnc1.46525798 * sex
                              + plxnc1.46521842 * sex
                              + plxnc1.46518352 * sex
                              + plxnc1.46526602 * sex
                              + sex  +
                      # proportional covariance matrix based on relatedness
                              propto(0 + sample.id | 1, relate_mat_plum_snps),
                                data = snps_plum_repro_df,
                                  family = gaussian, na.action = "na.fail")
      
    ## b) use dredge from MUMIn package to run all model combinations
      plxnc1_aic_models <- dredge(plxnc1_all_plum_mods)
      
    ## c) view the top models
      head(plxnc1_aic_models)
      # NOTE: model with single snp plxnc1.46524855 is best fit based on AICc
            # no other model has AICc < 1
      
    ## d) check if other snps independently predict breast brightness
      summary(plxnc1_all_plum_mods)
      # NOTE: snps that independently predict breast brightness:
        # plxnc1.465205432
        # sexm:plxnc1.46521842
  
    ## e) glmmTMB model for effect of ROR2 snps on breast brightness
        # ror2.44206228
        # ror2.44196487
      ror2_all_plum_mods <-  glmmTMB(breast.avg.bright ~ 
                                ror2.44206228
                              + ror2.44196487
                              + sex  +
                          # proportional covariance matrix based on relatedness
                              propto(0 + sample.id | 1, relate_mat_plum_snps),
                                      data = snps_plum_repro_df,
                                      family = gaussian, na.action = "na.fail")
      
    ## b) use dredge from MUMIn package to run all model combinations
      ror2_aic_models <- dredge(ror2_all_plum_mods)
      
    ## c) view the top models
      head(ror2_aic_models)
      # note model with both snps ror2.44206228 and ror2.44196487 is 
      # among the best fit based on AICc - difference in AICc < 1
    
    ## d) check if snps indendently predict breast brightness
      summary(ror2_all_plum_mods)
      # NOTE: both snps independently predict breast brightens, so pick 
            # top model with both snps and based on AICc
      
    
      
###############################################################################
##############       4. Step 1: Create genetic instrument        ##############
###############################################################################   
      
      # For Step 1, regress plumage traits on the candidate snp - the 
      # genetic instrument (e.g., plumage ~ SNPs for each plumage locus) 
      # to get predicted plumage values for use in Step 2.        
      # Plumage snps
        # kitlg.44277058
        # plxnc1.46524855
        # plxnc1.46520543
        # sexm:plxnc1.46521842
        # ror2.44206228
        # ror2.44196487
      
      
  ### 4.1 KITLG, PLXNC1 and ROR2 effects on breast brightness 
    ## a) glmmTMB model for effect of KITLG, PLXNC1, and ROR2 on breast bright
      snps_plum_pred_mod <-  glmmTMB(breast.avg.bright ~ 
                                kitlg.44277058
                              + plxnc1.46524855
                              + plxnc1.46520543
                              + plxnc1.46521842 * sex 
                              + ror2.44206228
                              + ror2.44196487
                               + sex +
                        # proportional covariance matrix based on relatedness
                               propto(0 + sample.id | 1, relate_mat_plum_snps),
                               data = snps_plum_repro_df,
                               family = gaussian)
      
    ## b) model summary  
      summary(snps_plum_pred_mod)
      
    ## c) tidy effect estimates
      snps_plum_pred_mod_tidy <- 
        tidy(snps_plum_pred_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = F) %>%
        filter(!(term == '(Intercept)'| term == 'sexm')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value'))  
      
    ## d) label model
      # snps_plum_mod_tidy$model <- 
      #   c('KITLG')  
      
    ## e) Get fitted values (response scale)
      #snps_plum_pred_val <- fitted.values(snps_plum_pred_mod)
      snps_plum_pred_val <- predict(snps_plum_pred_mod, 
                                           type = "response")
      
    ## f) Combine fitted values with original data
      snps_plum_repro_df$snps.plum.pred.val <- snps_plum_pred_val
      
      
  ### 4.2 Create genetic instrument index variables
    ## a) transform snps to numeric 
      snps_plum_repro_df <- snps_plum_repro_df %>%
        mutate(across(c(kitlg.44277058, plxnc1.46524855, plxnc1.46520543,
                        plxnc1.46521842, ror2.44206228, ror2.44196487), 
                      as.numeric))
  
    ## b) Use allele dosages to calculate an breast plumage brightness 
      # instrument based on the imputed snps that were externally identified 
      # in hybrid GWAS from Schield et al 2024
      snps_plum_repro_df <- snps_plum_repro_df  %>%
        rowwise() %>%
        mutate(snps.sum.breast = sum(c_across(c(kitlg.44277058, 
                                      plxnc1.46524855, plxnc1.46520543,
                                      plxnc1.46521842, ror2.44206228, 
                                      ror2.44196487)))) %>%
        ungroup()
      
    ## c) Use imputed allele dosages and effect sizes from (Schield et al 2024)
      # to calculate a polygenic score (PGS) for breast plumage snps
      # PGS = sum(beta1*dosage1 + beta2*dosage2 ...) where beta is a 
      # population level estimate, and dosage is a vector for each individual
      
      # kitlg.44277058 beta = -2.99372967
      # plxnc1.46524855 beta = -2.13646641
      # plxnc1.46520543 beta = -2.13646641
      # plxnc1.46521842 beta = -0.88620883
      # ror2.44206228 beta = -3.5401068
      # ror2.44196487 beta = -3.3972115        

      # constants
      betas <- c(-2.99372967, -2.13646641, -2.13646641, -0.88620883, 
                 -3.5401068, -3.3972115)
      
      snps_plum_repro_df <- snps_plum_repro_df  %>%
        rowwise() %>%
        mutate(pgs.breast = sum(c_across(c(c(kitlg.44277058, 
                                             plxnc1.46524855, plxnc1.46520543,
                                             plxnc1.46521842, ror2.44206228, 
                                             ror2.44196487))) * betas)) %>%
        ungroup()
      
     
  ### 4.3 KITLG, PLXNC1, and ROR2 snps sum dosage index effects on 
        # breast brightness 
    ## a) glmmTMB model for effect of KITLG, PLXNC1 and ROR2 sum index on 
        # breast brightness
      snps_sum_plum_pred_mod <-  glmmTMB(breast.avg.bright ~ 
                              snps.sum.breast
                              + sex +
                        # proportional covariance matrix based on relatedness
                              propto(0 + sample.id | 1, relate_mat_plum_snps),
                              data = snps_plum_repro_df,
                                  family = gaussian)
      
    ## b) model summary  
      summary(snps_sum_plum_pred_mod)
      
    ## c) tidy effect estimates
      snps_sum_plum_pred_mod_tidy <- 
        tidy(snps_sum_plum_pred_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = F) %>%
        filter(!(term == '(Intercept)'| term == 'sexm')) %>%
        select(c('term', 'estimate', #'std.error', 
                 'conf.low', 'conf.high', 
                 #'statistic', 
                 'p.value'))  
      
    ## d) label model
      snps_sum_plum_pred_mod_tidy$model <-
        c('SNP dosage index')
      
    ## e) tidy snps_sum_plum_pred_mod_tidy by making a new column '95% CI'
      snps_sum_plum_pred_mod_tidy$'95% CI' = 
        paste(round(snps_sum_plum_pred_mod_tidy$conf.low, 3),
              round(snps_sum_plum_pred_mod_tidy$conf.high, 3), 
              sep=", ")
      
    ## f) tidy snps_sum_plum_pred_mod_tidy by removing extra columns and rows
      snps_sum_plum_pred_mod_tidy <- snps_sum_plum_pred_mod_tidy %>%
        rename(c(`β` = `estimate`, `P-value` = `p.value`)) %>%
        select(c(`β`, `95% CI`, `P-value`, `model`))
               
    ## g) move 'model' to first column
      snps_sum_plum_pred_mod_tidy <- snps_sum_plum_pred_mod_tidy %>%
                 relocate(model, .before = `β`) 
      
    ## h) Get fitted values (response scale)
      snps_sum_plum_pred_val <- predict(snps_sum_plum_pred_mod, 
                                        type = "response")

    ## i) Combine fitted values with original data
      snps_plum_repro_df$snps.sum.plum.pred.val <- snps_sum_plum_pred_val
      
      
  ### 4.4 KITLG, PLXNC1, and ROR2 snps PGS index effects on breast brightness 
    ## a) glmmTMB model for effect of KITLG, PLXNC1 and ROR2 PGS on breast brightness
      pgs_plum_pred_mod <-  glmmTMB(breast.avg.bright ~ 
                              pgs.breast
                              + sex +
                         # proportional covariance matrix based on relatedness
                              propto(0 + sample.id | 1, relate_mat_plum_snps),
                              data = snps_plum_repro_df,
                                  family = gaussian)
      
    ## b) model summary  
      summary(pgs_plum_pred_mod)
      
    ## c) tidy effect estimates
      pgs_plum_pred_mod_tidy <- 
        tidy(pgs_plum_pred_mod, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = F) %>%
        filter(!(term == '(Intercept)'| term == 'sexm')) %>%
        select(c('term', 'estimate', #'std.error', 
                 'conf.low', 'conf.high', 
                 #'statistic', 
                 'p.value'))  
      
    ## d) label model
      pgs_plum_pred_mod_tidy$model <-
        c('polygenic score (PGS)')
      
    ## e) tidy pgs_plum_pred_mod_tidy by making a new column '95% CI'
      pgs_plum_pred_mod_tidy$'95% CI' = 
        paste(round(pgs_plum_pred_mod_tidy$conf.low, 3),
              round(pgs_plum_pred_mod_tidy$conf.high, 3), 
              sep=", ")
      
    ## f) tidy pgs_plum_pred_mod_tidy by removing extra columns and rows
      pgs_plum_pred_mod_tidy <- pgs_plum_pred_mod_tidy %>%
        rename(c(`β` = `estimate`, `P-value` = `p.value`)) %>%
        select(c(`β`, `95% CI`, `P-value`, `model`))
      
    ## g) move 'model' to first column
      pgs_plum_pred_mod_tidy <- pgs_plum_pred_mod_tidy %>%
        relocate(model, .before = `β`) 
      
      
    ## h) Get fitted values (response scale)
      pgs_plum_pred_val <- predict(pgs_plum_pred_mod, 
                                        type = "response")
      
    ## i) Combine fitted values with original data
      snps_plum_repro_df$pgs.plum.pred.val <- pgs_plum_pred_val
      
  
  ### 4.5 Make formatted tables of the associations between the instrument
        # breast brightness
    ## a) row bind snps_sum_plum_pred_mod_tidy and pgs_plum_pred_mod_tidy,
      instrument_mod_sum <- bind_rows(snps_sum_plum_pred_mod_tidy, 
                                   pgs_plum_pred_mod_tidy) 
      
    ## b) format table creating flextable object 
      instrument_mod_sum_frmt <- nice_table(instrument_mod_sum, 
                                         separate.header = T)
      
    ## c) preview the table as a word doc    
      #print(instrument_mod_sum_frmt, preview = 'docx')
      
    ## d) save instrument_mod_sum_frmt as a word doc
      flextable::save_as_docx(instrument_mod_sum_frmt, 
                        path = here('output/instrument_mod_sum_frmt.docx'))  
      
      
  ### 4.6 Tidy predicted breast brightness values from Step 1 for use in Step 2
    ## a) select the predicted values
      snps_pred_plum_repro_df <- snps_plum_repro_df %>%
        select(c(sample.id, snps.plum.pred.val, snps.sum.plum.pred.val,
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
      
      

    

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
    
  ### 3.1 Tidy snps_plum_repro_df
    ## a) transform snsp to factor 
      snps_plum_repro_df <- snps_plum_repro_df %>%
        mutate(across(c(kitlg.44277058:apc.46065597),
                      ~factor(., levels =  c('0', '1', '2'))))
      
    ## b) transform snsp to factor (female df)
      snps_plum_repro_f_df <- snps_plum_repro_f_df %>%
        mutate(across(c(kitlg.44277058:apc.46065597),
                      ~factor(., levels =  c('0', '1', '2'))))
      
    ## a) transform snsp to factor (male df)
      snps_plum_repro_m_df <- snps_plum_repro_m_df %>%
        mutate(across(c(kitlg.44277058:apc.46065597),
                      ~factor(., levels =  c('0', '1', '2'))))
      
     
  ### 3.2 Assess KITLG snp associations with breast brightness
    ## a) make a list of snp names from KITLG
      kitlg_snps <- snps_plum_repro_df %>%
        select(kitlg.44277058:kitlg.44285485) %>%
        colnames()
      
      kitlg_snps <- as.list(kitlg_snps)
      
    # ## b) Use map to run multiple glmmTMB models - one for each snp
    #   kitlg_snps_plum_mods <- map(kitlg_snps, ~ {
    #     # specify model formula including a proportional covariance matrix 
    #     # based on relatedness
    #     # test for sex x snp interactions - NO sig intx p < 0.1
    #     #form <- as.formula(paste('breast.avg.bright ~ sex *', .x,
    #     # if no intx, then control for sex  
    #     form <- as.formula(paste('breast.avg.bright ~ sex +', .x,
    #                     '+ propto(0 + sample.id | 1, relate_mat_plum_snps)'))
    #     glmmTMB(
    #       formula = form,
    #       data = snps_plum_repro_df,
    #       family = gaussian(),
    #     )
    #   })
    #   
    # ## c) View Summary of first model
    #   summary(kitlg_snps_plum_mods[[1]])
    #   
    # ## d) Test to make sure map working by comparing to a glmmTMB model for 
    #   # effect of kitlg.44277058 on breast brightness - shoud be the same
    #   # as the results from c
    #   snps_plum_mod <-  glmmTMB(breast.avg.bright ~ 
    #                              sex + kitlg.44277058 +
    #                     # proportional covariance matrix based on relatedness
    #                       propto(0 + sample.id | 1, relate_mat_plum_snps),
    #                             data = snps_plum_repro_df,
    #                             family = gaussian)
    #   summary(snps_plum_mod)
      
      
    ## b) Use map to run multiple glmmTMB models - one for each snp - females
      kitlg_snps_plum_mods_f <- map(kitlg_snps, ~ {
        form <- as.formula(paste('breast.avg.bright ~', .x
       ))
        glmmTMB(
          formula = form,
          data = snps_plum_repro_f_df,
          family = gaussian(),
        )
      })
      
    ## c) Use broom.mixed to tidy effect estimates for females
      kitlg_snps_plum_mods_tidy_f <- kitlg_snps_plum_mods_f %>%
        map_dfr(tidy, effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = F) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value'))
      # NOTE: NO significant snps
      
    ## d) Use map to run multiple glmmTMB models - one for each snp - males
      kitlg_snps_plum_mods_m <- map(kitlg_snps, ~ {
        form <- as.formula(paste('breast.avg.bright ~', .x
        ))
        glmmTMB(
          formula = form,
          data = snps_plum_repro_m_df,
          family = gaussian(),
        )
      })
      
    ## e) Use broom.mixed to tidy effect estimates for females
      kitlg_snps_plum_mods_tidy_m <- kitlg_snps_plum_mods_m %>%
        map_dfr(tidy, effects = 'fixed',
                conf.int = T,
                conf.method = 'Wald',
                conf.level = 0.95, 
                exponentiate = F) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value'))
      # NOTE: NO significant snps

      
  ### 3.3 Assess PLXNC1 snp associations with breast brightness
    ## a) make a list of snp names from PLXNC1
      plxnc1_snps <- snps_plum_repro_df %>%
        select(plxnc1.46512949:plxnc1.46529670) %>%
        colnames()
      
      plxnc1_snps <- as.list(plxnc1_snps)
      
    ## b) Use map to run multiple glmmTMB models - one for each snp - females
      plxnc1_snps_plum_mods_f <- map(plxnc1_snps, ~ {
        form <- as.formula(paste('breast.avg.bright ~', .x))
        glmmTMB(
          formula = form,
          data = snps_plum_repro_f_df,
          family = gaussian(),
        )
      })
      
     ## c) Use broom.mixed to tidy effect estimates - females
      plxnc1_snps_plum_mods_tidy_f <- plxnc1_snps_plum_mods_f %>%
        map_dfr(tidy, effects = 'fixed',
                conf.int = T,
                conf.method = 'Wald',
                conf.level = 0.95, 
                exponentiate = F) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value'))
      # NOTE: significant snps: 
      # plxnc1.
      
      
    ## d) Use map to run multiple glmmTMB models - one for each snp - males
      plxnc1_snps_plum_mods_m <- map(plxnc1_snps, ~ {
        form <- as.formula(paste('breast.avg.bright ~', .x))
        glmmTMB(
          formula = form,
          data = snps_plum_repro_m_df,
          family = gaussian(),
        )
      })
      
    ## c) Use broom.mixed to tidy effect estimates - males
      plxnc1_snps_plum_mods_tidy_m <- plxnc1_snps_plum_mods_m %>%
        map_dfr(tidy, effects = 'fixed',
                conf.int = T,
                conf.method = 'Wald',
                conf.level = 0.95, 
                exponentiate = F) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value'))
      # NOTE: significant snps: 
      # plxnc1.46520400
       
      # ## a) make a list of snp names from PLXNC1
      #   plxnc1_snps <- snps_plum_repro_df %>%
      #     select(plxnc1.46512949:plxnc1.46529670) %>%
      #     colnames()
      #   
      #   plxnc1_snps <- as.list(plxnc1_snps)
      #   
      # ## b) Use map to run multiple glmmTMB models - one for each snp
      #   plxnc1_snps_plum_mods <- map(plxnc1_snps, ~ {
      #     # specify model formula including a proportional covariance matrix 
      #     # based on relatedness
      #     # test for sex x snp interactions - record sig intx p < 0.1
      #     # form <- as.formula(paste('breast.avg.bright ~ sex *', .x,
      #     # if no intx, then control for sex  
      #      form <- as.formula(paste('breast.avg.bright ~ sex +', .x,
      #                     '+ propto(0 + sample.id | 1, relate_mat_plum_snps)'))
      #     glmmTMB(
      #       formula = form,
      #       data = snps_plum_repro_df,
      #       family = gaussian(),
      #     )
      #   })
      #   
      # ## c) Use broom.mixed to tidy effect estimates
      #   plxnc1_snps_plum_mods_tidy <- plxnc1_snps_plum_mods %>%
      #     map_dfr(tidy, effects = 'fixed',
      #             conf.int = T,
      #             conf.method = 'Wald',
      #             conf.level = 0.95, 
      #             exponentiate = F) %>%
      #     filter(!(term == '(Intercept)'| term == 'sexm')) %>%
      #     select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
      #              'statistic', 'p.value'))
      #   # NOTE: significant snps: 
      #     # sexm:plxnc1.46525798
      #     # sexm:plxnc1.46518352
      #     # plxnc1.46521340
      #     # plxnc1.46521344
      #     # plxnc1.46524855
      #     # plxnc1.46513206
      #     # plxnc1.46520400
    
        
  ### 3.4 SLC45A3 effect on breast brightness - both sexes model
      
    ## a) glmmTMB model for effect of SLC45A3 on breast brightness - female
      slc45a3_plum_mod_f <-  glmmTMB(breast.avg.bright ~ 
                                     slc45a3.19332900,
                         data = snps_plum_repro_f_df,
                         family = gaussian)
      
    ## b) model summary  
      summary(slc45a3_plum_mod_f)
      
    ## c) tidy effect estimates - female
      slc45a3_plum_mod_tidy_f <- 
        tidy(slc45a3_plum_mod_f, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = F) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value'))
      
    ## d) glmmTMB model for effect of SLC45A3 on breast brightness - male
      slc45a3_plum_mod_m <-  glmmTMB(breast.avg.bright ~ 
                                       slc45a3.19332900,
                                     data = snps_plum_repro_m_df,
                                     family = gaussian)
      
    ## e) model summary  
      summary(slc45a3_plum_mod_m)
      
    ## f) tidy effect estimates - male
      slc45a3_plum_mod_tidy_m <- 
        tidy(slc45a3_plum_mod_m, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = F) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value'))
      
    # ## a) glmmTMB model for effect of SLC45A3 on breast brightness
    #   slc45a3_plum_mod <-  glmmTMB(breast.avg.bright ~ 
    #                           slc45a3.19332900
    #                           # test for intx - no intx
    #                           # * sex +
    #                           + sex  +
    #                       # proportional covariance matrix based on relatedness
    #                           propto(0 + sample.id | 1, relate_mat_plum_snps),
    #                           data = snps_plum_repro_df,
    #                              family = gaussian)
    #   
    # ## b) model summary  
    #   summary(slc45a3_plum_mod)
    #   
    # ## c) tidy effect estimates
    #   slc45a3_plum_mod_tidy <- 
    #     tidy(slc45a3_plum_mod, 
    #          effects = 'fixed',
    #          conf.int = T,
    #          conf.method = 'Wald',
    #          conf.level = 0.95, 
    #          exponentiate = F) %>%
    #     filter(!(term == '(Intercept)'| term == 'sexm')) %>%
    #     select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
    #              'statistic', 'p.value'))
    #   # NOTE: NO significant snps 
      
    
  ### 3.5 BNC2 effect on breast brightness - both sexes model
      
    ## a) glmmTMB model for effect of BNC2 on breast brightness - female
      bnc2_plum_mod_f <-  glmmTMB(breast.avg.bright ~ 
                                  bnc2.33387499,
                                data = snps_plum_repro_f_df,
                                family = gaussian)
      
    ## b) model summary - female 
      summary(bnc2_plum_mod_f)
      
    ## c) tidy effect estimates - female
      bnc2_plum_mod_tidy_f <- 
        tidy(bnc2_plum_mod_f, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = F) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value')) 
      # NOTE: NO significant snps:
      
    ## d) glmmTMB model for effect of BNC2 on breast brightness - male
      bnc2_plum_mod_m <-  glmmTMB(breast.avg.bright ~ 
                                    bnc2.33387499,
                                  data = snps_plum_repro_m_df,
                                  family = gaussian)
      
    ## b) model summary - male
      summary(bnc2_plum_mod_m)
      
      ## c) tidy effect estimates - male
      bnc2_plum_mod_tidy_m <- 
        tidy(bnc2_plum_mod_m, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = F) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value')) 
      
    # ## a) glmmTMB model for effect of BNC2 on breast brightness
    #   bnc2_plum_mod <-  glmmTMB(breast.avg.bright ~ 
    #                            bnc2.33387499
    #                            # test for intx - no intx
    #                            #* sex +
    #                            + sex  +
    #                      # proportional covariance matrix based on relatedness
    #                            propto(0 + sample.id | 1, relate_mat_plum_snps),
    #                            data = snps_plum_repro_df,
    #                                family = gaussian)
    #   
    # ## b) model summary  
    #   summary(bnc2_plum_mod)
    #   
    # ## c) tidy effect estimates
    #   bnc2_plum_mod_tidy <- 
    #     tidy(bnc2_plum_mod, 
    #          effects = 'fixed',
    #          conf.int = T,
    #          conf.method = 'Wald',
    #          conf.level = 0.95, 
    #          exponentiate = F) %>%
    #     filter(!(term == '(Intercept)'| term == 'sexm')) %>%
    #     select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
    #              'statistic', 'p.value')) 
    #   # NOTE: NO significant snps:
      
      
  ### 3.6 Assess ROR2 snp associations with breast brightness
    ## a) make a list of snp names from ROR2
      ror2_snps <- snps_plum_repro_df %>%
        select(ror2.44196487:ror2.44214826) %>%
        colnames()
     
      ror2_snps <- as.list(ror2_snps)
      
    ## b) Use map to run multiple glmmTMB models - one for each snp - female
      ror2_snps_plum_mods_f <- map(ror2_snps, ~ {
        form <- as.formula(paste('breast.avg.bright ~', .x))
        glmmTMB(
          formula = form,
          data = snps_plum_repro_f_df,
          family = gaussian(),
        )
      })
      
    ## c) Use broom.mixed to tidy effect estimates - female
      ror2_snps_plum_mods_tidy_f <- ror2_snps_plum_mods_f %>%
        map_dfr(tidy, effects = 'fixed',
                conf.int = T,
                conf.method = 'Wald',
                conf.level = 0.95, 
                exponentiate = F) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value'))
      
    ## d) Use map to run multiple glmmTMB models - one for each snp - male
      ror2_snps_plum_mods_m <- map(ror2_snps, ~ {
        form <- as.formula(paste('breast.avg.bright ~', .x))
        glmmTMB(
          formula = form,
          data = snps_plum_repro_m_df,
          family = gaussian(),
        )
      })
      
    ## e) Use broom.mixed to tidy effect estimates - female
      ror2_snps_plum_mods_tidy_m <- ror2_snps_plum_mods_m %>%
        map_dfr(tidy, effects = 'fixed',
                conf.int = T,
                conf.method = 'Wald',
                conf.level = 0.95, 
                exponentiate = F) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value'))
      
      # NOTE: significant snps: 
      # ror2.44206228
      # ror2.44196487
      # ror2.44201599
      
    # ## b) Use map to run multiple glmmTMB models - one for each snp
    #   ror2_snps_plum_mods <- map(ror2_snps, ~ {
    #     # specify model formula including a proportional covariance matrix 
    #     # based on relatedness
    #     # test for sex x snp interactions - NO sig intx p < 0.1
    #     #form <- as.formula(paste('breast.avg.bright ~ sex *', .x,
    #     # if no intx, then control for sex  
    #     form <- as.formula(paste('breast.avg.bright ~ sex +', .x,
    #                       '+ propto(0 + sample.id | 1, relate_mat_plum_snps)'))
    #     glmmTMB(
    #       formula = form,
    #       data = snps_plum_repro_df,
    #       family = gaussian(),
    #     )
    #   })
    #   
    # ## c) Use broom.mixed to tidy effect estimates
    #   ror2_snps_plum_mods_tidy <- ror2_snps_plum_mods %>%
    #     map_dfr(tidy, effects = 'fixed',
    #             conf.int = T,
    #             conf.method = 'Wald',
    #             conf.level = 0.95, 
    #             exponentiate = F) %>%
    #     filter(!(term == '(Intercept)'| term == 'sexm')) %>%
    #     select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
    #              'statistic', 'p.value'))
    #   # NOTE: significant snps: 
    #     # ror2.44206228
    #     # ror2.44196487
    #     # ror2.44201599
      
      
  ### 3.7 APC effect on breast brightness - both sexes model 
    ## a) glmmTMB model for effect of APC on breast brightness - female
      apc_plum_mod_f <-  glmmTMB(breast.avg.bright ~ 
                                 apc.46065597,
                               data = snps_plum_repro_f_df,
                               family = gaussian)
      
    ## b) model summary - female
      summary(apc_plum_mod_f)
      
    ## c) tidy effect estimates - female
      apc_plum_mod_tidy_f <- 
        tidy(apc_plum_mod_f, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = F) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value'))  
      # NOTE: NO significant snps
      
    ## d) glmmTMB model for effect of APC on breast brightness - male
      apc_plum_mod_m <-  glmmTMB(breast.avg.bright ~ 
                                   apc.46065597,
                                 data = snps_plum_repro_m_df,
                                 family = gaussian)
      
    ## e) model summary - male
      summary(apc_plum_mod_m)
      
    ## f) tidy effect estimates - male
      apc_plum_mod_tidy_m <- 
        tidy(apc_plum_mod_m, 
             effects = 'fixed',
             conf.int = T,
             conf.method = 'Wald',
             conf.level = 0.95, 
             exponentiate = F) %>%
        filter(!(term == '(Intercept)')) %>%
        select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
                 'statistic', 'p.value'))  
      
    # ## a) glmmTMB model for effect of APC on breast brightness
    #   apc_plum_mod <-  glmmTMB(breast.avg.bright ~ 
    #                              apc.46065597
    #                            # test for intx - no intx
    #                            #* sex +
    #                            + sex  +
    #                       # proportional covariance matrix based on relatedness
    #                             propto(0 + sample.id | 1, relate_mat_plum_snps),
    #                             data = snps_plum_repro_df,
    #                               family = gaussian)
    #  
    # ## b) model summary  
    #   summary(apc_plum_mod)
    #   
    # ## c) tidy effect estimates
    #   apc_plum_mod_tidy <- 
    #     tidy(apc_plum_mod, 
    #          effects = 'fixed',
    #          conf.int = T,
    #          conf.method = 'Wald',
    #          conf.level = 0.95, 
    #          exponentiate = F) %>%
    #     filter(!(term == '(Intercept)'| term == 'sexm')) %>%
    #     select(c('term', 'estimate', 'std.error', 'conf.low', 'conf.high', 
    #              'statistic', 'p.value'))  
    #   # NOTE: NO significant snps
   
      
  ### 3.8 Make formatted tables of the associations between the each 
      # candidate SNP and plumage color in females and males 
    ## a) row bind female models
      snps_plum_mod_sum_f <- bind_rows(kitlg_snps_plum_mods_tidy_f, 
                                      plxnc1_snps_plum_mods_tidy_f,
                                      slc45a3_plum_mod_tidy_f, 
                                      bnc2_plum_mod_tidy_f,
                                      ror2_snps_plum_mods_tidy_f,
                                      apc_plum_mod_tidy_f) 
      
      
    ## b) row bind male models
      snps_plum_mod_sum_m <- bind_rows(kitlg_snps_plum_mods_tidy_m, 
                                      plxnc1_snps_plum_mods_tidy_m,
                                      slc45a3_plum_mod_tidy_m, 
                                      bnc2_plum_mod_tidy_m,
                                      ror2_snps_plum_mods_tidy_m,
                                      apc_plum_mod_tidy_m) 
      
    ## c) tidy snps_plum_mod_sum_f
      snps_plum_mod_sum_f <- snps_plum_mod_sum_f %>%
        rename(c('se' = 'std.error')) %>%
        select(c(`term`, `estimate`, `se`, `p.value`))
      
    ## d) tidy snps_plum_mod_sum_m
      snps_plum_mod_sum_m <- snps_plum_mod_sum_m %>%
        rename(c('se' = 'std.error')) %>%
        select(c(`term`, `estimate`, `se`, `p.value`))
      
    ## e) format table creating flextable object - female
      snps_plum_mod_sum_frmt_f <- nice_table(snps_plum_mod_sum_f, 
                                            separate.header = T)
      
    ## f) format table creating flextable object - male
      snps_plum_mod_sum_frmt_m <- nice_table(snps_plum_mod_sum_m, 
                                             separate.header = T)
    
    ## g) save snps_plum_mod_sum_frmt_f as a word doc
      flextable::save_as_docx(snps_plum_mod_sum_frmt_f, 
                      path = here('output/snps_plum_mod_sum_frmt_f.docx'))  
      
    ## h) save snps_plum_mod_sum_frmt_m as a word doc
      flextable::save_as_docx(snps_plum_mod_sum_frmt_m, 
                      path = here('output/snps_plum_mod_sum_frmt_m.docx')) 
  
 
    
      
###############################################################################
##############       4. Step 1: Create genetic instrument        ##############
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
     
      
      # Female Plumage snps
        # kitlg.44285279
        # plxnc1.46520400
        # slc45a3.19332900
        # bnc2.33387499
        # ror2.44206228
        # apc.46065597
      
      # Male Plumage snps
        # kitlg.44285485
        # plxnc1.46520543
        # slc45a3.19332900
        # bnc2.33387499
        # ror2.44206228
        # apc.46065597
      
      
      # + plxnc1.46524855
      # + ror2.44206228
      # + ror2.44196487
      
  ### 4.1 PLXNC1 SNPs predictors of breast brightness - males
    ## a) glmmTMB models of PLXNC1 candidate SNPs on breast bright - males
      plxnc1_plum_pred_mods_m <-  glmmTMB(breast.avg.bright ~
                                    + plxnc1.46520543
                                    + plxnc1.46521340 # same
                                    + plxnc1.46521344 # same
                                    + plxnc1.46518352
                                    + plxnc1.46524855
                                    + plxnc1.46525798 
                                    + plxnc1.46526602
                                    + plxnc1.46512949
                                    + plxnc1.46528880 
                                     # +
                                     ,
                        # proportional covariance matrix based on relatedness
                        # propto(0 + sample.id | 1, relate_mat_plum_snps),
                         data = snps_plum_repro_m_df,
                         family = gaussian)
      
    ## b) use dredge from MUMIn package to run all model combinations
      options(na.action = "na.fail")
      aic_plxnc1_mods_m <- dredge(plxnc1_plum_pred_mods_m)
      aic_plxnc1_mods_m <- as.data.frame(aic_plxnc1_mods_m)
      write_csv(aic_plxnc1_mods_m, here('output/aic_plxnc1_mods_m.csv'))
      
    ## c) test independent effects of multiple PLXNC1 proxy SNPs - male
      plxnc1_plum_mods_m <-  glmmTMB(breast.avg.bright ~
                                       + plxnc1.46521340 # same
                                      # + plxnc1.46521344 # same
                                      # + plxnc1.46524855
                                       + plxnc1.46525798 
                                      # + plxnc1.46526602
                                      + plxnc1.46528880 
                                        # +
                                        ,
                      # proportional covariance matrix based on relatedness
                      # propto(0 + sample.id | 1, relate_mat_plum_snps),
                              data = snps_plum_repro_m_df,
                              family = gaussian)
      
    ## d) model summary
      summary(plxnc1_plum_mods_m)
      
      
  ### 4.2  ROR2 SNPs predictors of breast brightness - males
    ## a) glmmTMB models of ROR2 candidate SNPs on breast bright - males
      ror2_plum_pred_mods_m <-  glmmTMB(breast.avg.bright ~
                            + ror2.44196487
                            + ror2.44201599
                            + ror2.44206228
                  # +
                    ,
                  # proportional covariance matrix based on relatedness
                  # propto(0 + sample.id | 1, relate_mat_plum_snps),
                   data = snps_plum_repro_m_df,
                  family = gaussian)
      
    ## b) use dredge from MUMIn package to run all model combinations
      aic_ror2_mods_m <- dredge(ror2_plum_pred_mods_m)
      aic_ror2_mods_m <- as.data.frame(aic_ror2_mods_m)
      write_csv(aic_ror2_mods_m, here('output/aic_ror2_mods_m.csv'))
      
    ## c) test independent effects of multiple ROR2 proxy SNPs - male
      ror2_plum_mods_m <-  glmmTMB(breast.avg.bright ~
                                        + ror2.44196487
                                        + ror2.44206228
                                        # +
                                        ,
                      # proportional covariance matrix based on relatedness
                      # propto(0 + sample.id | 1, relate_mat_plum_snps),
                        data = snps_plum_repro_m_df,
                        family = gaussian)
      
      ## d) model summary
        summary(ror2_plum_mods_m)
        
      
  ## 4.3 Locus specific F-tests
        
      # calculate an F-statistic based on locus specif model and the
        #  snps_plum_pred_reduced_mod sums of squares
        
      ## a) update the name of the locus specific model to snps_plum_pred_mod
        snps_plum_pred_mod <- plxnc1_plum_mods_m
        #snps_plum_pred_mod <- ror2_plum_mods_m
        
      ## b) run a reduced model with no snps
        snps_plum_pred_reduced_mod <-  glmmTMB(breast.avg.bright ~
                                                 + 1,
                      # proportional covariance matrix based on relatedness
                      # propto(0 + sample.id | 1, relate_mat_plum_snps),
                        data = snps_plum_repro_m_df,
                        family = gaussian)
      ## c) Determine the model df from reduced model 
        df_mod_snps_r <- get_df(snps_plum_pred_reduced_mod, type = "model") 
      
      ## d) Determine residual and model degrees of freedom from full model
        df_resid_snps <- get_df(snps_plum_pred_mod) # same as df.residual(model)
        df_mod_snps <- get_df(snps_plum_pred_mod, type = "model") 
        
      ## e) extract residual sums of squares from full model- snps_plum_pred_mod
        resid_snps <- residuals(snps_plum_pred_mod, type = "response")
        resid_ss_snps <- sum(resid_snps^2)
        
      ## f) calculate the Mean Square error from full model - snps_plum_pred_mo
        ms_error_snps = resid_ss_snps/df_resid_snps
        
      ## g) extract full model predicted values
        y_hat_snps <- predict(snps_plum_pred_mod,
                              type = "response")
        
      ## h) extract reduced model predicted values
        y_hat_snps_r <- predict(snps_plum_pred_reduced_mod,
                                type = "response")
        
      ## i) extract actual response values - same for models
        y_snps <- model.response(model.frame(snps_plum_pred_mod))
        
      ## j) Calculate Model Sum of Squares from full model
        mss_snps <- sum((y_hat_snps - mean(y_snps))^2)
        
      ## k) Calculate Model Sum of Squares from reduced model
        mss_snps_r <- sum((y_hat_snps_r - mean(y_snps))^2)
        
      ## l) calculate the extra-model mean square
        ms_model_snps = (mss_snps - mss_snps_r)/(df_mod_snps - df_mod_snps_r)
        
      ## m) calculate partial F-value that is attirbutable to instrument
        f_snps = ms_model_snps/ms_error_snps
        
      ## n) estimate the p-value for the F-value
        p_value_snps <- pf(f_snps, (df_mod_snps - df_mod_snps_r) , df_resid_snps,
                           lower.tail = FALSE)
        
      # ## o) estimate the marginal and conditional R-square
      #   r.squaredGLMM(snps_plum_pred_mod)
        
        
        
        
        
      # install.packages("tidymodels")
      # library(tidymodels)
      # install.packages("glmnet")
      # library(glmnet)
      # 
      # 
      # lasso_spec <- linear_reg(penalty = 0.1, mixture = 1) %>% # mixture = 1 means LASSO
      #   set_engine("glmnet")
      # 
      # # Fit the model using a standard formula interface
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
      
    
        
    
    # 6 create genetic instrument    
RUN AIC model combinations with Plxnc1 and ror2 - find best using aictab or
F
Then add single SNPs loci and compare
pick final isntrument model - get predicted values
create pgs 
regress fitness on both snps and pgs
      
      snps_instrument_mod_m <-  glmmTMB(breast.avg.bright ~
                                      # + kitlg.44285485
                                       #+ plxnc1.46520543
                                     # + plxnc1.46521340 # same
                                     # + plxnc1.46521344 # same
                                     # + plxnc1.46518352
                                     # + plxnc1.46524855
                                     # + plxnc1.46525798 # f-test
                                      + plxnc1.46526602
                                     # + plxnc1.46512949
                                     + plxnc1.46528880 # f-test
                                      # + slc45a3.19332900
                                      # + bnc2.33387499
                                        + ror2.44206228
                                          # + ror2.44201599
                                           + ror2.44196487
                                   #    + apc.46065597
                              # + sex 
                             # +
                             ,
                        # proportional covariance matrix based on relatedness
                              # propto(0 + sample.id | 1, relate_mat_plum_snps),
                               data = snps_plum_repro_m_df,
                               family = gaussian)

    ## b) model summary
      summary(snps_instrument_mod_m)

    ## c) run a reduced model with no snps
      snps_plum_pred_reduced_mod <-  glmmTMB(breast.avg.bright ~
                              #sex 
                                + 1,
                          # proportional covariance matrix based on relatedness
                              #propto(0 + sample.id | 1, relate_mat_plum_snps),
                              data = snps_plum_repro_m_df,
                                family = gaussian)

    ## d) use stats::anova to compare full and reduced model
      snps_mod_comparison <- tibble::tibble(anova(snps_plum_pred_mod,
                                                  snps_plum_pred_reduced_mod))

    ## e) add model variable
      snps_mod_comparison$model <-  c('')

    ## f) move 'model' to first column
      snps_mod_comparison <- snps_mod_comparison %>%
        relocate(model, .before = `Df`)

    ## g) update values to indicate model
      snps_mod_comparison[1,1] <- 'reduced SNPs model'
      snps_mod_comparison[2,1] <- 'full SNPs model'


    ## h) Get fitted values (response scale)
      #snps_plum_pred_val <- fitted.values(snps_plum_pred_mod)
      snps_plum_pred_val <- predict(snps_plum_pred_mod,
                                           type = "response")

    ## i) Combine fitted values with original data
      snps_plum_repro_df$snps.plum.pred.val <- snps_plum_pred_val


  ### 4.2 calculate an F-statistic based on snps_plum_pred_mod and the
      #  snps_plum_pred_reduced_mod sums of squares


    ## a) Determine residual and model degrees of freedom from full model
      df_resid_snps <- get_df(snps_plum_pred_mod) # same as df.residual(model)
      df_mod_snps <- get_df(snps_plum_pred_mod, type = "model") # same as attr(logLik(model), "df")

   ## b) Determine the model df from reduced model
      df_mod_snps_r <- get_df(snps_plum_pred_reduced_mod, type = "model") # same as attr(logLik(model), "df")

    ## c) extract residual sums of squares from full model
      resid_snps <- residuals(snps_plum_pred_mod, type = "response")
      resid_ss_snps <- sum(resid_snps^2)

    ## d) calculate the Mean Square error from full model
      ms_error_snps = resid_ss_snps/df_resid_snps

    ## e) extract full model predicted values
      y_hat_snps <- predict(snps_plum_pred_mod,
                       type = "response")

    ## f) extract reduced model predicted values
      y_hat_snps_r <- predict(snps_plum_pred_reduced_mod,
                              type = "response")

    ## g) extract actual response values - same for models
      y_snps <- model.response(model.frame(snps_plum_pred_mod))

    ## h) Calculate Model Sum of Squares from full model
      mss_snps <- sum((y_hat_snps - mean(y_snps))^2)

    ## i) Calculate Model Sum of Squares from reduced model
      mss_snps_r <- sum((y_hat_snps_r - mean(y_snps))^2)

    ## j) calculate the extra-model mean square
      ms_model_snps = (mss_snps - mss_snps_r)/(df_mod_snps - df_mod_snps_r)

    ## k) calculate partial F-value that is attirbutable to instrument
      f_snps = ms_model_snps/ms_error_snps

    ## l) estimate the p-value for the F-value
      p_value_snps <- pf(f_snps, (df_mod_snps - df_mod_snps_r) , df_resid_snps,
                         lower.tail = FALSE)

    ## m) estimate the marginal and conditional R-square
      r.squaredGLMM(snps_plum_pred_mod)


  ### 4.3 Test if any single snp is a better predictor of breast brightness
    ## a) glmmTMB model for effect of KITLG on breast bright
      sing_snp_plum_pred_mod <-  glmmTMB(breast.avg.bright ~
                                #plxnc1.46524855
                                #ror2.44206228
                                ror2.44196487
                               + sex +
                         # proportional covariance matrix based on relatedness
                               propto(0 + sample.id | 1, relate_mat_plum_snps),
                               data = snps_plum_repro_df,
                                 family = gaussian)

    ## b) model summary
      summary(sing_snp_plum_pred_mod)

    ## c) use stats::anova to compare full and reduced model
      sing_snp_mod_comparison <- tibble::tibble(anova(snps_plum_pred_mod,
                                                   sing_snp_plum_pred_mod))

      
  ### 4.4 Create genetic instrument - polygenic score
    ## a) transform snps to numeric 
      snps_plum_repro_df <- snps_plum_repro_df %>%
        mutate(across(c(kitlg.44277058:apc.46065597), 
                      as.numeric))
      
    ## b) transform snps to numeric - female df
      snps_plum_repro_f_df <- snps_plum_repro_f_df %>%
        mutate(across(c(kitlg.44277058:apc.46065597), 
                      as.numeric))
    
    ## b) transform snps to numeric - male df
      snps_plum_repro_m_df <- snps_plum_repro_m_df %>%
        mutate(across(c(kitlg.44277058:apc.46065597), 
                      as.numeric))
      
  
    ## b) Use imputed allele dosages and effect sizes from (Schield et al 2024)
      # to calculate a polygenic score (PGS) for breast plumage snps
      # PGS = sum(beta1*dosage1 + beta2*dosage2 ...) where beta is a 
      # population level estimate, and dosage is a vector for each individual
      
      # Female Plumage snps
        # kitlg.44285279 beta = -2.37519737
        # plxnc1.46520400 beta = -2.44366682
        # slc45a3.19332900 beta = -3.227701
        # bnc2.33387499 beta = -2.83339322
        # ror2.44206228 beta = -3.5401068
        # apc.46065597 beta = -4.7054829
      
      # Male Plumage snps
        # kitlg.44285485 beta = -4.5346084
        # plxnc1.46520543 beta = -2.13646641
        # slc45a3.19332900 beta = -3.227701
        # bnc2.33387499 beta = -2.83339322
        # ror2.44206228 beta = -3.5401068
        # apc.46065597 beta = -4.7054829

      # constants females
      betas_f <- c(-2.37519737, -2.44366682, -3.227701, -2.83339322, 
                 -3.5401068, -4.7054829)
      
      snps_plum_repro_f_df <- snps_plum_repro_f_df  %>%
        rowwise() %>%
        mutate(pgs.breast = sum(c_across(c(c(kitlg.44285279, plxnc1.46520400,
                                             slc45a3.19332900, bnc2.33387499,
                                             ror2.44206228, apc.46065597))) 
                                * betas_f)) %>%
        ungroup()
      
      # constants males
      betas_m <- c(-4.5346084, -2.13646641, -3.227701, -2.83339322, 
                 -3.5401068, -4.7054829)
      
      
      snps_plum_repro_m_df <- snps_plum_repro_m_df  %>%
        rowwise() %>%
        mutate(pgs.breast = sum(c_across(c(c(kitlg.44285485, plxnc1.46520543,
                                             slc45a3.19332900, bnc2.33387499,
                                             ror2.44206228, apc.46065597))) 
                                * betas_m)) %>%
        ungroup()
    
      
  ### 4.5 KITLG, PLXNC1, and ROR2 snps PGS index effects on breast brightness 
    ## a) glmmTMB model for effect of KITLG, PLXNC1 and ROR2 PGS on 
      # breast brightness
      pgs_plum_pred_mod <-  glmmTMB(breast.avg.bright ~ 
                              pgs.breast +
                              #sex +
                         # proportional covariance matrix based on relatedness
                              propto(0 + sample.id | 1, relate_mat_plum_snps),
                              data = snps_plum_repro_m_df,
                                  family = gaussian)
      
      
    ## b) model summary  
      summary(pgs_plum_pred_mod)
      
    ## c) run a reduced model with no snps 
      pgs_plum_pred_reduced_mod <- glmmTMB(breast.avg.bright ~ 1+ 
                            #sex +
                        # proportional covariance matrix based on relatedness
                            propto(0 + sample.id | 1, relate_mat_plum_snps),
                            data = snps_plum_repro_m_df,
                              family = gaussian)
      
      summary(pgs_plum_pred_reduced_mod)
      
    ## d) use stats::anova to compare full and reduced model  
      pgs_mod_comparison <- tibble::tibble(anova(pgs_plum_pred_mod, 
                                                 pgs_plum_pred_reduced_mod))
      
    ## e) add model variable
      pgs_mod_comparison$model <-  c('')
      
    ## f) move 'model' to first column
      pgs_mod_comparison <- pgs_mod_comparison %>%
        relocate(model, .before = `Df`) 
      
    ## g) update values to indicate model
      pgs_mod_comparison[1,1] <- 'reduced PGS model'
      pgs_mod_comparison[2,1] <- 'full PGS model'
      
    ## h) Get fitted values (response scale)
      pgs_plum_pred_val <- predict(pgs_plum_pred_mod, 
                                   type = "response")
      
    ## i) Combine fitted values with original data
      snps_plum_repro_df$pgs.plum.pred.val <- pgs_plum_pred_val
      
      
  ### 4.6 calculate an F-statistic based on pgs_plum_pred_mod sums of squares
      # F = (model sums squares / model df) / (resid sums squares/ resid df)
      
    ## a) Determine residual and model degrees of freedom
      df_resid_pgs <- get_df(pgs_plum_pred_mod) # same as df.residual(model)
      df_mod_pgs <- get_df(pgs_plum_pred_mod, type = "model") # same as attr(logLik(model), "df")
      
    ## b) extract residual sums of squares  
      resid_pgs <- residuals(pgs_plum_pred_mod, type = "response")
      resid_ss_pgs <- sum(resid_pgs^2)
      
    ## c) calculate the Mean Square error
      ms_error_pgs = resid_ss_pgs/df_resid_pgs
      
    ## d) extract model predicted values
      y_hat_pgs <- predict(pgs_plum_pred_mod, 
                       type = "response")
    
    ## e) extract actual response values  
      y_pgs <- model.response(model.frame(pgs_plum_pred_mod))
      
    ## f) Calculate Model Sum of Squares error
      mss_pgs <- sum((y_hat_pgs - mean(y_pgs))^2)
      ms_model_pgs = mss_pgs/df_mod_pgs
      
    ## g) calculate F-value
      f_pgs = ms_model_pgs/ms_error_pgs
      
    ## h) estimate the p-value for the F-value  
      p_value_pgs <- pf(f_pgs, df_mod_pgs, df_resid_pgs, 
                        lower.tail = FALSE)
      
    ## i) estimate the marginal and conditional R-square
      r.squaredGLMM(pgs_plum_pred_mod)
      
  
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
      
      

    

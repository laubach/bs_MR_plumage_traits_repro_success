################################################################################
#############                                                      #############
#############     Assessing the causal effect of plumage traits    #############
#############     on reproduction using Mendelian randomization    #############   
#############                                                      #############
#############            3. Descriptive stats summaries            #############
#############                                                      #############
#############                   By: Zach Laubach                   #############
#############                 created: 9 Sept 2025                 #############
#############              last updated: 10 June 2026              #############
################################################################################


  ### PURPOSE: Assess descriptive statistics and summaries of the data 
  
  # Code Blocks
    # 1. Configure work space
    # 2. Import data 
    # 3. Descriptive stats traits
    # 4. Summarize snps and relatedness 
    # 5. Descriptive stats reproduction
    # 6. Model natural plumage traits and reproduction
    # 7. Export data
  


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
##############            3. Descriptive stats traits            ##############
############################################################################### 
      
  ### 3.1 Descriptive stats for breast brightness 
    ## a) Total sample size breast brightness
      n.tot.plum <- as.numeric(sum(!is.na(snps_plum_repro_df$sample.id)))
      
    ## b) Total sample size TS   
      n.tot.TS <- as.numeric(sum(!is.na(snps_TS_repro_df$sample.id)))
      
    ## c) Descriptive stats for sex  
      snps_pheno_repro_df %>%
        group_by(sex) %>%
        summarise (n = sum(!is.na(sample.id)))%>%
        mutate(n.plum = n / n.tot.plum) %>%
        ungroup() 
      
     ## d) Descriptive stats for breast.avg.bright  
      breast_sex_sum <- snps_plum_repro_df %>%
        group_by(sex) %>%
        summarise(n.breast.avg.bright = sum(!is.na(breast.avg.bright)),
                avg.breast.avg.bright = format(round(mean(breast.avg.bright, 
                                            na.rm = T), 2), nsmall = 2),
                stdev.breast.avg.bright = format(round(sd(breast.avg.bright, 
                                            na.rm = T), 2), nsmall = 2),
                med.breast.avg.bright = format(round(median(breast.avg.bright,
                                            na.rm = T), 2), nsmall = 2)) %>%
        ungroup() %>%
        mutate(freq.breast.avg.bright = (n.breast.avg.bright / n.tot.plum))
      
    ## e) tidy breast_sex_sum by making a new column 'mean  +/- SD'
      breast_sex_sum$'(mean +/- SD)' = paste(breast_sex_sum$avg.breast.avg.bright, 
                                         breast_sex_sum$stdev.breast.avg.bright, 
                                        sep = ' +/- ')
      
    ## f) tidy breast_sex_sum by removing extra columns and rows
      breast_sex_sum <- breast_sex_sum %>%
        filter(!is.na(sex)) %>%
        rename(`n` = `n.breast.avg.bright`,) %>%
        select(c(n, '(mean +/- SD)')) 
      
    ## g) transpose the data frame for easier viewing
      breast_sex_sum <- as.data.frame(t(breast_sex_sum))
      
    ## h) rename variables
      breast_sex_sum <- breast_sex_sum %>%
        rename(c(`female` = V1,
                 `male` = V2))
      
    ## f) add a variable indicating the summarized measure
      breast_sex_sum$'variable' <- 'breast brightness'
      
    ## g) tidy breast_sex_sum
      breast_sex_sum <- breast_sex_sum %>%
        rownames_to_column() %>%
        relocate('rowname', .before = female) %>%
        relocate('variable', .before = rowname) %>%
        unite(col = ' ', variable, rowname, sep = ' ', remove = T) 
      
      
  ### 3.1 Descriptive stats for right tail streamer length
    ## a) Descriptive stats for rt.intact.avg.ts  
      TS_sex_sum <- snps_TS_repro_df %>%
        group_by(sex) %>%
        summarise(n.rt.TS = sum(!is.na(rt.intact.avg.ts)),
                  avg.rt.TS = format(round(mean(rt.intact.avg.ts, 
                                              na.rm = T), 2), nsmall = 2),
                  stdev.rt.TS = format(round(sd(rt.intact.avg.ts, 
                                              na.rm = T), 2), nsmall = 2),
                  med.rt.TS = format(round(median(rt.intact.avg.ts,
                                              na.rm = T), 2), nsmall = 2)) %>%
        ungroup() %>%
        mutate(freq.rt.TS = (n.rt.TS / n.tot.TS))      

    ## b) tidy TS_sex_sum by making a new column 'mean ± SD'
      TS_sex_sum$'(mean +/- SD)' = paste(TS_sex_sum$avg.rt.TS, 
                                     TS_sex_sum$stdev.rt.TS, 
                                        sep = ' +/- ')
      
    ## c) tidy TS_sex_sum by removing extra columns and rows
      TS_sex_sum <- TS_sex_sum %>%
        filter(!is.na(sex)) %>%
        rename(`n` = `n.rt.TS`,) %>%
        select(c(n, '(mean +/- SD)')) 
      
    ## d) transpose the data frame for easier viewing
      TS_sex_sum <- as.data.frame(t(TS_sex_sum))
      
    ## e) rename variables
      TS_sex_sum <- TS_sex_sum %>%
        rename(c(`female` = V1,
                 `male` = V2))
      
    ## f) add a variable indicating the summarized measure
      TS_sex_sum$'variable' <- 'tail streamer (TS)'
      
    ## g) tidy TS_sex_sum
      TS_sex_sum <- TS_sex_sum %>%
        rownames_to_column() %>%
        relocate('rowname', .before = female) %>%
        relocate('variable', .before = rowname) %>%
        unite(col = ' ', variable, rowname, sep = ' ', remove = T) 
      
      
  
###############################################################################
##############         4. Summarize snps and relatedness         ##############
############################################################################### 
  
  ### 4.1 Check for monomorphic snps
    ## a) assess monorphic breast brightness snps among females
      df_f_plum_snp_summary <- snps_plum_repro_df %>%
        filter(sex == 'f') %>%
        pivot_longer(cols = c(kitlg.44277058:apc.46065597), 
                     names_to = "snp", values_to = "dosage") %>%
        count(snp, dosage)
        
    ## b) create a list of monomorphic breast brightness snps for females
      drop_female_plum <- 
        df_f_plum_snp_summary$snp[df_f_plum_snp_summary$n == 83]
    
    ## c) assess monorphic breast brightness snps among males    
      df_m_plum_snp_summary <- snps_plum_repro_df %>%
        filter(sex == 'm') %>%
        pivot_longer(cols = c(kitlg.44277058:apc.46065597), 
                     names_to = "snp", values_to = "dosage") %>%
        count(snp, dosage)
      
    ## d) create a list of monomorphic breast brightness snps for males  
        drop_male_plum <- 
          df_m_plum_snp_summary$snp[df_m_plum_snp_summary$n == 93]
        
    ## e) union creating a combined list of breast brightness snps that are 
        # monomorphic in either sex
        drop_plum <- dplyr::union(drop_female_plum, drop_male_plum)
        
    ## f) assess monomorphic tail streamer snps for among females   
      df_f_TS_snp_summary <- snps_TS_repro_df %>%
        filter(sex == 'f') %>%
        pivot_longer(cols = c(ice1.98076904:ped1c.100834905), 
                     names_to = "snp", values_to = "dosage") %>%
        count(snp, dosage)
    
    ## g) create a list of monomorphic tail streamer snps for females  
      drop_female_TS <- df_f_TS_snp_summary$snp[df_f_TS_snp_summary$n == 83]
     
    ## h) assess monomorphic tail streamer snps for among males    
      df_m_TS_snp_summary <- snps_TS_repro_df %>%
        filter(sex == 'm') %>%
        pivot_longer(cols = c(ice1.98076904:ped1c.100834905), 
                     names_to = "snp", values_to = "dosage") %>%
        count(snp, dosage)
    
    ## i) create a list of monomorphic tail streamer snps for males       
      drop_male_TS <- df_m_TS_snp_summary$snp[df_m_TS_snp_summary$n == 92]
      
    ## j) union creating a combined list of tail streamer snps that are 
      # monomorphic in either sex  
      drop_TS <- dplyr::union(drop_female_TS, drop_male_TS)
     
    ## k) union creating a combined list of breast brightness and TS snps 
      # that are monomorphic in either sex   
      drop_all <- dplyr::union(drop_plum, drop_TS)
    
    ## l) remove monomorphic snps from snps_pheno_repro_df
      snps_pheno_repro_df <- snps_pheno_repro_df %>%
        select(!c(drop_all))
        
    ## m) remove monomorphic snps from snps_plum_repro_df  
      snps_plum_repro_df <- snps_plum_repro_df %>%
        select(!c(drop_all))
    
    ## n) remove monomorphic snps from snps_plum_repro_f_df      
      snps_plum_repro_f_df <- snps_plum_repro_f_df %>%
        select(!c(drop_all))
        
    ## o) remove monomorphic snps from snps_plum_repro_m_df   
      snps_plum_repro_m_df <- snps_plum_repro_m_df %>%
        select(!c(drop_all))
      
    ## p) remove monomorphic snps from snps_TS_repro_df      
      snps_TS_repro_df <- snps_TS_repro_df %>%
        select(!c(drop_all))
      
    ## q) remove monomorphic snps from snps_TS_repro_f_df    
      snps_plum_repro_m_df 
        snps_TS_repro_f_df <- snps_TS_repro_f_df %>%
          select(!c(drop_all))
      
    ## r) remove monomorphic snps from snps_TS_repro_m_df   
        snps_TS_repro_m_df <- snps_TS_repro_m_df %>%
          select(!c(drop_all))
        
  
  ### 4.2 Calculate genotype frequencies
  
       ## a) update genotype frequencies 
        geno_freq_summary <- snps_pheno_repro_df %>%
          pivot_longer(cols =  c(kitlg.44277058:apc.46065597), 
                       names_to = 'Variable', 
                       values_to = 'Response') %>%
          group_by(Variable, Response) %>%
          summarise(count = n(), .groups = 'drop_last') %>%
          mutate(proportion = count / sum(count)) 
        
      ## b) update genotype frequencies - females
        geno_freq_summary_f <- snps_plum_repro_f_df %>%
          pivot_longer(cols =  c(kitlg.44277058:apc.46065597), 
                       names_to = 'Variable', 
                       values_to = 'Response') %>%
          group_by(Variable, Response) %>%
          summarise(count = n(), .groups = 'drop_last') %>%
          mutate(proportion = count / sum(count)) 
        
      ## c) update genotype frequencies - males
        geno_freq_summary_m <- snps_plum_repro_m_df %>%
          pivot_longer(cols =  c(kitlg.44277058:apc.46065597), 
                       names_to = 'Variable', 
                       values_to = 'Response') %>%
          group_by(Variable, Response) %>%
          summarise(count = n(), .groups = 'drop_last') %>%
          mutate(proportion = count / sum(count)) 

        
  ### 4.3 Tidy and summarize relatedness in the sample population
    ## a) filter relate_coef_df to match individual in snps_plum_repro_df
        relate_coef_plum_df <- relate_coef_df %>%
          filter(ind1 %in% snps_plum_repro_df$sample.id & 
                   ind2 %in% snps_plum_repro_df$sample.id)
        
    ## b) Descriptive stats for reproductive success based on paternity  
      relate_sum <- relate_coef_plum_df %>%
        summarise(n.coef.relate = sum(!is.na(coef.relate)),
                  avg.coef.relate = format(round(mean(coef.relate, 
                                            na.rm = T), 2), nsmall = 2),
                  stdev.coef.relate = format(round(sd(coef.relate, 
                                            na.rm = T), 2), nsmall = 2),
                  med.coef.relate = format(round(median(coef.relate,
                                            na.rm = T), 2), nsmall = 2),
                  min.coef.relate = format(round(min(coef.relate,
                                            na.rm = T), 2), nsmall = 2),
                  max.coef.relate = format(round(max(coef.relate,
                                            na.rm = T), 2), nsmall = 2))%>%
        ungroup()
    
      

###############################################################################
##############        5. Descriptive stats reproduction          ##############
############################################################################### 
     
  ### 5.1 Descriptive stats for reproductive
    ## a) Descriptive stats for total reproductive success based on paternity 
        # from complete nests attempts 1 and 2 
      tot_repro_sex_sum <- snps_pheno_repro_df %>%
        group_by(sex) %>%
        summarise(n.tot.repro = sum(!is.na(tot.repro)),
                  avg.tot.repro = format(round(mean(tot.repro, 
                                              na.rm = T), 2), nsmall = 2),
                  stdev.tot.repro = format(round(sd(tot.repro, 
                                              na.rm = T), 2), nsmall = 2),
                  med.tot.repro = format(round(median(tot.repro,
                                              na.rm = T), 2), nsmall = 2)) %>%
      ungroup() %>%
      mutate(freq.tot.repro = (n.tot.repro / n.tot.plum))
      
    ## b) tidy tot_repro_sex_sum by making a new column 'mean  +/- SD'
      tot_repro_sex_sum$'(mean +/- SD)' = paste(tot_repro_sex_sum$avg.tot.repro, 
                                              tot_repro_sex_sum$stdev.tot.repro, 
                                           sep = ' +/- ')
      
    ## c) tidy tot_repro_sex_sum by removing extra columns and rows
      tot_repro_sex_sum <- tot_repro_sex_sum %>%
        filter(!is.na(sex)) %>%
        rename(`n` = `n.tot.repro`,) %>%
        select(c(n, '(mean +/- SD)')) 
      
    ## d) transpose the data frame for easier viewing
      tot_repro_sex_sum <- as.data.frame(t(tot_repro_sex_sum))
      
    ## e) rename variables
      tot_repro_sex_sum <- tot_repro_sex_sum %>%
        rename(c(`female` = V1,
                 `male` = V2))
      
    ## f) add a variable indicating the summarized measure
      tot_repro_sex_sum$'variable' <- 'total reproduction'
      
    ## g) tidy tot_repro_sex_sum
      tot_repro_sex_sum <- tot_repro_sex_sum %>%
        rownames_to_column() %>%
        relocate('rowname', .before = female) %>%
        relocate('variable', .before = rowname) %>%
        unite(col = ' ', variable, rowname, sep = ' ', remove = T) 
          
    ## h) Descriptive stats for extra pair (EP) reproductive success 
        # based on paternity from complete nests attempts 1 and 2 
      ep_repro_sex_sum <- snps_pheno_repro_df %>%
        group_by(sex) %>%
        summarise(avg.ep.repro = format(round(mean(ep.repro, 
                                                    na.rm = T), 2), 
                                          nsmall = 2),
                  stdev.ep.repro = format(round(sd(ep.repro, 
                                                    na.rm = T), 2), 
                                          nsmall = 2),
                  med.ep.repro = format(round(median(ep.repro,
                                                      na.rm = T), 2), 
                                          nsmall = 2)) %>%
        ungroup()
      
    ## i) tidy ep_repro_sex_sum by making a new column 'mean  +/- SD'
      ep_repro_sex_sum$'(mean +/- SD)' = paste(ep_repro_sex_sum$avg.ep.repro, 
                                             ep_repro_sex_sum$stdev.ep.repro, 
                                              sep = ' +/- ')
      
    ## j) tidy ep_repro_sex_sum by removing extra columns and rows
      ep_repro_sex_sum <- ep_repro_sex_sum %>%
        filter(!is.na(sex)) %>%
        select(c('(mean +/- SD)')) 
      
    ## k) transpose the data frame for easier viewing
      ep_repro_sex_sum <- as.data.frame(t(ep_repro_sex_sum))
      
    ## e) rename variables
      ep_repro_sex_sum <- ep_repro_sex_sum %>%
        rename(c(`female` = V1,
                 `male` = V2))
      
    ## l) add a variable indicating the summarized measure
      ep_repro_sex_sum$'variable' <- 'EP reproduction'
      
    ## m) tidy ep_repro_sex_sum
      ep_repro_sex_sum <- ep_repro_sex_sum %>%
        rownames_to_column() %>%
        relocate('rowname', .before = female) %>%
        relocate('variable', .before = rowname) %>%
        unite(col = ' ', variable, rowname, sep = ' ', remove = T) 
    
    ## n) Descriptive stats for nest 1 and 2 total fecundity (females only)   
      fecund_sex_sum <- snps_pheno_repro_df %>%
        group_by(sex) %>%
        summarise(avg.fecund = format(round(mean(sum.fecundity, 
                                                    na.rm = T), 2), 
                                        nsmall = 2),
                  stdev.fecund = format(round(sd(sum.fecundity, 
                                                    na.rm = T), 2), 
                                        nsmall = 2),
                  med.fecund = format(round(median(sum.fecundity,
                                                      na.rm = T), 2), 
                                        nsmall = 2)) %>%
        ungroup()
      
    ## o) tidy fecund_sex_sum by making a new column 'mean  +/- SD'
      fecund_sex_sum$'(mean +/- SD)' = paste(fecund_sex_sum$avg.fecund, 
                                           fecund_sex_sum$stdev.fecund, 
                                             sep = ' +/- ')
      
    ## p) tidy fecund_sex_sum by removing extra columns and rows
      fecund_sex_sum <- fecund_sex_sum %>%
        filter(!is.na(sex)) %>%
        select(c('(mean +/- SD)')) 
      
    ## q) transpose the data frame for easier viewing
      fecund_sex_sum <- as.data.frame(t(fecund_sex_sum))
      
    ## r) rename variables
      fecund_sex_sum <- fecund_sex_sum %>%
        rename(c(`female` = V1,
                 `male` = V2))
      
    ## s) add a variable indicating the summarized measure
      fecund_sex_sum$'variable' <- 'fecundity'
      
    ## t) tidy fecund_sex_sum
      fecund_sex_sum <- fecund_sex_sum %>%
        rownames_to_column() %>%
        relocate('rowname', .before = female) %>%
        relocate('variable', .before = rowname) %>%
        unite(col = ' ', variable, rowname, sep = ' ', remove = T)
      
    ## u) set male fecundity to NA
      fecund_sex_sum$male <- 'NA'
      
     
  ### 5.2 Combine traits and reproduction and make formatted summary tables
    ## a) row bind breast_sex_sum, TS_sex_sum, tot_repro_sex_sum, 
        # ep_repro_sex_sum, and fecund_sex_sum
      trait_repro_sum <- bind_rows(breast_sex_sum, TS_sex_sum, 
                                   tot_repro_sex_sum, ep_repro_sex_sum, 
                                   fecund_sex_sum) 
      
    ## c) update values to note n
      trait_repro_sum[1,1] <- 'breast brightness (n)'
      trait_repro_sum[3,1] <- 'avg. right TS (n)'
      trait_repro_sum[5,1] <- 'reproduction (n)'
    
    ## c) format table creating flextable object 
      trait_repro_sum_frmt <- nice_table(trait_repro_sum, 
                                               separate.header = T)
      
    ## c) preview the table as a word doc    
      #print(trait_repro_sum_frmt, preview = 'docx')
      
    ## d) save trait_repro_sum_frmt as a word doc
      flextable::save_as_docx(trait_repro_sum_frmt, 
                              path = here('output/trait_repro_sum_frmt.docx'))  
      
      
      

###############################################################################
##############                  6. Export data                   ##############
###############################################################################
      
  ### 6.1 Export data to an RData file 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      
    ## a) Save and export snps_plum_repro data data
      save(file = here('data/5_MR_data.RData'), 
           list = c('snps_plum_repro_df', 'snps_plum_repro_f_df',
                    'snps_plum_repro_m_df', 'geno_freq_summary',
                    'geno_freq_summary_f', 'geno_freq_summary_m'))
   
  

################################################################################
#############                                                      #############
#############     Assessing the causal effect of plumage traits    #############
#############     on reproduction using Mendelian randomization    #############   
#############                                                      #############
#############     1. Tidy snps phenotype and reproduction data     #############
#############                                                      #############
#############                   By: Zach Laubach                   #############
#############                 created: 9 Sept 2025                 #############
#############             last updated: 27 April 2026              #############
################################################################################


  ### PURPOSE: Tidy the 
  
  # Code Blocks
    # 1. Configure work space
    # 2. Import data 
    # 3. Tidy attribute data frame
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
    
    # load lubridate package
      library ('lubridate')

    
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
    ## a) Load phenotype and reproductive data 
      load(here('data/1_pheno_reprod_df.RData'))
      
    ## b) Load imputed snp for candidate loci 
      load(here('data/1_snps_df.RData'))
      
    ## c) Load genome meta data 
      load(here('data/1_genome_meta_df.RData'))

      

###############################################################################
##############          3. Tidy attribute data frame             ##############
############################################################################### 
      
  ### 3.1 Basic tidying pheno_reprod_df
    ## a) Format data to all lower case
      pheno_reprod_df <- AllCharactersToLower(pheno_reprod_df)
      
    ## b) Format variable names
      pheno_reprod_df <- FormatVarNames(pheno_reprod_df)
      
    ## c) Filter adults from pheno_reprod_df
      pheno_reprod_df <- pheno_reprod_df  %>%
        filter(!is.na(year)) %>%
        filter(year >= 2021 & year <= 2023) 
    
    ## d) check for repeats
      repeated_ids <- pheno_reprod_df %>%
        pull(band) %>%
        duplicated() %>%
        {pheno_reprod_df$band[.]}
      
    ## e) remove older repeats and keep most recent sample
      # for the genomics data we removed the older sample for repeats
      pheno_reprod_df <- pheno_reprod_df %>%
        arrange(desc(year)) %>% # Reverses the order
        distinct(band, .keep_all = TRUE)
      
    ## f) create a new column to match sample.id in genomics samples
      pheno_reprod_df$sample.id <- gsub(".*-", "CO_", pheno_reprod_df$band) 
      
    ## g) move sample.id to first column position
      pheno_reprod_df <-  pheno_reprod_df %>% relocate(sample.id)
    
    ## h) remove unused columns
      pheno_reprod_df <- pheno_reprod_df %>%
        select(-c(lts1:notes.1 ,additional.data:mass.2.date)) %>%
        select(-matches(c('initials', 'mb.side', 'old.new', 'hatched',
                          'recap', 'rwl', 'how.sexed', 'body.feathers', 
                          'success', 'julian', 'notes', 'mass', 'manipulated',
                          'manipulation.date', 'chroma', 'hue', 
                          'total.bright')))
        
      
  ### 3.2 Basic tidying imputed snp for candidate loci
    ## a) Format variable names
      snps_df <- FormatVarNames(snps_df)
      
    ## b) Use allele dosages and effect sizes from (Schield et al 2025)
        # to calculate a polygenic risk score (PRS)
        # PRS = sum(beta1*dosage1 + beta2*dosage2 ...) where beta is a 
        # population level estimate, and dosage is a vector for each individual
     
      # Example constants
      betas <- c(2, 2, 2)
      
      snps_df <- snps_df  %>%
        rowwise() %>%
        mutate(prs = sum(c_across(c(kitlg:ror2)) * betas)) %>%
        ungroup()
    
    ## c) Convert snps to categorical variables - homozygous 
        # reference 0|0 as reference category 0. Heterozygous 0|1 or 1|0 = 1,
        # and homozygous alternative allele 1|1 = 2. These values convert to 
        # allele dosages if converted back to numeric.
      snps_df <- transform(snps_df, 
                           kitlg = factor(kitlg, levels = c('0', '1', '2')))
      snps_df <- transform(snps_df, 
                           slc45a2 = factor(slc45a2, levels = c('0', '1', '2')))
      snps_df <- transform(snps_df, 
                           ror2 = factor(ror2, levels = c('0', '1', '2')))
      snps_df <- transform(snps_df, 
                           pde1c = factor(pde1c, levels = c('0', '1', '2')))
    
   
  ### 3.3 Create and tidy snps_pheno_repro_df
    ## a) left join snps_df to pheno_reprod_df 
      snps_pheno_repro_df <- pheno_reprod_df  %>%
        left_join(snps_df,
                  by = c('sample.id' = 'sample.id'), 
                  copy = F) 
    
    ## b) create a variable max.avg.ts that is the longest tail streamer 
        # measure (rts and lts) for each bird
      snps_pheno_repro_df <- snps_pheno_repro_df %>% 
        mutate(max.avg.ts = pmax(mean.rts, mean.lts, na.rm = TRUE))
      
    ## c) create a variable rt.intact.avg.ts that is the average right tail 
      # streamer length measure unless the right 
      snps_pheno_repro_df <- snps_pheno_repro_df %>% 
        mutate(rt.intact.avg.ts = case_when(intact.rts == 'y'
                                                 ~ mean.rts,
                                            intact.rts == 'n'
                                                 ~ mean.lts)) 
      
    # ## c) create a variable ventral.avg.bright that is the mean of the throat
    #   # breast, belly, and vent avg. brightness values
    #   snps_pheno_repro_df <- snps_pheno_repro_df %>% 
    #     mutate(ventral.avg.bright = rowMeans(across(c(throat.avg.bright, 
    #                                                   breast.avg.bright,
    #                                                   belly.avg.bright,
    #                                                   vent.avg.bright)),
    #                                          na.rm = TRUE))
      
      
    ## d) coerce non-numeric values to NA for eggs.1, eggs.2, eggs.3 and 
        # eggs.4 
      snps_pheno_repro_df <- snps_pheno_repro_df %>% 
        mutate(across(c(eggs.1, eggs.2, eggs.3, eggs.4), ~as.numeric(.)))
      
    ## e) create a variable for total fecundity
      snps_pheno_repro_df <- snps_pheno_repro_df  %>%
        rowwise()%>%
        mutate(sum.fecundity = sum(eggs.1, eggs.2, 
                                   #eggs.3, eggs.4, 
                                     na.rm = T))
      
    ## f) coerce non-numeric values to NA for day.12.1, day.12.2, day.12.3, 
      # and day.12.4 fledgling
      snps_pheno_repro_df <- snps_pheno_repro_df %>% 
        mutate(across(c(day.12.1, day.12.2, day.12.3, day.12.4), ~as.numeric(.)))
      
    ## g) create a variable for total fledglings from each adults social pair
        # nest
      snps_pheno_repro_df <- snps_pheno_repro_df  %>%
        rowwise()%>%
        mutate(sum.fledge = sum(day.12.1, day.12.2, 
                                #day.12.3, day.12.4,
                                     na.rm = T))
  
      
  ### 3.4 Subset snps_pheno_repro_df for analyses basesd on missing data 
    ## a) filter samples with no snps data
      snps_pheno_repro_df <- snps_pheno_repro_df %>%
        filter(!is.na(kitlg) | !is.na(slc45a2) | !is.na(ror2) |
                 !is.na(pde1c))  
      
    # ## b) filter samples with no reproductive data
    #   snps_pheno_repro_df <- snps_pheno_repro_df %>%
    #     filter(!is.na(eggs.1) | !is.na(eggs.2) | !is.na(eggs.3) |
    #              !is.na(eggs.4))
    #   
    # ## c) filter samples with no ventral plumage brightness data
    #   snps_pheno_repro_df <- snps_pheno_repro_df %>%
    #     filter(ventral.avg.bright != 'NaN')
      
      

  ### 3.5 Tidy paternity data
    ## a) Format variable names
      paternity_df <- FormatVarNames(paternity_df)
      
    ## b) rename variables
      paternity_df <- paternity_df %>%
        rename(nest = clutch.id.ind2)
      
    ## c) simplify fert.type labels
      paternity_df$fert.type <- gsub("ep.*", "ep", paternity_df$fert.type) 
     
    ## d) create a new columns to match sample.id in genomics samples
      paternity_df$mom.id <- gsub(".*-", "CO_", paternity_df$band.mom) 
      paternity_df$dad.id <- gsub(".*-", "CO_", paternity_df$band.dad) 
     
    ## e) move sample.id to first column position
      paternity_df <- paternity_df %>% relocate(c(mom.id, dad.id))
     
    ## f) filter to include just broods 1-2 and exclude eggs collected
      paternity_df <- paternity_df %>%
        filter(brood == 1 | brood == 2)
      
    ## g) create separate mom data frame
      mom_repro_df <- paternity_df %>%
        select(-c(dad.id, band.dad)) %>%
        filter(!is.na(mom.id))
      
    ## h) create separate dad data frame
      dad_repro_df <- paternity_df %>%
        select(-c(mom.id, band.mom)) %>%
        filter(!is.na(dad.id))
      
    ## i) Calculate female reproductive success (sum of her wp and ep nestlings
        # that survive to ~12 days in broods 1 and 2)
      mom_repro_df <- mom_repro_df%>%
        group_by(mom.id, year) %>%
        summarise(tot.repro = sum(fert, na.rm = TRUE)) %>%
        rename(sample.id = mom.id)  %>%
        ungroup()
      
      
    ## j) Calculate female reproductive success (sum of her wp and ep nestlings
      # that survive to ~12 days in broods 1 and 2)
      dad_repro_df <- dad_repro_df%>%
        group_by(dad.id, year) %>%
        summarise(tot.repro = sum(fert, na.rm = T)) %>%
        rename(sample.id = dad.id) %>%
        ungroup()
      
    ## k) combine mom_repro_df and dad_repro_df
      tot_repro_df <- rbind(mom_repro_df, dad_repro_df)
      
    ## l) left join tot_repro_df to snps_pheno_repro_df
      snps_pheno_repro_df <- snps_pheno_repro_df  %>%
        left_join(tot_repro_df,
                  by = c('sample.id' = 'sample.id', 'year' = 'year'), 
                  copy = F) 
      
    ## m) fill in tot.repro zero according to sum.fledge
      snps_pheno_repro_df$tot.repro <- ifelse(
        is.na(snps_pheno_repro_df$tot.repro) & 
          snps_pheno_repro_df$sum.fledge == 0, 0, snps_pheno_repro_df$tot.repro)
     

  ### 3.6  Basic tidying snps_pheno_repro_df
    ## a) need to reorder the snps_pheno_repro_df by sample.id to match 
      # order of the relate_mat
      # Ascending order using arrange
      snps_pheno_repro_df <- snps_pheno_repro_df %>% 
        arrange(sample.id) 
    
    ## b) filter and subset by missing tot.repro
      snps_pheno_repro_df <- snps_pheno_repro_df %>% 
        filter((!is.na(tot.repro)))
    
    ## c) filter and subset by missing breast brightness plumage
      snps_plum_repro_df <- snps_pheno_repro_df %>% 
        filter((!is.na(breast.avg.bright)))
      
    ## d) filter and subset by misting TS 
      snps_TS_repro_df <- snps_pheno_repro_df %>% 
        filter((!is.na(rt.intact.avg.ts)))
      
      
  ### 3.2 Subset data for downstream analysis
    ## a) make a female breast plumage subset
      snps_plum_repro_f_df <- snps_plum_repro_df %>%
        filter(sex == 'f')
      
    ## b) make a male breast plumage subset
      snps_plum_repro_m_df <- snps_plum_repro_df %>%
        filter(sex == 'm')
      
      
    ## c) make a female tail streamer subset
      snps_TS_repro_f_df <- snps_TS_repro_df %>%
        filter(sex == 'f')
      
    ## d) make a male tail streamer subset
      snps_TS_repro_m_df <- snps_TS_repro_df %>%
        filter(sex == 'm')
      
      
      levels(snps_pheno_repro_df$breast.avg.bright)
      
      

###############################################################################
##############                  4. Export data                   ##############
###############################################################################
      
  ### 4.1 Export data to an RData file 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      
    ## a) Save and export data for CHR 2015 pre-manipulation >= 20 RSSI
      save(file = here('data/2-4_snps_pheno_repro_df.RData'), 
           list = c('snps_pheno_repro_df', 'snps_plum_repro_df',
                    'snps_TS_repro_df', 'snps_plum_repro_f_df', 
                    'snps_plum_repro_m_df', 'snps_TS_repro_f_df',
                    'snps_TS_repro_m_df'))
      
      

    

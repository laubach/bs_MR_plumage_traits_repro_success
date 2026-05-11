################################################################################
#############                                                      #############
#############     Assessing the causal effect of plumage traits    #############
#############     on reproduction using Mendelian randomization    #############   
#############                                                      #############
#############     1. Tidy snps phenotype and reproduction data     #############
#############                                                      #############
#############                   By: Zach Laubach                   #############
#############                 created: 9 Sept 2025                 #############
#############               last updated: 7 May 2026               #############
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
      
    ## i) convert sex to factor
      pheno_reprod_df <- transform(pheno_reprod_df, 
                           sex = factor(sex, levels = c('f', 'm')))
    
    # ## j) create a variable max.avg.ts that is the longest tail streamer 
    #     # measure (rts and lts) for each bird
    #   pheno_reprod_df <- pheno_reprod_df %>% 
    #     mutate(max.avg.ts = pmax(mean.rts, mean.lts, na.rm = TRUE))
      
    ## j) create a variable rt.intact.avg.ts that is the average right tail 
      # streamer length measure unless the right 
      pheno_reprod_df <- pheno_reprod_df %>% 
        mutate(rt.intact.avg.ts = case_when(intact.rts == 'y'
                                                 ~ mean.rts,
                                            intact.rts == 'n'
                                                 ~ mean.lts)) 
      
    # ## k) create a variable ventral.avg.bright that is the mean of the throat
    #   # breast, belly, and vent avg. brightness values
    #   pheno_reprod_df <- pheno_reprod_df %>% 
    #     mutate(ventral.avg.bright = rowMeans(across(c(throat.avg.bright, 
    #                                                   breast.avg.bright,
    #                                                   belly.avg.bright,
    #                                                   vent.avg.bright)),
    #                                          na.rm = TRUE))
      
      
    # ## k) coerce non-numeric values to NA for eggs.1, eggs.2, eggs.3 and 
    #     # eggs.4 
    #   pheno_reprod_df <- pheno_reprod_df %>% 
    #     mutate(across(c(eggs.1, eggs.2, eggs.3, eggs.4), ~as.numeric(.)))
      
    # ## k) create a variable for nest attempts 1&2 fecundity
    #   pheno_reprod_df <- pheno_reprod_df  %>%
    #     rowwise()%>%
    #     mutate(sum.fecundity = sum(eggs.1, eggs.2, 
    #                                #eggs.3, eggs.4, 
    #                                  na.rm = T))
      
    ## k) coerce non-numeric values to NA for day.12.1, day.12.2, day.12.3,
      # and day.12.4 fledgling
      pheno_reprod_df <- pheno_reprod_df %>%
        mutate(across(c(day.12.1, day.12.2, day.12.3, day.12.4), ~as.numeric(.)))
      
    ## l) create a variable for total fledglings from each adults social pair
        # nest
      pheno_reprod_df <- pheno_reprod_df  %>%
        rowwise()%>%
        mutate(sum.fledge = sum(day.12.1, day.12.2,
                                #day.12.3, day.12.4,
                                     na.rm = T))
      
      
  ### 3.2 Tidy paternity data
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
      
    ## l) left join tot_repro_df to pheno_reprod_df
      pheno_reprod_df <- pheno_reprod_df  %>%
        left_join(tot_repro_df,
                  by = c('sample.id' = 'sample.id', 'year' = 'year'), 
                  copy = F) 
      
    ## m) fill in tot.repro zero according to sum.fledge
      pheno_reprod_df$tot.repro <- ifelse(
        is.na(pheno_reprod_df$tot.repro) & 
          pheno_reprod_df$sum.fledge == 0, 0, pheno_reprod_df$tot.repro)
      
      
  ### 3.3 Basic tidying imputed snps for plumage and TS candidate loci
    ## a) Format variable names
      snps_df <- FormatVarNames(snps_df)
      
    ## b) left join snps_df to pheno_reprod_df 
      snps_pheno_repro_df <- pheno_reprod_df  %>%
        left_join(snps_df,
                  by = c('sample.id' = 'sample.id'), 
                  copy = F) 
  
      
  ### 3.4 Subset snps_pheno_repro_df for analyses basesd on missing data 
    ## a) filter samples with no snps data
      snps_pheno_repro_df <- snps_pheno_repro_df %>%
        filter(!is.na(kitlg.44277058))  
      
      # ## b) filter samples with no reproductive data
      #   snps_pheno_repro_df <- snps_pheno_repro_df %>%
      #     filter(!is.na(eggs.1) | !is.na(eggs.2) | !is.na(eggs.3) |
      #              !is.na(eggs.4))
      #   
      # ## c) filter samples with no ventral plumage brightness data
      #   snps_pheno_repro_df <- snps_pheno_repro_df %>%
      #     filter(ventral.avg.bright != 'NaN')
      
     

  ### 3.5  Basic tidying snps_pheno_repro_df
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
      
      
  ### 3.6 Subset data for downstream analysis
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
      
      

###############################################################################
##############                  4. Export data                   ##############
###############################################################################
      
  ### 4.1 Export data to an RData file 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      
    ## a) Save and export data for CHR 2015 pre-manipulation >= 20 RSSI
      save(file = here('data/2-3_snps_pheno_repro_df.RData'), 
           list = c('snps_pheno_repro_df', 'snps_plum_repro_df',
                    'snps_TS_repro_df', 'snps_plum_repro_f_df', 
                    'snps_plum_repro_m_df', 'snps_TS_repro_f_df',
                    'snps_TS_repro_m_df'))
      
      

    

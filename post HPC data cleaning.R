# Mixtures Analysis 
library(purrr)
library(tidyverse)
# Read in data that will be loaded on the HPC 

# Read and restructure results from mixtures analysis performed on HPC

# The issue: Each node writes results to a single row, so data is not rectangular 
# (because each node does not run exactly the same number of models, as 23173 is
# not divisible by 128 nodes). 
source(here::here("!directories.r"))
source(here::here("!set_exposure_outcome_vars.r"))
source(here::here("!load_data.r"))

# SOLAR -------------------------------------------------------------------
# Read in data, without headers.
sol_og_for_colnames <- read.table(fs::path(dir_results, 
                                           "PFAS_mixtures",
                                           "solar_pfas_mixtures_mwas_from_HPC_v3.csv"),
                                  sep = ",", na.strings = "",
                                  as.is = TRUE, fill = TRUE, header = FALSE)

# Get temp col names. Col names are not indexed correctly though- they 
# Need to be shifted to the right by 1, and we need to 
row1 <- sol_og_for_colnames[1,] %>% as.character(.)
col_names <- c("beta",
               row1[-length(row1)], 
               "var.names.9999", 
               "Mean.9999", 
               "SD.9999",         
               "X2.5..9999",      
               "X97.5..9999",    
               "p.val.9999",      
               "metabolite.9999")


# Read in final data with headers
sol_og <- read.table(fs::path(dir_results, 
                                           "PFAS_mixtures",
                                           "solar_pfas_mixtures_mwas_from_HPC_v3.csv"),
                                  sep = ",",col.names = col_names,
                                  na.strings = "",
                                  as.is = TRUE, fill = TRUE, header = TRUE) %>% 
  dplyr::select(-beta)


# Clean Col names
sol_1 <- sol_og %>% 
  janitor::clean_names() %>% 
  rename_all(~str_replace(., "x2_5", "lcl") %>% 
               str_replace(., "x97_5", "ucl") %>% 
               str_replace(., "p_val", "p") %>% 
               str_replace(., "var_names", "var"))

# Get new column names in a dataframe
new_colnames <- tibble(cnms = colnames(sol_1)) %>% 
  mutate(variable = str_split_fixed(cnms, "_", n = 2)[,1], 
         group = str_split_fixed(cnms, "_", n = 2)[,2] %>% 
           if_else(. == "", "0", .) %>% 
           as.numeric)

# create a list of colnames by column group
cnms_bygroup <- split(new_colnames, new_colnames$group)


# create a list of sol result by column groups
sol_2 <- cnms_bygroup %>% 
  modify(~dplyr::select(sol_1, all_of(.$cnms)))

# rename all cols to match names, then cbind
sol_3 <- sol_2 %>% 
  modify(~rename_all(., ~str_split_fixed(., "_", 2)[,1])) %>% 
  bind_rows() %>% 
  filter(!is.na(metabolite)) %>% 
  mutate(var = str_remove(var, ".beta") %>% 
           str_replace("psi", "Mixture effect")) %>% 
  rename(exposure = var, 
         estimate = mean) 


# Calculate new p values
sol_4 <- sol_3 %>% 
  mutate(wald = (estimate/sd),
         p = 2*(1-pnorm(abs(wald),0,1)),
         p_value = pchisq(wald^2, df = 1,lower.tail = FALSE),
         p_value = if_else(wald^2 > 2000, 4.6*(10^ (-256)), p_value),
         neg_log_p = -log10(p_value)) %>% 
  group_by(exposure) %>% 
  mutate(q_value = p.adjust(p_value, method = "fdr"), 
         significance = if_else(p_value < 0.05, "p < 0.05", "Not Sig."), 
         significancefdr = if_else(q_value < 0.05, "q < 0.05", "Not Sig.")) %>% 
  ungroup()


# Join with ft_metadata
sol_final <- left_join(ft_metadata, sol_4, by = c("feature" = "metabolite"))



# Save 
write_csv(sol_final, 
          file = fs::path(dir_results, 
                          'PFAS_Mixtures', 
                          "sol_pfas_mixtures_results_final_v3.csv"))

rm(cnms_bygroup, new_colnames, row1, col_names, 
   sol_1, sol_2,sol_3, sol_4, sol_final, sol_og, sol_og_for_colnames)

# CHS -------------------------------------------------------------------
# Read in data, without headers.
chs_og_for_colnames <- read.table(fs::path(dir_results, 
                                           "PFAS_mixtures",
                                           "chs_pfas_mixtures_mwas_from_HPC_v3.csv"),
                                  sep = ",", na.strings = "",
                                  as.is = TRUE, fill = TRUE, header = FALSE)

# Get temp col names. Col names are not indexed correctly though- they 
# Need to be shifted to the right by 1, and we need to 
row1 <- chs_og_for_colnames[1,] %>% as.character(.)
col_names <- c("beta",
               row1[-length(row1)], 
               "var.names.9999", 
               "Mean.9999", 
               "SD.9999",         
               "X2.5..9999",      
               "X97.5..9999",    
               "p.val.9999",      
               "metabolite.9999")

# Read in final data with headers
chs_og <- read.table(fs::path(dir_results, 
                              "PFAS_mixtures",
                              "chs_pfas_mixtures_mwas_from_HPC_v3.csv"),
                     sep = ",",
                     col.names = col_names,
                     na.strings = "",
                     as.is = TRUE, fill = TRUE, header = TRUE) %>% 
  dplyr::select(-beta)


# Clean Col names
chs_1 <- chs_og %>% 
  janitor::clean_names() %>% 
  rename_all(~str_replace(., "x2_5", "lcl") %>% 
               str_replace(., "x97_5", "ucl") %>% 
               str_replace(., "p_val", "p") %>% 
               str_replace(., "var_names", "var"))

# Get new column names in a dataframe
new_colnames <- tibble(cnms = colnames(chs_1)) %>% 
  mutate(variable = str_split_fixed(cnms, "_", n = 2)[,1], 
         group = str_split_fixed(cnms, "_", n = 2)[,2] %>% 
           if_else(. == "", "0", .) %>% 
           as.numeric)

# create a list of colnames by column group
cnms_bygroup <- split(new_colnames, new_colnames$group)


# create a list of chs result by column groups
chs_2 <- cnms_bygroup %>% 
  modify(~dplyr::select(chs_1, all_of(.$cnms)))

# rename all cols to match names, then cbind
chs_3 <- chs_2 %>% 
  modify(~rename_all(., ~str_split_fixed(., "_", 2)[,1])) %>% 
  bind_rows() %>% 
  filter(!is.na(metabolite))  %>% 
  mutate(var = str_remove(var, ".beta") %>% 
           str_replace("psi", "Mixture effect")) %>% 
  rename(exposure = var, 
         estimate = mean)



# Calculate new p values
chs_4 <- chs_3 %>% 
  mutate(wald = (estimate/sd),
         p = 2*(1-pnorm(abs(wald),0,1)),
         p_value = pchisq(wald^2, df = 1,lower.tail = FALSE),
         p_value = if_else(wald^2 > 2000, 4.6*(10^ (-256)), p_value),
         neg_log_p = -log10(p_value)) %>% 
  group_by(exposure) %>% 
  mutate(q_value = p.adjust(p_value, method = "fdr"), 
         significance = if_else(p_value < 0.05, "p < 0.05", "Not Sig."), 
         significancefdr = if_else(q_value < 0.05, "q < 0.05", "Not Sig.")) %>% 
  ungroup()

# Check the number of unique metabolites (should be 23,176)
length(unique(chs_4$metabolite))


# Join with ft_metadata
chs_final <- left_join(ft_metadata, chs_4, by = c("feature" = "metabolite"))

# Save 
write_csv(chs_final, 
          file = fs::path(dir_results, 
                          'PFAS_Mixtures', 
                          "chs_pfas_mixtures_results_final_v3.csv"))


rm(cnms_bygroup, new_colnames, row1, col_names, 
   chs_1, chs_2,chs_3, chs_4, chs_final, chs_og, chs_og_for_colnames)
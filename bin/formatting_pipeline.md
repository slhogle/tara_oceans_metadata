Tara Oceans Formatting, imputing, and finalizing data
================

  - [INTRO](#intro)
  - [DOWNLOAD AND FORMAT DARWIN DATA](#download-and-format-darwin-data)
      - [Tara and Simons generic
        metadata](#tara-and-simons-generic-metadata)
      - [Synechococcus](#synechococcus)
      - [Prochlorococcus](#prochlorococcus)
      - [Nitrate](#nitrate)
      - [Phosphate](#phosphate)
      - [Dissolved Organic Phosphorus](#dissolved-organic-phosphorus)
      - [Dissolved Total Iron](#dissolved-total-iron)
  - [READ DATA](#read-data)
  - [FORMAT DATA](#format-data)
  - [OUTPUT](#output)

``` r
set.seed(472)

library(here)
library(tidyverse)
library(readxl)
library(janitor)
library(missRanger)
```

# INTRO

We want to combine environmental, biogeographic, biogeochemical, and
hydrographic data into a single dataset that can be loaded with ease.
The final product will be a single R object that we can load and use
right away with statistical models and visualization approaches.

For cases with sparsely missing data (\~25% of observations) we will
impute values using a Random Forest with the missRanger package. The
missRanger package uses the ranger package to do fast missing value
imputation by chained random forest prediction. As such, it serves as an
alternative implementation of the ‘MissForest’ algorithm, see vignette.
In my experience this imputation strategy produces reasonable results.

# DOWNLOAD AND FORMAT DARWIN DATA

Here we use the cmap4r package to get the Darwin model nutrient
climatology. We will use this modeled data since it is highly complete
and represents more of an “time-averaged” picture rather than the point
source observations we have from the GEOTRACES samples.

We will download modeled: - Synechococcus (mmol C/L) =
prokaryote\_c01\_darwin\_clim - Prochlorococcus (mmol C/L) =
prokaryote\_c02\_darwin\_clim - Dissolved Nitrate (mmol N/L) =
NO3\_darwin\_clim - Dissolved Phosphate (mmol P/L) = PO4\_darwin\_clim -
Dissolved Organic Phosphate (mmol P/L) = DOP\_darwin\_clim - Dissolved
Total Fe (mmol Fe/L) = FeT\_darwin\_clim

## Tara and Simons generic metadata

Read this to reduce file size for saving

``` r
mg <- read_tsv(here::here("data", "tara-geotraces-coord-depth.tsv")) 
```

## Synechococcus

``` r
darwinSyn.clim <- get_spacetime(tableName = 'tblDarwin_Plankton_Climatology',
              varName = 'prokaryote_c01_darwin_clim',
              dt1='2008-01-01',
              dt2='2008-12-31',
              lat1=-63,
              lat2=50,
              lon1=-180,
              lon2=180,
              depth1=0,
              depth2=500)

darwinSyn.clim %<>% group_by(lat, lon, depth) %>%
  summarize(value=mean(prokaryote_c01_darwin_clim, na.rm=T)) %>%
  ungroup() %>%
  mutate(variable="prokaryote_c01_darwin_clim")

darwinSyn.clim.joined <- geo_left_join(select(mg, sampleID, lat, lon, depth), darwinSyn.clim, 
                        by = c("lat","lon"),
                        distance_col="geojoin.dist.km",
                        max_dist = 65)

darwinSyn.clim.joined.final <- darwinSyn.clim.joined %>%
  mutate(depth.diff=abs(depth.x-depth.y)) %>%
  group_by(sampleID) %>%
  filter(geojoin.dist.km==min(geojoin.dist.km)) %>%
  filter(depth.diff==min(depth.diff)) %>%
  slice(1) %>%
  ungroup() 

coord.info <- select(darwinPro.clim.joined.final, sampleID, lat=lat.y, lon=lon.y, depth=depth.y, geojoin.dist.km)

saveRDS(coord.info, here::here("data", "darwin_clim", "darwin.coords.depth"))

saveRDS(select(darwinSyn.clim.joined.final, sampleID, value, variable),
        here::here("data", "darwin_clim", "darwinSyn.clim"))
```

## Prochlorococcus

``` r
coord.info <- readRDS(here::here("data", "darwin_clim", "darwin.coords.depth"))

darwinPro.clim <- get_spacetime(tableName = 'tblDarwin_Plankton_Climatology',
              varName = 'prokaryote_c02_darwin_clim',
              dt1='2008-01-01',
              dt2='2008-12-31',
              lat1=-63,
              lat2=50,
              lon1=-180,
              lon2=180,
              depth1=0,
              depth2=500)

darwinPro.clim %<>% group_by(lat, lon, depth) %>%
  summarize(value=mean(prokaryote_c02_darwin_clim, na.rm=T)) %>%
  ungroup() %>%
  mutate(variable="prokaryote_c02_darwin_clim")

darwinPro.clim.joined <- left_join(coord.info, darwinPro.clim) %>%
  select(sampleID, value, variable)

saveRDS(darwinPro.clim.joined, here::here("data", "darwin_clim", "darwinPro.clim"))
```

## Nitrate

``` r
coord.info <- readRDS(here::here("data", "darwin_clim", "darwin.coords.depth"))

darwinNO3.clim <- get_spacetime(tableName = 'tblDarwin_Nutrient_Climatology',
              varName = 'NO3_darwin_clim',
              dt1='2008-01-01',
              dt2='2008-12-31',
              lat1=-63,
              lat2=50,
              lon1=-180,
              lon2=180,
              depth1=0,
              depth2=500)

darwinNO3.clim %<>% group_by(lat, lon, depth) %>%
  summarize(value=mean(NO3_darwin_clim, na.rm=T)) %>%
  ungroup() %>%
  mutate(variable="NO3_darwin_clim")

darwinNO3.clim.joined <- left_join(coord.info, darwinNO3.clim) %>%
  select(sampleID, value, variable)

saveRDS(darwinNO3.clim.joined, here::here("data", "darwin_clim", "darwinNO3.clim"))
```

## Phosphate

``` r
coord.info <- readRDS(here::here("data", "darwin_clim", "darwin.coords.depth"))

darwinPO4.clim <- get_spacetime(tableName = 'tblDarwin_Nutrient_Climatology',
              varName = 'PO4_darwin_clim',
              dt1='2008-01-01',
              dt2='2008-12-31',
              lat1=-63,
              lat2=50,
              lon1=-180,
              lon2=180,
              depth1=0,
              depth2=500)

darwinPO4.clim %<>% group_by(lat, lon, depth) %>%
  summarize(value=mean(PO4_darwin_clim, na.rm=T)) %>%
  ungroup() %>%
  mutate(variable="PO4_darwin_clim")

darwinPO4.clim.joined <- left_join(coord.info, darwinPO4.clim) %>%
  select(sampleID, value, variable)

saveRDS(darwinPO4.clim.joined, here::here("data", "darwin_clim", "darwinPO4.clim"))
```

## Dissolved Organic Phosphorus

``` r
coord.info <- readRDS(here::here("data", "darwin_clim", "darwin.coords.depth"))

darwinDOP.clim <- get_spacetime(tableName = 'tblDarwin_Nutrient_Climatology',
              varName = 'DOP_darwin_clim',
              dt1='2008-01-01',
              dt2='2008-12-31',
              lat1=-63,
              lat2=50,
              lon1=-180,
              lon2=180,
              depth1=0,
              depth2=500)

darwinDOP.clim %<>% group_by(lat, lon, depth) %>%
  summarize(value=mean(DOP_darwin_clim, na.rm=T)) %>%
  ungroup() %>%
  mutate(variable="DOP_darwin_clim")

darwinDOP.clim.joined <- left_join(coord.info, darwinDOP.clim) %>%
  select(sampleID, value, variable)

saveRDS(darwinDOP.clim.joined, here::here("data", "darwin_clim", "darwinDOP.clim"))
```

## Dissolved Total Iron

``` r
coord.info <- readRDS(here::here("data", "darwin_clim", "darwin.coords.depth"))

darwinFe.clim <- get_spacetime(tableName = 'tblDarwin_Nutrient_Climatology',
              varName = 'FeT_darwin_clim',
              dt1='2008-01-01',
              dt2='2008-12-31',
              lat1=-63,
              lat2=50,
              lon1=-180,
              lon2=180,
              depth1=0,
              depth2=500)

darwinFe.clim %<>% group_by(lat, lon, depth) %>%
  summarize(value=mean(FeT_darwin_clim, na.rm=T)) %>%
  ungroup() %>%
  mutate(variable="FeT_darwin_clim")

darwinFe.clim.joined <- left_join(coord.info, darwinFe.clim) %>%
  select(sampleID, value, variable)

saveRDS(darwinFe.clim.joined, here::here("data", "darwin_clim", "darwinFe.clim"))
```

Read in if not already in environment

``` r
darwinSyn.clim <- readRDS(here::here("data", "darwin_clim", "darwinSyn.clim"))
darwinPro.clim <- readRDS(here::here("data", "darwin_clim", "darwinPro.clim"))
darwinNO3.clim <- readRDS(here::here("data", "darwin_clim", "darwinNO3.clim"))
darwinPO4.clim <- readRDS(here::here("data", "darwin_clim", "darwinPO4.clim"))
darwinDOP.clim <- readRDS(here::here("data", "darwin_clim", "darwinDOP.clim"))
darwinFe.clim  <- readRDS(here::here("data", "darwin_clim", "darwinFe.clim"))
```

Bind DARWIN tables and pivot

``` r
darwin.bound.wide <- bind_rows(darwinSyn.clim, darwinPro.clim, darwinNO3.clim,
                          darwinPO4.clim, darwinDOP.clim, darwinFe.clim) %>%
  pivot_wider(names_from = "variable", values_from = "value")
```

``` r
write_tsv(darwin.bound.wide, here::here("data", "darwin2metagenomes.tsv"))
```

# READ DATA

metagenome sequence metadata

``` r
tara.samp <- read_tsv(here::here("data", "tara_metadata.tsv"))
```

Environmental data from the PANGAEA dataset. IMPORTANT: the
ENV-WATERCOLUMN appears to be summary values derived from the
ENV\_DEPTH\_SENSORS product and maybe remote observations (eg MODIS).
The table we use primarly is the DEPTH\_SENSORS table. Also there is a
lot of missing and incomplete nutrient data in the PANGAEA dataset. We
will try to fill in using the supplementary data from Sunugawa 2015
and/or by imputing.

All data was downloaded from here:
<https://doi.pangaea.de/10.1594/PANGAEA.875576>

``` r
tara.meta.a <- read_tsv(here::here("data", "pangea-metadata", "TARA_ENV_DEPTH_SENSORS.tsv")) %>%
  janitor::clean_names() %>%
  janitor::remove_empty(which = "cols") %>%
  dplyr::rename(pangaea_id=sample_id, tara_station=station,
         ebi_sample_id=sample_id_ena) %>%
  mutate(tara_station=ifelse(str_detect(tara_station, "TARA_148"), "TARA_148", tara_station)) %>%
  left_join(select(tara.samp, sampleID, ebi_sample_id), .) %>% 
  distinct()

colnames(tara.meta.a)
```

    ##  [1] "sampleID"                       "ebi_sample_id"                 
    ##  [3] "pangaea_id"                     "sample_id_biosamples"          
    ##  [5] "basis"                          "campaign"                      
    ##  [7] "tara_station"                   "method_device"                 
    ##  [9] "event"                          "date_time"                     
    ## [11] "latitude"                       "longitude"                     
    ## [13] "env_feature"                    "depth"                         
    ## [15] "depth_top"                      "depth_bot"                     
    ## [17] "fraction_lower"                 "fraction_upper"                
    ## [19] "sample_material"                "sample_method"                 
    ## [21] "sample_label"                   "distance_km"                   
    ## [23] "duration"                       "file_name"                     
    ## [25] "nobs"                           "temp"                          
    ## [27] "temp_deg_c_q25"                 "temp_deg_c_q50"                
    ## [29] "temp_deg_c_q75"                 "temp_deg_c_max"                
    ## [31] "cond_m_s_cm_min"                "cond_m_s_cm_q25"               
    ## [33] "cond_m_s_cm_q50"                "cond_m_s_cm_q75"               
    ## [35] "cond_m_s_cm_max"                "sal_min"                       
    ## [37] "sal_q25"                        "sal_q50"                       
    ## [39] "sal_q75"                        "sal_max"                       
    ## [41] "tpot_deg_c_min"                 "tpot_deg_c_q25"                
    ## [43] "tpot_deg_c_q50"                 "tpot_deg_c_q75"                
    ## [45] "tpot_deg_c_max"                 "sigma_theta_kg_m3_min"         
    ## [47] "sigma_theta_kg_m3_q25"          "sigma_theta_kg_m3_q50"         
    ## [49] "sigma_theta_kg_m3_q75"          "sigma_theta_kg_m3_max"         
    ## [51] "chl_a_mg_m3_min"                "chl_a_mg_m3_q25"               
    ## [53] "chl_a_mg_m3_q50"                "chl_a_mg_m3_q75"               
    ## [55] "chl_a_mg_m3_max"                "chl_a_mg_m3_min_1"             
    ## [57] "chl_a_mg_m3_q25_1"              "chl_a_mg_m3_q50_1"             
    ## [59] "chl_a_mg_m3_q75_1"              "chl_a_mg_m3_max_1"             
    ## [61] "chl_a_mg_m3_min_2"              "chl_a_mg_m3_q25_2"             
    ## [63] "chl_a_mg_m3_q50_2"              "chl_a_mg_m3_q75_2"             
    ## [65] "chl_a_mg_m3_max_2"              "par2_day_mol_quanta_m2_day_min"
    ## [67] "par2_day_mol_quanta_m2_day_q50" "par2_day_mol_quanta_m2_day_max"
    ## [69] "par2_perc_min"                  "par2_perc_q50"                 
    ## [71] "par2_perc_max"                  "par3_day_mol_quanta_m2_day_min"
    ## [73] "par3_day_mol_quanta_m2_day_q50" "par3_day_mol_quanta_m2_day_max"
    ## [75] "par3_perc_min"                  "par3_perc_q50"                 
    ## [77] "par3_perc_max"                  "par4_day_mol_quanta_m2_day_min"
    ## [79] "par4_day_mol_quanta_m2_day_q50" "par4_day_mol_quanta_m2_day_max"
    ## [81] "par4_perc_min"                  "par4_perc_q50"                 
    ## [83] "par4_perc_max"

``` r
tara.meta.b <- read_tsv(here::here("data", "pangea-metadata", "TARA_SAMPLES_CONTEXT_ENV-WATERCOLUMN_2.tab"), 
                     skip = 2572) %>%
  janitor::clean_names() %>%
  janitor::remove_empty(which = "cols") %>%
  dplyr::rename(pangaea_id=sample_id_tara_barcode_number_registered_at, tara_station=station_tara_station_number_registered_at,
         ebi_sample_id=sample_id_ena_sample_accession_number) %>%
  mutate(tara_station=ifelse(str_detect(tara_station, "TARA_148"), "TARA_148", tara_station)) %>%
  left_join(select(tara.samp, sampleID, ebi_sample_id), .)

colnames(tara.meta.b)
```

    ##   [1] "sampleID"                                                    
    ##   [2] "ebi_sample_id"                                               
    ##   [3] "pangaea_id"                                                  
    ##   [4] "sample_id_bio_samples_accession_number"                      
    ##   [5] "basis"                                                       
    ##   [6] "campaign"                                                    
    ##   [7] "tara_station"                                                
    ##   [8] "method_device"                                               
    ##   [9] "event"                                                       
    ##  [10] "date_time"                                                   
    ##  [11] "latitude"                                                    
    ##  [12] "longitude"                                                   
    ##  [13] "env_feature_abbreviation_full_name_en"                       
    ##  [14] "depth_nominal_from_which_this_sample_was_co"                 
    ##  [15] "depth_top_m_from_which_this_sample_was_co"                   
    ##  [16] "depth_bot_m_from_which_this_sample_was_co"                   
    ##  [17] "fraction_lower_µm_used_on_board_to_prepare_samp"             
    ##  [18] "fraction_upper_µm_used_on_board_to_prepare_samp"             
    ##  [19] "sample_material_tara_station_number_environmental_f"         
    ##  [20] "sample_method_short_label_describing_the_ta"                 
    ##  [21] "sample_label_tara_event_datetime_station_number"             
    ##  [22] "env_feature_about_the_thickness_of_the_la"                   
    ##  [23] "env_feature_about_the_pelagic_zone_e_g"                      
    ##  [24] "env_feature_about_the_light_environment"                     
    ##  [25] "env_feature_about_the_oxygen_environment"                    
    ##  [26] "env_feature_about_the_nutrient_environmen"                   
    ##  [27] "env_feature_about_the_aggregation_of_phot"                   
    ##  [28] "env_feature_about_vertical_stratification"                   
    ##  [29] "env_feature_about_connectivity_with_the_b"                   
    ##  [30] "env_feature_about_other_information_provi"                   
    ##  [31] "env_feature_about_the_oxygen_environment_1"                  
    ##  [32] "env_feature_about_the_nutrient_environmen_1"                 
    ##  [33] "env_feature_about_the_aggregation_of_phot_1"                 
    ##  [34] "kd490_1_m_at_the_sampling_location_and"                      
    ##  [35] "kd490_1_m_at_the_sampling_location_for"                      
    ##  [36] "kd490_1_m_at_the_sampling_location_for_1"                    
    ##  [37] "kd490_1_m_at_the_sampling_location_for_2"                    
    ##  [38] "kd490_1_m_at_the_sampling_location_for_3"                    
    ##  [39] "bbp470_1_m_443_nm_at_the_sampling_loca"                      
    ##  [40] "ac_cdom_1_m_443_nm_at_the_sampling_loca"                     
    ##  [41] "kd_par_1_m_kd_par_1_is_calculated_for_a"                     
    ##  [42] "kd_par_1_m_kd_par_2_is_calculated_for_a"                     
    ##  [43] "z_eu_m_at_the_sampling_location_for"                         
    ##  [44] "z_eu_m_zeu_0_415_is_the_depth_of_th_1"                       
    ##  [45] "z_eu_m_zeu_0_415_is_the_depth_of_th_2"                       
    ##  [46] "umld_m_based_on_sigma_theta_calcula"                         
    ##  [47] "umld_m_based_on_temperature_calcula"                         
    ##  [48] "d_chl_m_m_calculated_from_in_situ_senso"                     
    ##  [49] "depth_max_brunt_vaisala_freq_m_calculated_from_in_situ_senso"
    ##  [50] "temp_c_at_a_depth_of_10_m_below_the"                         
    ##  [51] "temp_c_at_the_base_of_the_euphotic_z"                        
    ##  [52] "temp_c_at_the_base_of_the_euphotic_z_1"                      
    ##  [53] "temp_c_at_the_depth_of_the_mixed_lay"                        
    ##  [54] "temp_c_at_the_depth_of_the_mixed_lay_1"                      
    ##  [55] "temp_c_at_the_depth_of_maximum_chlor"                        
    ##  [56] "temp_c_at_the_depth_of_maximum_brunt"                        
    ##  [57] "temp_c_at_the_depth_of_maximum_oxyge"                        
    ##  [58] "temp_c_at_the_depth_of_minimum_oxyge"                        
    ##  [59] "temp_c_at_the_depth_of_the_nitraclin"                        
    ##  [60] "sal_at_a_depth_of_10_m_below_the"                            
    ##  [61] "sal_at_the_base_of_the_euphotic_z"                           
    ##  [62] "sal_at_the_base_of_the_euphotic_z_1"                         
    ##  [63] "sal_at_the_depth_of_the_mixed_lay"                           
    ##  [64] "sal_at_the_depth_of_the_mixed_lay_1"                         
    ##  [65] "sal_at_the_depth_of_maximum_chlor"                           
    ##  [66] "sal_at_the_depth_of_maximum_brunt"                           
    ##  [67] "sal_at_the_depth_of_maximum_oxyge"                           
    ##  [68] "sal_at_the_depth_of_minimum_oxyge"                           
    ##  [69] "sal_at_the_depth_of_the_nitraclin"                           
    ##  [70] "sigma_theta_kg_m_3_at_a_depth_of_10_m_below_the"             
    ##  [71] "sigma_theta_kg_m_3_at_the_base_of_the_euphotic_z"            
    ##  [72] "sigma_theta_kg_m_3_at_the_base_of_the_euphotic_z_1"          
    ##  [73] "sigma_theta_kg_m_3_at_the_depth_of_the_mixed_lay"            
    ##  [74] "sigma_theta_kg_m_3_at_the_depth_of_the_mixed_lay_1"          
    ##  [75] "sigma_theta_kg_m_3_at_the_depth_of_maximum_chlor"            
    ##  [76] "sigma_theta_kg_m_3_at_the_depth_of_maximum_brunt"            
    ##  [77] "sigma_theta_kg_m_3_at_the_depth_of_maximum_oxyge"            
    ##  [78] "sigma_theta_kg_m_3_at_the_depth_of_minimum_oxyge"            
    ##  [79] "sigma_theta_kg_m_3_at_the_depth_of_the_nitraclin"            
    ##  [80] "n_2_1_s_2_at_a_depth_of_10_m_below_the"                      
    ##  [81] "n_2_1_s_2_at_the_base_of_the_euphotic_z"                     
    ##  [82] "n_2_1_s_2_at_the_base_of_the_euphotic_z_1"                   
    ##  [83] "n_2_1_s_2_at_the_depth_of_the_mixed_lay"                     
    ##  [84] "n_2_1_s_2_at_the_depth_of_the_mixed_lay_1"                   
    ##  [85] "n_2_1_s_2_at_the_depth_of_maximum_chlor"                     
    ##  [86] "n_2_1_s_2_at_the_depth_of_maximum_brunt"                     
    ##  [87] "n_2_1_s_2_at_the_depth_of_maximum_oxyge"                     
    ##  [88] "n_2_1_s_2_at_the_depth_of_minimum_oxyge"                     
    ##  [89] "n_2_1_s_2_at_the_depth_of_the_nitraclin"                     
    ##  [90] "chl_a_mg_m_3_at_a_depth_of_10_m_below_the"                   
    ##  [91] "chl_a_mg_m_3_at_the_base_of_the_euphotic_z"                  
    ##  [92] "chl_a_mg_m_3_at_the_base_of_the_euphotic_z_1"                
    ##  [93] "chl_a_mg_m_3_at_the_depth_of_the_mixed_lay"                  
    ##  [94] "chl_a_mg_m_3_at_the_depth_of_the_mixed_lay_1"                
    ##  [95] "chl_a_mg_m_3_at_the_depth_of_maximum_chlor"                  
    ##  [96] "chl_a_mg_m_3_at_the_depth_of_maximum_brunt"                  
    ##  [97] "chl_a_mg_m_3_at_the_depth_of_maximum_oxyge"                  
    ##  [98] "chl_a_mg_m_3_at_the_depth_of_minimum_oxyge"                  
    ##  [99] "chl_a_mg_m_3_at_the_depth_of_the_nitraclin"                  
    ## [100] "oxygen_µmol_kg_at_the_depth_of_the_mixed_lay_1"              
    ## [101] "oxygen_µmol_kg_at_the_depth_of_maximum_chlor"                
    ## [102] "oxygen_µmol_kg_at_the_depth_of_maximum_brunt"                
    ## [103] "oxygen_µmol_kg_at_the_depth_of_maximum_oxyge"                
    ## [104] "oxygen_µmol_kg_at_the_depth_of_minimum_oxyge"                
    ## [105] "oxygen_µmol_kg_at_the_depth_of_the_nitraclin"                
    ## [106] "no3_µmol_l_at_the_depth_of_the_mixed_lay_1"                  
    ## [107] "no3_µmol_l_at_the_depth_of_maximum_chlor"                    
    ## [108] "no3_µmol_l_at_the_depth_of_maximum_brunt"                    
    ## [109] "no3_µmol_l_at_the_depth_of_maximum_oxyge"                    
    ## [110] "no3_µmol_l_at_the_depth_of_minimum_oxyge"                    
    ## [111] "no3_µmol_l_at_the_depth_of_the_nitraclin"

Companion table from Sunigawa et al 2015 Science.

``` r
taraW8 <- read_excel(here::here("data", "sunigawa_et_al_CompanionTables.xlsx"), sheet = "Table W8") %>%
  janitor::clean_names() %>%
  janitor::remove_empty(which = "cols") %>%
  select(pangaea_id=pangaea_sample_id, 
         mean_depth_m, mean_temperature_deg_c, mean_salinity_psu, mean_oxygen_umol_kg,
         po4_umol_l, si_umol_l, mean_depth_max_fluo_m)
```

Longhurst codes/identifiers

``` r
longhursts <- read_tsv(here::here("data", "longhurst_codes.tsv"), 
                       col_names = c("sampleID", "lat", "lon", "LongCode", "LongWind", "LongDesc"))
```

Prochlorococcus Ecotype abundances

``` r
ecotype <- read_tsv(here::here("data", "ecotype_read_abundance.tsv")) %>%
  select(sampleID=sample, reads_total, pro_all, pro_HLI, pro_HLII, pro_HLIII_HLIV, pro_LLI, pro_LLII_LLIII, pro_LLIV) %>%
  pivot_longer(pro_HLI:pro_LLIV, names_to = "ecotype", values_to = "count") %>%
  mutate(ra=count/pro_all) %>%
  select(-count) %>%
  pivot_wider(names_from = "ecotype", values_from ="ra") %>%
  mutate(pro=pro_all/reads_total) %>% select(-pro_all, -reads_total)
```

DARWIN modeled data

``` r
darwin.joined.final <- read_tsv(here::here("data", "darwin2metagenomes.tsv")) %>%
  dplyr::rename(synDarwin_umolC.kg=prokaryote_c01_darwin_clim,
                proDarwin_umolC.kg=prokaryote_c02_darwin_clim,
                nitrateDarwin_dissolved_umol.kg=NO3_darwin_clim,
                phosphateDarwin_dissolved_umol.kg=PO4_darwin_clim,
                DOPDarwin_dissolved_umol.kg=DOP_darwin_clim) %>%
  mutate(ironDarwin_dissolved_nmol.kg=FeT_darwin_clim*1000,
         synDarwin_umolC.kg = ifelse(synDarwin_umolC.kg < 0, 0, synDarwin_umolC.kg),
         proDarwin_umolC.kg = ifelse(proDarwin_umolC.kg < 0, 0, proDarwin_umolC.kg)) %>%
  select(-FeT_darwin_clim)
```

# FORMAT DATA

``` r
tara.meta.red.a <- tara.meta.a %>%
  select(sampleID, pangaea_id, tara_station,
         ebi_sample_id, env_feature,
         temp_deg_c_q50,
         sal_q50,
         sigma_theta_kg_m3_q50,
         chl_a_mg_m3_q50_2)

tara.meta.red.b <- tara.meta.b %>%
  select(sampleID, pangaea_id, tara_station,
         ebi_sample_id,
         chlmax_depth=d_chl_m_m_calculated_from_in_situ_senso,
         no3_µmol_l_at_the_depth_of_the_mixed_lay_1,
         no3_µmol_l_at_the_depth_of_maximum_chlor,
         oxygen_µmol_kg_at_the_depth_of_the_mixed_lay_1,
         oxygen_µmol_kg_at_the_depth_of_maximum_chlor) %>%
  group_by(tara_station) %>%
  mutate(chlmax_depth1=mean(chlmax_depth)) %>%
  ungroup()
```

Hmmm… now it appears there is no useable nitrate/NH4/PO4 data in the
PANGEA dataset… Strange because earlier versions had nitrate and
phosphate values (albeit they were often negative). Maybe the quality of
this data is now suspect? The nitrate values from the Sunagawa 2015
supplement seem quite high to me… Will proceed with omitting nitrate
values in Tara oceans and will try and either interpolate them or add
them from the DARWIN model later on.

``` r
tara.samp.meta <- left_join(tara.samp, tara.meta.red.a) %>%
  group_by(tara_station) %>%
  mutate(lat=round(mean(lat), 3), lon=round(mean(lon), 3)) %>%
  ungroup()
```

    ## Joining, by = c("sampleID", "tara_station", "ebi_sample_id")

Final joining

1.  chla\_norm is every value divided by chla at the SCML layer
2.  dcm.layer is if `chla_norm > 0.5*chla_norm`
3.  dcm.max is if `chla_norm >= 0.9*chla_norm`
4.  hdb\_cl is whether the depth of maximum fluorescence falls into the
    basic depth ranges determined for the SCML clusters derived from the
    GEOTRACES dataset

<!-- end list -->

``` r
tara.samp.meta.w8.wide.darwin <- left_join(tara.samp.meta, taraW8, by="pangaea_id") %>%
  left_join(., select(longhursts, -lat, -lon)) %>%
  left_join(., darwin.joined.final) %>%
  left_join(., ecotype) %>%
  left_join(., select(tara.meta.red.b, sampleID, chlmax_depth1)) %>%
  mutate(temperature_dissolved_deg.c = temp_deg_c_q50,
         salinity_dissolved_PSS.1978 = sal_q50,
         sigmatheta_dissolved_kg.m3 = sigma_theta_kg_m3_q50,
         chl_dissolved_mg.kg = chl_a_mg_m3_q50_2/1000/1.025,
         oxygen_dissolved_umol.kg = mean_oxygen_umol_kg,
         phosphate_dissolved_umol.kg = po4_umol_l/1.025,
         silicate_dissolved_umol.kg = si_umol_l/1.025) %>%
  group_by(tara_station, lat, lon) %>%
  mutate(chla_norm=chl_dissolved_mg.kg/max(chl_dissolved_mg.kg), 
         dcm.layer=ifelse(chla_norm>0.5*max(chla_norm), 1, 0),
         dcm.max=ifelse(chla_norm>=0.9*max(chla_norm), 1, 0),
         hdb_cl=case_when(chlmax_depth1 <= 21.4 ~ 1,
                          chlmax_depth1 > 21.4 & chlmax_depth1 <=62.5 ~ 2,
                          chlmax_depth1 > 62.5 & chlmax_depth1 <=84.9 ~ 3,
                          chlmax_depth1 > 84.9 ~ 4)) %>%
  select(sampleID, date, replicate, tara_station, environment, depth, lat, lon,
         ocean, LongWind, LongCode, LongDesc, temperature_dissolved_deg.c,
         salinity_dissolved_PSS.1978, sigmatheta_dissolved_kg.m3, chl_dissolved_mg.kg,
         chla_norm, dcm.layer, dcm.max, hdb_cl, pro_HLI, pro_HLII, pro_HLIII_HLIV, pro_LLI,
         pro_LLII_LLIII, pro_LLIV, pro, oxygen_dissolved_umol.kg, phosphate_dissolved_umol.kg,
         silicate_dissolved_umol.kg, synDarwin_umolC.kg, proDarwin_umolC.kg,
         nitrateDarwin_dissolved_umol.kg, phosphateDarwin_dissolved_umol.kg, 
         DOPDarwin_dissolved_umol.kg, ironDarwin_dissolved_nmol.kg, 
         ebi_run_id, release_version, date_acquired_from_ebi,
         full_sample_name, pangaea_id, env_feature)
```

    ## Joining, by = "sampleID"
    ## Joining, by = "sampleID"
    ## Joining, by = "sampleID"
    ## Joining, by = "sampleID"

# OUTPUT

``` r
write_tsv(tara.samp.meta.w8.wide.darwin, here::here("output", "tara_biogeochem_metadata.tsv"))
```

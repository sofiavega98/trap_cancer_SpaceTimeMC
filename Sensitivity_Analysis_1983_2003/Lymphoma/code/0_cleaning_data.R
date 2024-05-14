## This script cleans SEER lymphoma incidence data 

# Define years of interest
years <- c(1983:2003)

# Define treated year
treated_year <- 1990

# Load packages
library(readr)
library(tidyverse)

# Load cleaned data
SEER <- read.csv("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/data/Lymphoma_11062022.txt")
load("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/data/county_hisp_pop_lymphoma_0_29.RData")

# Recode data
SEER_1 <- SEER %>%
  mutate(cancer_type = recode(`Site.recode.ICD.O.3.WHO.2008...Lymphoma`, "0" = "Lymphoma"),
         state_county_recode = recode(`State.county`,
                                      "0"="San Francisco-Oakland SMSA Registry",
                                      "1"="  CA: Alameda County (06001)",
                                      "2"="  CA: Contra Costa County (06013)",
                                      "3"="  CA: Marin County (06041)",
                                      "4"="  CA: San Francisco County (06075)",
                                      "5"="  CA: San Mateo County (06081)",
                                      "6"="Connecticut Registry",
                                      "7"="  CT: Fairfield County (09001)",
                                      "8"="  CT: Hartford County (09003)",
                                      "9"="  CT: Litchfield County (09005)",
                                      "10"="  CT: Middlesex County (09007)",
                                      "11"="  CT: New Haven County (09009)",
                                      "12"="  CT: New London County (09011)",
                                      "13"="  CT: Tolland County (09013)",
                                      "14"="  CT: Windham County (09015)",
                                      "15"="  CT: Unknown (09999)",
                                      "16"="Atlanta (Metropolitan) Registry",
                                      "17"="  GA: Clayton County (13063)",
                                      "18"="  GA: Cobb County (13067)",
                                      "19"="  GA: DeKalb County (13089)",
                                      "20"="  GA: Fulton County (13121)",
                                      "21"="  GA: Gwinnett County (13135)",
                                      "22"="Hawaii",
                                      "23"="  HI: Hawaii County (15001) - 2000+",
                                      "24"="  HI: Honolulu County (15003) - 2000+",
                                      "25"="  HI: Kalawao County (15005) - 2000+",
                                      "26"="  HI: Kauai County (15007) - 2000+",
                                      "27"="  HI: Maui County (15009) - 2000+",
                                      "28"="  HI: Hawaii County (15911) - Cases Only - pre-2000",
                                      "29"="  HI: Honolulu County (15912) - Cases Only - pre-2000",
                                      "30"="  HI: Kalawao County (15913) - Cases Only - pre-2000",
                                      "31"="  HI: Kauai County (15914) - Cases Only - pre-2000",
                                      "32"="  HI: Maui County (15915) - Cases Only - pre-2000",
                                      "33"="  HI: Hawaii (15900) - Populations Only - pre-2000",
                                      "34"="  HI: Unknown (15999)",
                                      "35"="Iowa Registry",
                                      "36"="  IA: Adair County (19001)",
                                      "37"="  IA: Adams County (19003)",
                                      "38"="  IA: Allamakee County (19005)",
                                      "39"="  IA: Appanoose County (19007)",
                                      "40"="  IA: Audubon County (19009)",
                                      "41"="  IA: Benton County (19011)",
                                      "42"="  IA: Black Hawk County (19013)",
                                      "43"="  IA: Boone County (19015)",
                                      "44"="  IA: Bremer County (19017)",
                                      "45"="  IA: Buchanan County (19019)",
                                      "46"="  IA: Buena Vista County (19021)",
                                      "47"="  IA: Butler County (19023)",
                                      "48"="  IA: Calhoun County (19025)",
                                      "49"="  IA: Carroll County (19027)",
                                      "50"="  IA: Cass County (19029)",
                                      "51"="  IA: Cedar County (19031)",
                                      "52"="  IA: Cerro Gordo County (19033)",
                                      "53"="  IA: Cherokee County (19035)",
                                      "54"="  IA: Chickasaw County (19037)",
                                      "55"="  IA: Clarke County (19039)",
                                      "56"="  IA: Clay County (19041)",
                                      "57"="  IA: Clayton County (19043)",
                                      "58"="  IA: Clinton County (19045)",
                                      "59"="  IA: Crawford County (19047)",
                                      "60"="  IA: Dallas County (19049)",
                                      "61"="  IA: Davis County (19051)",
                                      "62"="  IA: Decatur County (19053)",
                                      "63"="  IA: Delaware County (19055)",
                                      "64"="  IA: Des Moines County (19057)",
                                      "65"="  IA: Dickinson County (19059)",
                                      "66"="  IA: Dubuque County (19061)",
                                      "67"="  IA: Emmet County (19063)",
                                      "68"="  IA: Fayette County (19065)",
                                      "69"="  IA: Floyd County (19067)",
                                      "70"="  IA: Franklin County (19069)",
                                      "71"="  IA: Fremont County (19071)",
                                      "72"="  IA: Greene County (19073)",
                                      "73"="  IA: Grundy County (19075)",
                                      "74"="  IA: Guthrie County (19077)",
                                      "75"="  IA: Hamilton County (19079)",
                                      "76"="  IA: Hancock County (19081)",
                                      "77"="  IA: Hardin County (19083)",
                                      "78"="  IA: Harrison County (19085)",
                                      "79"="  IA: Henry County (19087)",
                                      "80"="  IA: Howard County (19089)",
                                      "81"="  IA: Humboldt County (19091)",
                                      "82"="  IA: Ida County (19093)",
                                      "83"="  IA: Iowa County (19095)",
                                      "84"="  IA: Jackson County (19097)",
                                      "85"="  IA: Jasper County (19099)",
                                      "86"="  IA: Jefferson County (19101)",
                                      "87"="  IA: Johnson County (19103)",
                                      "88"="  IA: Jones County (19105)",
                                      "89"="  IA: Keokuk County (19107)",
                                      "90"="  IA: Kossuth County (19109)",
                                      "91"="  IA: Lee County (19111)",
                                      "92"="  IA: Linn County (19113)",
                                      "93"="  IA: Louisa County (19115)",
                                      "94"="  IA: Lucas County (19117)",
                                      "95"="  IA: Lyon County (19119)",
                                      "96"="  IA: Madison County (19121)",
                                      "97"="  IA: Mahaska County (19123)",
                                      "98"="  IA: Marion County (19125)",
                                      "99"="  IA: Marshall County (19127)",
                                      "100"="  IA: Mills County (19129)",
                                      "101"="  IA: Mitchell County (19131)",
                                      "102"="  IA: Monona County (19133)",
                                      "103"="  IA: Monroe County (19135)",
                                      "104"="  IA: Montgomery County (19137)",
                                      "105"="  IA: Muscatine County (19139)",
                                      "106"="  IA: OBrien County (19141)",
                                      "107"="  IA: Osceola County (19143)",
                                      "108"="  IA: Page County (19145)",
                                      "109"="  IA: Palo Alto County (19147)",
                                      "110"="  IA: Plymouth County (19149)",
                                      "111"="  IA: Pocahontas County (19151)",
                                      "112"="  IA: Polk County (19153)",
                                      "113"="  IA: Pottawattamie County (19155)",
                                      "114"="  IA: Poweshiek County (19157)",
                                      "115"="  IA: Ringgold County (19159)",
                                      "116"="  IA: Sac County (19161)",
                                      "117"="  IA: Scott County (19163)",
                                      "118"="  IA: Shelby County (19165)",
                                      "119"="  IA: Sioux County (19167)",
                                      "120"="  IA: Story County (19169)",
                                      "121"="  IA: Tama County (19171)",
                                      "122"="  IA: Taylor County (19173)",
                                      "123"="  IA: Union County (19175)",
                                      "124"="  IA: Van Buren County (19177)",
                                      "125"="  IA: Wapello County (19179)",
                                      "126"="  IA: Warren County (19181)",
                                      "127"="  IA: Washington County (19183)",
                                      "128"="  IA: Wayne County (19185)",
                                      "129"="  IA: Webster County (19187)",
                                      "130"="  IA: Winnebago County (19189)",
                                      "131"="  IA: Winneshiek County (19191)",
                                      "132"="  IA: Woodbury County (19193)",
                                      "133"="  IA: Worth County (19195)",
                                      "134"="  IA: Wright County (19197)",
                                      "135"="Detroit (Metropolitan) Registry",
                                      "136"="  MI: Macomb County (26099)",
                                      "137"="  MI: Oakland County (26125)",
                                      "138"="  MI: Wayne County (26163)",
                                      "139"="New Mexico Registry",
                                      "140"="  NM: Bernalillo County (35001)",
                                      "141"="  NM: Catron County (35003)",
                                      "142"="  NM: Chaves County (35005)",
                                      "143"="  NM: Colfax County (35007)",
                                      "144"="  NM: Curry County (35009)",
                                      "145"="  NM: De Baca County (35011)",
                                      "146"="  NM: Dona Ana County (35013)",
                                      "147"="  NM: Eddy County (35015)",
                                      "148"="  NM: Grant County (35017)",
                                      "149"="  NM: Guadalupe County (35019)",
                                      "150"="  NM: Harding County (35021)",
                                      "151"="  NM: Hidalgo County (35023)",
                                      "152"="  NM: Lea County (35025)",
                                      "153"="  NM: Lincoln County (35027)",
                                      "154"="  NM: Los Alamos County (35028)",
                                      "155"="  NM: Luna County (35029)",
                                      "156"="  NM: McKinley County (35031)",
                                      "157"="  NM: Mora County (35033)",
                                      "158"="  NM: Otero County (35035)",
                                      "159"="  NM: Quay County (35037)",
                                      "160"="  NM: Rio Arriba County (35039)",
                                      "161"="  NM: Roosevelt County (35041)",
                                      "162"="  NM: Sandoval County (35043)",
                                      "163"="  NM: San Juan County (35045)",
                                      "164"="  NM: San Miguel County (35047)",
                                      "165"="  NM: Santa Fe County (35049)",
                                      "166"="  NM: Sierra County (35051)",
                                      "167"="  NM: Socorro County (35053)",
                                      "168"="  NM: Taos County (35055)",
                                      "169"="  NM: Torrance County (35057)",
                                      "170"="  NM: Union County (35059)",
                                      "171"="  NM: Cibola/Valencia",
                                      "172"="    NM: Cibola County (35006) - 1982+",
                                      "173"="    NM: Valencia County (35061) - 1982+",
                                      "174"="    NM: Cibola/Valencia (35910) - pre-1982",
                                      "175"="  NM: Unknown (35999)",
                                      "176"="Utah Registry",
                                      "177"="  UT: Beaver County (49001)",
                                      "178"="  UT: Box Elder County (49003)",
                                      "179"="  UT: Cache County (49005)",
                                      "180"="  UT: Carbon County (49007)",
                                      "181"="  UT: Daggett County (49009)",
                                      "182"="  UT: Davis County (49011)",
                                      "183"="  UT: Duchesne County (49013)",
                                      "184"="  UT: Emery County (49015)",
                                      "185"="  UT: Garfield County (49017)",
                                      "186"="  UT: Grand County (49019)",
                                      "187"="  UT: Iron County (49021)",
                                      "188"="  UT: Juab County (49023)",
                                      "189"="  UT: Kane County (49025)",
                                      "190"="  UT: Millard County (49027)",
                                      "191"="  UT: Morgan County (49029)",
                                      "192"="  UT: Piute County (49031)",
                                      "193"="  UT: Rich County (49033)",
                                      "194"="  UT: Salt Lake County (49035)",
                                      "195"="  UT: San Juan County (49037)",
                                      "196"="  UT: Sanpete County (49039)",
                                      "197"="  UT: Sevier County (49041)",
                                      "198"="  UT: Summit County (49043)",
                                      "199"="  UT: Tooele County (49045)",
                                      "200"="  UT: Uintah County (49047)",
                                      "201"="  UT: Utah County (49049)",
                                      "202"="  UT: Wasatch County (49051)",
                                      "203"="  UT: Washington County (49053)",
                                      "204"="  UT: Wayne County (49055)",
                                      "205"="  UT: Weber County (49057)",
                                      "206"="  UT: Unknown (49999)",
                                      "207"="Seattle (Puget Sound) Registry",
                                      "208"="  WA: Clallam County (53009)",
                                      "209"="  WA: Grays Harbor County (53027)",
                                      "210"="  WA: Island County (53029)",
                                      "211"="  WA: Jefferson County (53031)",
                                      "212"="  WA: King County (53033)",
                                      "213"="  WA: Kitsap County (53035)",
                                      "214"="  WA: Mason County (53045)",
                                      "215"="  WA: Pierce County (53053)",
                                      "216"="  WA: San Juan County (53055)",
                                      "217"="  WA: Skagit County (53057)",
                                      "218"="  WA: Snohomish County (53061)",
                                      "219"="  WA: Thurston County (53067)",
                                      "220"="  WA: Whatcom County (53073)"),
         fips = str_extract(state_county_recode, "([:digit:]{5})"),
         year_recode = case_when(`Year.of.diagnosis` == 0 ~ NA_real_,
                                 TRUE ~ `Year.of.diagnosis` + 1974)
  )


# Subset to be ages less than 29
SEER_final <- SEER_1 %>% select(fips, year_recode, `Age.recode.with..1.year.olds`, Count) %>%
  filter(`Age.recode.with..1.year.olds` <= 6,
         !is.na(fips),
         !is.na(year_recode)) %>%
  mutate(age_recode = recode(`Age.recode.with..1.year.olds`,
                             "0"="00 years",
                             "1"="01-04 years",
                             "2"="05-09 years",
                             "3"="10-14 years",
                             "4"="15-19 years",
                             "5"="20-24 years",
                             "6"="25-29 years")
  ) %>%
  select(fips, year_recode, age_recode, Count) %>%
  group_by(fips, year_recode) %>%
  summarise(count = sum(Count)) %>%
  rename(year = year_recode) %>%
  ungroup()

# Rename columns
colnames(SEER_final) <- c("FIPS","YEAR_DX","CL_CASES") 

# Subset SEER to start at start year and end at end year
data_full <- SEER_final
data_full <- data_full %>% filter(!YEAR_DX < years[1] & !YEAR_DX > years[length(years)])

# For Hawaii counties, combine appropriate FIPS
HI_pre2000 <- c(15911, 15912, 15913, 15914, 15915)
HI_post2000 <- c(15001, 15003, 15005, 15007, 15009)

for (i in 1:length(HI_pre2000)){
  # Find the indices prior to 2000 in post HI FIPS to be replaced
  ind_rep_post2000 <- which((data_full$FIPS == HI_post2000[i]) & (data_full$YEAR_DX < 2000))
  # Find the indices prior 2000 in the pre HI FIPS to replace 
  ind_rep_pre2000 <- which((data_full$FIPS == HI_pre2000[i]) & (data_full$YEAR_DX < 2000))
  
  print(data_full[ind_rep_pre2000,2:3])
  print(data_full[ind_rep_post2000,2:3])
  
  # Replace every column but FIPS column (this keeps current FIPS)
  data_full[ind_rep_post2000,2:3] <- data_full[ind_rep_pre2000,2:3]
  print(data_full[ind_rep_post2000,2:3])
}

# Remove old HI counties 
data_full <- data_full %>% filter(!(FIPS %in% HI_pre2000))

## Set Up Data

# Specify treated CT and CA Counties
data_full <- data_full %>% mutate(
  # Specify CT as treated after 2003 (should be this but check)
  C = case_when(
    str_starts(FIPS, "09") & YEAR_DX >= treated_year ~ 1,
    TRUE ~ 0),
  # Specify CA as treated after 2003 (should be this but check)
  C = case_when(
    str_starts(FIPS, "06") & YEAR_DX >= treated_year ~ 1, 
    TRUE ~ C)
)

# Add Hispanic Data replace population data with Census data

## dont have census data from 1983 to 1989 so copy 1990
data_full_hisp <- left_join(data_full,county_hisp_pop, by = c("FIPS" = "fips", "YEAR_DX" = "year"))

data_full_hisp <- data_full_hisp %>%
  fill(c(pct_hispanic_pop, total_pop), .direction = "up") 

# Assign column names
colnames(data_full_hisp) <- c("FIPS","YEAR_DX","CL_CASES","C","pct_hispanic_pop",
                              "POP")
# FIPS 15005 was incorporated into 15009 (HI) I'm going to remove
data_full_hisp <- data_full_hisp %>% filter(!(FIPS==15005))

# Remove unknown counties from the dataset
data_full_hisp<- data_full_hisp %>% filter(! (FIPS %in% c("09999", "15999", "35999", "49999")) )

# Identify if a row (county) has all 0s
data_full_hisp %>%
  group_by(FIPS) %>%
  summarise(total_cases = sum(CL_CASES)) %>%
  filter(total_cases == 0) %>%
  pull(FIPS) -> zero_case_FIPS

# Remove counties with all 0s
data_full_hisp_noZero <- data_full_hisp %>%
  filter(! FIPS %in% zero_case_FIPS)

# Take out counties if there is only 1 case or no cases in the pre-treatment years
data_full_hisp_noZero %>%
  group_by(FIPS) %>%
  summarise(total_cases = sum(CL_CASES[which(YEAR_DX %in% seq(years[1],(trt_year-1)))])) %>%
  filter(total_cases <= 1) %>%
  pull(FIPS) -> zero_case_FIPS
data_full_hisp_noZero <- data_full_hisp_noZero %>%
  filter(! FIPS %in% zero_case_FIPS)

# Save final data
save(data_full_hisp_noZero ,file = paste0('/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/data/Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData'))

##########################
## Get adjacency matrix ##
##########################

# Load full adjacency matrix
load("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/data/us_adj_mat.RData")

# Subset full adjacency matrix based on data_full_hisp_noZero
partial_us_adj_mat <- us_adj_mat[rownames(us_adj_mat)  %in%  
                                        unique(data_full_hisp_noZero$FIPS),colnames(us_adj_mat)  %in%  unique(data_full_hisp_noZero$FIPS)]

# Re-order columns and rows to be in order of the fips in the dataset
partial_us_adj_mat <- partial_us_adj_mat[match(unique(data_full_hisp_noZero$FIPS), colnames(partial_us_adj_mat)),
                                         match(unique(data_full_hisp_noZero$FIPS), colnames(partial_us_adj_mat))]

# Save adjacency matrix
save(partial_us_adj_mat,file=paste0('/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/data/lymphoma_adj_mat_',years[1],'_',years[length(years)],'_0-29.RData'))



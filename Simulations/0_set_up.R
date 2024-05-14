# This script sets up the simulations
# For the paper, we use New Mexico as a base for the simulations
# This script identifies the NM counties in the SEER data,
# finds the partial adjacency matrix, and identifies 6 as treated

trt_year <- 1996

# Load packages
#library(tidycensus)
#library(spdep)
source("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/functions/Functions.R")

# Load data
load("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Simulations/data/data_full.RData")
data <- data_full
rm(data_full)

# Pull New Mexico
NM <- data[which(startsWith(data$FIPS,"35")),]

# Remove county 35999 as it doesn't really exist
NM_new <- subset(NM, NM$FIPS !="35999")

# Load adjacency matrices
load("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/data/partial_W_nm.RData") # NM spatial adj matrix
load("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Simulations/data/temp_adj_matrix_nm.RData") # NM temporal adj matrix

### Artificially block out time and counties as treated ###
# Choosing artificial treated counties
trt_fips <- c("35001", "35005", "35006", "35007", "35009", "35013")

# Assigning a variable indicating whether treated or not
# The treated counties are treated at trt_year and after
NM_new$C[which((NM_new$YEAR_DX>=trt_year) & (NM_new$FIPS %in% trt_fips))] <- 1

# Get population matrix
pop_mat <- get_pop_mat(NM_new)

#-------------------------------------------------------------------------------
# Get adjacency matrix
#nm <- get_acs(geography = "county", 
#              variables = c(medincome = "B19013_001"), 
#              state = "NM",
#              geometry = TRUE,
#              year = 2010)
#W_nm <-nb2mat(poly2nb(nm, queen=TRUE),style='B') #full adjacency matrix

#since we don't have every county in New Mexico we need to get a partial adj matrix.
#partial_W_nm <- W_nm[which(nm$GEOID  %in%  
#                             unique(NM$FIPS)),which(nm$GEOID  %in%  unique(NM$FIPS))] #partial adjacency matrix

# Save adjacency matrix
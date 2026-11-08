# Jan Blackburn
# Math 257 Multivariate Analysis Project
# Code for data staging and cleaning
# Initial code generated using Claude based on prompts and modified as needed

########### House Cleaning

rm(list=ls())


########### Linking census and congressional roster data 

# =============================================================================
# Math 257 Project 2: build the 2024 House district analysis dataset
#
# Output: one row per U.S. House district (435), 119th Congress lines
#   Outcomes:   turnout (House votes cast / CVAP), dshare (D two-party share)
#   Predictors: race/ethnicity of CVAP, log density, % bachelor's+, median
#               household income, % of adults 18-29 and 65+, incumbency
#
# Data sources
#   Local files:  1976-2024-house.xlsx (MIT Election Data + Science Lab)
#                 Geo_Data.xlsx        (2024 Census Gazetteer, 119th CDs)
#   tidycensus:   ACS 5-year 2020-2024 detailed tables B01003, B05003 (+ B, D,
#                 H, I iterations), B15003, B19013, B01001
# =============================================================================

library(tidyverse)
library(readxl)
library(tidycensus)

if (Sys.getenv("CENSUS_API_KEY") == "") {
  stop("No Census API key found. Run census_api_key('YOUR_KEY', install = TRUE) once in the console.")
}

house_file <- "1976-2024-house.xlsx"
gaz_file   <- "Geo_Data.xlsx"
acs_year   <- 2024                                  # = ACS 5-year 2020-2024

# -----------------------------------------------------------------------------
# 1. Census data via tidycensus
# -----------------------------------------------------------------------------

# Citizen voting-age population (CVAP) = citizens 18+ = native + naturalized,
# male (009, 011) and female (020, 022). The same line numbers apply to the
# race/ethnicity iterations of B05003:
#   ""  = all     B = Black alone     D = Asian alone
#   H = White alone, not Hispanic     I = Hispanic or Latino
cvap_lines <- c("009", "011", "020", "022")
cvap_groups <- c(all = "", black = "B", asian = "D", white_nh = "H", hisp = "I")

cvap_vars <- imap(cvap_groups, ~ setNames(paste0("B05003", .x, "_", cvap_lines),
                                          paste0("cvap_", .y, "_", seq_along(cvap_lines)))) %>%
  unname() %>% unlist()

acs_vars <- c(
  tot_pop  = "B01003_001",                         # total population
  adult_m  = "B05003_008", adult_f = "B05003_019", # population 18+
  edu_tot  = "B15003_001",                         # population 25+
  edu_ba   = "B15003_022", edu_ma = "B15003_023",
  edu_prof = "B15003_024", edu_phd = "B15003_025",
  med_inc  = "B19013_001",                         # median household income
  lf_civ   = "B23025_003",   # civilian labor force
  unemp    = "B23025_005",   # unemployed
  setNames(sprintf("B01001_%03d", c(7:11, 31:35)),  paste0("a1829_", 1:10)),
  setNames(sprintf("B01001_%03d", c(20:25, 44:49)), paste0("a65_",   1:12)),
  hh_tot   = "B19001_001",
  setNames(sprintf("B19001_%03d", 2:5),   paste0("hh_lo_", 1:4)),   # under $25k
  setNames(sprintf("B19001_%03d", 16:17), paste0("hh_hi_", 1:2)),   # $150k and over
  gini     = "B19083_001",
  pov_tot  = "B17001_001", pov_below = "B17001_002",
  cvap_vars
)

acs_raw <- get_acs(geography = "congressional district",
                   variables = acs_vars,
                   year = acs_year, survey = "acs5",
                   output = "wide")

census <- acs_raw %>%
  rename(key = GEOID) %>%
  filter(!str_detect(key, "ZZ$|98$")) %>%   # drop undefined areas, DC, PR
  transmute(
    key,
    tot_pop      = tot_popE,
    adults       = adult_mE + adult_fE,
    cvap         = rowSums(across(matches("^cvap_all_\\d+E$"))),
    pct_black    = rowSums(across(matches("^cvap_black_\\d+E$")))    / cvap,
    pct_asian    = rowSums(across(matches("^cvap_asian_\\d+E$")))    / cvap,
    pct_white_nh = rowSums(across(matches("^cvap_white_nh_\\d+E$"))) / cvap,
    pct_hisp     = rowSums(across(matches("^cvap_hisp_\\d+E$")))     / cvap,
    pct_ba       = (edu_baE + edu_maE + edu_profE + edu_phdE) / edu_totE,
    med_inc      = med_incE,
    pct_18_29    = rowSums(across(matches("^a1829_\\d+E$"))) / adults,
    pct_female   = (cvap_all_3E + cvap_all_4E) / cvap,
    pct_low_inc  = rowSums(across(matches("^hh_lo_\\d+E$"))) / hh_totE,
    pct_high_inc = rowSums(across(matches("^hh_hi_\\d+E$"))) / hh_totE,
    pct_unemp = unempE / lf_civE,
    gini         = giniE,
    pct_poverty  = pov_belowE / pov_totE,
    pct_65plus   = rowSums(across(matches("^a65_\\d+E$")))   / adults
  )

stopifnot(nrow(census) == 435)

# -----------------------------------------------------------------------------
# 2. Gazetteer: land area -> population density (local file)
# -----------------------------------------------------------------------------

geo <- read_excel(gaz_file) %>%
  rename_with(str_trim) %>%                                # INTPTLONG has trailing spaces
  mutate(key = str_pad(as.character(GEOID), 4, pad = "0")) %>%  # Excel dropped leading 0
  select(key, ALAND_SQMI, INTPTLAT, INTPTLONG)

# -----------------------------------------------------------------------------
# 3. MIT House returns (local file): 2022 and 2024 general elections
# -----------------------------------------------------------------------------

# Rows that are not votes for a candidate (blank, over-, under-votes, void,
# Maine ranked-choice exhausted ballots). Write-ins and "SCATTERING" are kept:
# they are real votes.
nonvote_rows <- c("BLANK", "BLANKS", "BLANK VOTES", "OVER VOTES", "OVERVOTES",
                  "UNDER VOTES", "UNDERVOTES", "VOID",
                  "EXHAUSTED BALLOT", "EXHAUSTED BALLOTS")

house_raw <- read_excel(house_file, na = c("", "NA")) %>%   # older rows store "NA" as text
  filter(year %in% c(2022, 2024), stage == "GEN", !special,
         state_po != "DC") %>%
  filter(!candidate %in% nonvote_rows) %>%
  mutate(party = recode(party, "DEMOCRATIC-FARMER-LABOR" = "DEMOCRAT"),
         key   = paste0(str_pad(state_fips, 2, pad = "0"),
                        str_pad(district,   2, pad = "0")))

# One row per candidate per district: sum fusion lines, main party = largest line
cands <- house_raw %>%
  group_by(year, state_po, district, key, candidate) %>%
  summarise(votes      = sum(candidatevotes),
            party      = party[which.max(candidatevotes)],
            totalvotes = first(totalvotes),
            .groups = "drop")

# ---- 3a. 2024 outcomes ------------------------------------------------------

house24 <- cands %>%
  filter(year == 2024) %>%
  group_by(key, state_po, district) %>%
  summarise(dem        = sum(votes[party %in% "DEMOCRAT"]),
            rep        = sum(votes[party %in% "REPUBLICAN"]),
            votes_cast = sum(votes),          # consistent across states
            totalvotes = first(totalvotes),   # -1 = race not on ballot (FL-20, OK-3)
            .groups = "drop") %>%
  mutate(votes_cast = if_else(totalvotes < 0, NA_real_, votes_cast),
         contested  = dem > 0 & rep > 0)

stopifnot(nrow(house24) == 435)

# ---- 3b. Incumbency ---------------------------------------------------------
# Incumbent = a 2024 candidate who won in 2022 (same state, party, surname,
# first initial), plus members who entered through 2023-24 special elections,
# which the MIT file does not contain. Matching by state rather than district
# handles incumbents who changed districts after 2024 redistricting.

name_parts <- function(x) {
  tokens <- x %>%
    str_to_upper() %>%
    str_remove_all('"[^"]*"') %>%            # quoted nicknames
    str_replace_all("[^A-Z\\- ]", " ") %>%   # punctuation and garbled accents
    str_split("\\s+")
  tokens <- map(tokens, ~ .x[nchar(.x) > 1 & !.x %in% c("JR", "SR", "II", "III", "IV")])
  tibble(first = map_chr(tokens, ~ if (length(.x)) .x[1] else NA_character_),
         sur   = map_chr(tokens, ~ if (length(.x)) .x[length(.x)] else NA_character_))
}

cands <- bind_cols(cands, name_parts(cands$candidate)) %>%
  mutate(init = str_sub(first, 1, 1))

winners22 <- cands %>%
  filter(year == 2022) %>%
  group_by(key) %>%
  slice_max(votes, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(state_po, party, sur, init, district22 = district, winner22 = candidate)

inc_auto <- cands %>%
  filter(year == 2024) %>%
  inner_join(winners22, by = c("state_po", "party", "sur", "init")) %>%
  group_by(winner22) %>%                      # if a 2022 winner matches two 2024
  filter(n() == 1 | district == district22) %>%   # candidates, keep same district
  ungroup() %>%
  select(key, candidate, party)

# Manual additions (state, district, surname as it appears in the 2024 file)
inc_manual <- tribble(
  ~state_po, ~district, ~sur_pattern,   ~reason,
  "CA", 37, "KAMLAGER",  "Name changed to Kamlager-Dove",
  "IL",  4, "GARC",      "Accented name garbled in file",
  "VA",  4, "MCCLELLAN", "Special election Feb 2023",
  "RI",  1, "AMO",       "Special election Nov 2023",
  "UT",  2, "MALOY",     "Special election Nov 2023",
  "NY",  3, "SUOZZI",    "Special election Feb 2024",
  "NY", 26, "KENNEDY",   "Special election Apr 2024",
  "CA", 20, "FONG",      "Special election May 2024",
  "OH",  6, "RULLI",     "Special election Jun 2024",
  "NJ", 10, "MCIVER",    "Special election Sep 2024"
)

inc_extra <- cands %>%
  filter(year == 2024) %>%
  inner_join(inc_manual, by = c("state_po", "district"),
             relationship = "many-to-many") %>%
  filter(str_detect(str_to_upper(candidate), sur_pattern),
         party %in% c("DEMOCRAT", "REPUBLICAN")) %>%
  select(key, candidate, party)

stopifnot(nrow(inc_extra) == nrow(inc_manual))   # every manual fix found

incumbents <- bind_rows(inc_auto, inc_extra) %>% distinct(key, .keep_all = TRUE)

# +1 = Republican incumbent running, -1 = Democratic incumbent, 0 = open seat
incumbency <- house24 %>%
  select(key) %>%
  left_join(incumbents %>% select(key, inc_party = party), by = "key") %>%
  mutate(inc = case_when(inc_party == "REPUBLICAN" ~  1,
                         inc_party == "DEMOCRAT"   ~ -1,
                         TRUE                      ~  0),
         inc_dem = as.integer(inc == -1),
         inc_rep = as.integer(inc ==  1))

# -----------------------------------------------------------------------------
# 4. Combine and create analysis variables
# -----------------------------------------------------------------------------

analysis <- house24 %>%
  left_join(incumbency %>% select(key, inc, inc_dem, inc_rep), by = "key") %>%
  left_join(census, by = "key") %>%
  left_join(geo,    by = "key") %>%
  mutate(turnout     = votes_cast / cvap,
         dshare      = if_else(contested, dem / (dem + rep), NA_real_),
         density     = tot_pop / ALAND_SQMI,
         log_density = log10(density),
         inc_10k     = med_inc / 10000)

stopifnot(nrow(analysis) == 435,
          !anyNA(analysis$cvap), !anyNA(analysis$ALAND_SQMI))

# Modeling sample: contested races with a valid vote count (expect 396)
model_data <- analysis %>% filter(contested, !is.na(turnout))

# -----------------------------------------------------------------------------
# 5. Checks
# -----------------------------------------------------------------------------

cat("Districts:", nrow(analysis), "| contested:", sum(analysis$contested),
    "| modeling sample:", nrow(model_data), "\n")
print(count(analysis, inc))                       # expect 382 incumbents in total
print(summary(select(analysis, turnout, dshare, log_density, pct_ba,
                     med_inc, pct_18_29, pct_65plus, pct_black, pct_hisp)))

write_csv(analysis, "analysis_2024.csv")

View(analysis)

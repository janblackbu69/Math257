# Math257: Project 1
This project looks at what drives turnout and Democratic two-party vote share in the 2024 U.S. House elections, using the 119th Congress district lines. 
It also compares separate univariate regressions with a single multivariate multiple regression that models both outcomes jointly. 
The analysis uses 396 contested districts; the 39 uncontested races were excluded. 
The predictors cover race and ethnicity, education, income and inequality, age, sex, population density, and incumbency. 
The data come from three public sources:

MIT Election Data and Science Lab. 2017. "U.S. House 1976–2024." Harvard Dataverse.
https://doi.org/10.7910/DVN/IG0UN2

This dataset is released under the Creative Commons CC0 1.0 Universal Public Domain
Dedication (https://creativecommons.org/publicdomain/zero/1.0/). CC0 waives copyright, so
the data may be reused without permission. CC0 does not require attribution, but the
dataset is cited here as standard scholarly practice.

US Census Data:
- American Community Survey (ACS) 5-Year Estimates, 2020–2024, detailed tables B01001,
  B01003, B05003 (with B05003B, D, H and I), B15003, B17001, B19001, B19013, B19083 and
  B23025, for 119th Congress districts. Retrieved through the Census Data API with the
  tidycensus R package.
- 2024 Gazetteer Files: 119th Congressional Districts (land area).
  https://www.census.gov/geographies/reference-files/time-series/geo/gazetteer-files.html
- Citizen Voting Age Population (CVAP) Special Tabulation, 2020–2024, and its technical
  documentation.

Census Bureau data products are works of the U.S. federal government. They are not subject
to copyright in the United States (17 U.S.C. § 105) and may be freely reused and
redistributed.

- Data_Extraction.R extracts and formats the data storing the data in 'analysis_2024.csv'.
- Data_Analysis.R performs modeling and generates the exhibits stored in the figures and output folders.
- Draft reports are stored in Report_Draft_v# for major revisions.
- Jan_Blackburn_Report_Final_qmd is the final version of the submitted report.

LLM Use: Claude Opus 5.5 was used in developing the initial R-code for data extraction, data analysis and formatting of exhibits in Quarto.

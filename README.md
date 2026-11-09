# Math257: Project 1
This project looks at what drives turnout and Democratic two-party vote share in the 2024 U.S. House elections, using the 119th Congress district lines. 
It also compares separate univariate regressions with a single multivariate multiple regression that models both outcomes jointly. 
The analysis uses 396 contested districts; the 39 uncontested races were excluded. 
The predictors cover race and ethnicity, education, income and inequality, age, sex, population density, and incumbency. 
The data come from three public sources:
  https://doi.org/10.7910/DVN/IG0UN2
  https://www.census.gov/programs-surveys/acs
  https://www.census.gov/data/developers.html
  https://www.census.gov/geographies/reference-files/time-series/geo/gazetteer-files.html

- Data_Extraction.R extracts and formats the data storing the data in 'analysis_2024.csv'.
- Data_Analysis.R performs modeling and generates the exhibits stored in the figures and output folders.
- Draft reports are stored in Report_Draft_v# for major revisions.
- Jan_Blackburn_Report_Final_qmd is the final version of the submitted report.

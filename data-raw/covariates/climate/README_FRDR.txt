This README.txt file was generated on 2021-03-15 by Frédéric Cyr

--------------------
GENERAL INFORMATION
--------------------

1. Title of Dataset: Newfoundland and Labrador climate index

2. Author Information
	A. Principal Investigator Contact Information
		Name: Frédéric Cyr
		Institution: Fisheries and Oceans Canada
		Email: Frederic.Cyr@dfo-mpo.gc.ca

3. Date of data collection: 
1951-01-01/2020-12-31

4. Geographic location of data collection: 
Newfoundland and Labrador, Northwest Atlantic Ocean

5. This dataset is regularly updated, please make sure to access the most recent version.

---------------------------
SHARING/ACCESS INFORMATION
---------------------------

1. Licenses/restrictions placed on the data: 
These data are available under a CC BY 4.0 license <https://creativecommons.org/licenses/by/4.0/> 

2. Links to publications that cite or use the data:
Cyr, F. and Galbraith, P. S.: A climate index for the Newfoundland and Labrador shelf, Earth Syst. Sci. Data Discuss. [preprint], https://doi.org/10.5194/essd-2020-350, in review, 2020.

3. Recommended citation for this dataset: 
Cyr, Frédéric, Galbraith, Peter S. (2020). Newfoundland and Labrador climate index. Federated Research Data Repository. doi:10.20383/101.0301

---------------------
DATA & FILE OVERVIEW
---------------------

1. File List

   A. Filename: NL_climate_index.csv
      Short description: The Newfoundland and Labrador (NL) climate index aims to describe the environmental conditions on the NL shelf and in the Northwest Atlantic as a whole. It consists of the arithmetic average of 10 annual normalized anomalies (or subindices).     

   B. Filename: NL_climate_index_all_fields.csv
      Short description: Individual subindices used to construct the Newfoundland and Labrador climate index. Each subindex is normalized using the 1991-2020 climatological period. Some signs were reversed so that positive anomalies are indicative of warmer conditions (e.g. sea ice)

   C. Filename: NL_climate_index_natural_signs.csv
      Short description: Individual subindices of the Newfoundland and Labrador climate index in their natural sign. Each subindex is normalized using the 1991-2020 climatological period.

-----------------------------------------------------------------
DATA-SPECIFIC INFORMATION FOR: NL_climate_index_all_fields.csv
-----------------------------------------------------------------

1. Number of variables: 11

2. Number of cases/rows: 69

3. Missing data codes:
        Code/symbol        Definition

4. Variable List:
   	Year
	Winter NAO	
	Air Temp	
	Sea Ice	Icebergs	
	SST	
	S27 Temp	
	S27 Sal	
	S27 CIL	
	CIL area	
	Bottom Temp

    A. Name: Year 
       Description: The year.

    B. Name: Winter NAO
       Description: Average North Atlantic Oscillation over the months of December to March (sign reversed).

    C. Name: Air Temp
       Description: Mean normalized anomalies of annual air temperature at Nuuk (Greenland), Iqaluit (Baffin Island), Cartwright (Labrador), Bonavista (Newfoundland) and St. John's (Newfoundland). 

    D. Name: Sea Ice
       Description: Mean normalized anomalies of Sea ice maximum area and season duration for Northern Labrador, Southern Labrador and Newfoundland shelves (sign reversed).
    
    E. Name: Icebergs
       Description: Normalized anomalies of the number of icebergs crossing 48degN on the Grand Banks (sign reversed).

    F. Name: SST
       Description: Mean normalized anomalies of Sea Surface Temperature over NAFO divisions 2HJ3KLNOP.

    G. Name: S27 Temp
       Description: Normalized anomalies of the vertically-averaged temperature at Station 27.

    H. Name: S27 Sal
       Description: Normalized anomalies of the vertically-averaged salinity at Station 27.

    I. Name: S27 CIL
       Description: Normalized anomalies of the summer (June-August) cold intermediate layer core temperature at Station 27.

    J. Name: CIL area
       Description: Mean normalized anomalies of the summer cold intermediate layer area over hydrographic sections Seal Island, Bonavista and Flemish Cap on the Newfoundland and Labrador shelf (sign reversed).

    K. Name: Bottom Temp
       Description: Mean normalized anomalies of the bottom temperature during spring (NAFO divisions 3LNOPs) and fall (NAFO divisions 2HJ3KLNO). 

-----------------------------------------------------------------
DATA-SPECIFIC INFORMATION FOR: NL_climate_index_all_fields_natural_signs.csv
-----------------------------------------------------------------

1. Number of variables: 11

2. Number of cases/rows: 69

3. Missing data codes:
        Code/symbol        Definition
	
4. Variable List:
   	Year
	Winter NAO	
	Air Temp	
	Sea Ice	Icebergs	
	SST	
	S27 Temp	
	S27 Sal	
	S27 CIL	
	CIL area	
	Bottom Temp

    A. Name: Year 
       Description: The year.

    B. Name: Winter NAO
       Description: Average North Atlantic Oscillation over the months of December to March.

    C. Name: Air Temp
       Description: Mean normalized anomalies of annual air temperature at Nuuk (Greenland), Iqaluit (Baffin Island), Cartwright (Labrador), Bonavista (Newfoundland) and St. John's (Newfoundland). 

    D. Name: Sea Ice
       Description: Mean normalized anomalies of Sea ice maximum area and season duration for Northern Labrador, Southern Labrador and Newfoundland shelves.
    
    E. Name: Icebergs
       Description: Normalized anomalies of the number of icebergs crossing 48degN on the Grand Banks.

    F. Name: SST
       Description: Mean normalized anomalies of Sea Surface Temperature over NAFO divisions 2HJ3KLNOP.

    G. Name: S27 Temp
       Description: Normalized anomalies of the vertically-averaged temperature at Station 27.

    H. Name: S27 Sal
       Description: Normalized anomalies of the vertically-averaged salinity at Station 27.

    I. Name: S27 CIL
       Description: Normalized anomalies of the summer (June-August) cold intermediate layer core temperature at Station 27.

    J. Name: CIL area
       Description: Mean normalized anomalies of the summer cold intermediate layer area over hydrographic sections Seal Island, Bonavista and Flemish Cap on the Newfoundland and Labrador shelf.

    K. Name: Bottom Temp
       Description: Mean normalized anomalies of the bottom temperature during spring (NAFO divisions 3LNOPs) and fall (NAFO divisions 2HJ3KLNO). 

-----------------------------------------------------------------
DATA-SPECIFIC INFORMATION FOR: NL_climate_index.csv
-----------------------------------------------------------------

1. Number of variables: 2

2. Number of cases/rows: 69

3. Missing data codes:
        Code/symbol        Definition

4. Variable List:
   	Year
	Climate index
	
    A. Name: Year
       Description: The year. 

    B. Name: Climate index.
       Description: Average of the 10 sub-indices of the NL climate index (found in file NL_climate_index_all_fields.csv).

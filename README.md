Sealice_Project

Task 1: Epidemiological study of sea lice on Altantic salmon farms in Ireland

Questions:
1 Effect of environmental factors on sea lice abundance
2 Effect of geographical factors on sea lice abundance
3 Effect of temporal factors on sea lice abundance

1 Effect of environmental factors on sea lice abundance
Using MI Individual Database 
Looking at the effect of environmental factors (temperature, salinity, oxygen) on sea lice abundance
Using glmmTMB
Response variable: count of total mobile sea lice
Random effect: Site identity
Fixed effect: temperature, oxygen, salinity, year, time in seawater, class of fish (S0/S1/2),
company (three levels: Mowi, Bradan, others), pen (Random, Standard) 
Detection of temporal autocorrelation

2 Effect of geographical factors on sea lice abundance
Using MI Individual Database 
Looking at the effect of geographical factors (Depth, Distance from mouth, Length, Orientation, Width of mouth) on sea lice abundance
Using glmmTMB
Response variable: count of total mobile sea lice
Random effect: Pen identity
Fixed effect: Depth, Distance from mouth, Length, Orientation, Width of mouth, Year, Month, Time in sea, Class, Company, Fallow period
Detection of temporal autocorrelation

3 Effect of temporal factors on sea lice abundance
Using MI Mean Database 
Looking at the effect of temporaal factors (Year, Month) on sea lice abundance
Using glmmTMB
Response variable: mean of total mobile sea lice per pen
Random effect: Site identity
Fixed effect: Year, Month, Class, Time in sea, Pen, Bay
Detection of temporal autocorrelation

Exploration of temporal autocorrelation in datasets with decreased variation
Mi Individual Database: exploration of data on the level of fish, pen, farm (count per fish, mean per pen, mean per farm)
Mi Mean Database: exploration of data on the level of pen, farm (mean per pen, mean per farm)
There is decrease of temporal autocorrelation and its complexity with the simplification of the data. It seems that the temporal autocorrelation is caused by the complex structure of the data (30 fish, 2 pens per farm from one data point)




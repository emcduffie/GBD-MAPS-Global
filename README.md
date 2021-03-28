# GBD-MAPS-Global
Analysis scripts for the 2021 Global Burden of Disease - Major Air Pollution Sources (GBD-MAPS) - Global project.

The Global Burden of Disease-Major Air Pollution Sources (GBD-MAPS) – Global project, funded by the Health Effects Institute, is a joint collaboration led by researchers at Washington University in St. Louis, the University of British Columbia, the University of Washington, and Dalhousie University. Expanding upon a similar approach used in two previous GBD-MAPS studies for India and China, the GBD-MAPS-Global project employs an integrated modeling approach to identify dominant sources of ambient PM<sub>2.5</sub> pollution and to quantify the associated health impacts at the national level, for all 204 countries and territories currently included in the GBD project. 

The analysis and results are described in detail in the following manuscript:
McDuffie, E. E., Martin, R. V., Spadaro, J. V., Burnett, R., Smith, S. J., O'Rourke, P., Hammer, M., van Donkelaar, A., Bindle, L., Shah, V., Jaegle, L., Luo, G., Yu, F., Adeniran, J., Lin, J., Brauer, M. Source Sector and Fuel Contributions to Ambient PM2.5 Attributable Mortality Across Multiple Spatial Scales, In Review

Additional GBD-MAPS Project Details are at: https://sites.wustl.edu/acag/datasets/gbd-maps/

# Overview

All code is written in MATLAB. 

## Input Files

Required input files are included in the Input_Data dataset (available for download at Zenodo at: 10.5281/zenodo.4642700) or are described in the headers of each script file. 

*NOTE: these scripts do not contain the raw GEOS-Chem model results.* Due to the size of the GEOS-Chem results dataset (~1 Tb), please contact erin.mcduffie (@) wustl.edu for inquiries about access to the raw GEOS-Chem model output dataset. 

The GEOS-Chem source code used to produce the results is avialble here: https://github.com/emcduffie/GC_v12.1.0_EEM. 

The emission input dataset used for the sensitivity simulations is avialable for download at Zenodo: https://zenodo.org/record/3754964
and is described in detail in McDuffie et al., ESSD 2020: https://essd.copernicus.org/articles/12/3413/2020/essd-12-3413-2020.html

## General Description 
All analysis scripts are called from the 'GBD_MAPS_Global_Analysis_Main.m' file and are called in order from A to E. 

A - Make_EmisSensSim_ResultsFileFn.m - reads in raw GEOS-Chem output files (monthly resolution) from each sensitivity simulation and calculates the PM<sub>2.5</sub> concentrations from each (ug/m<sup>-3</sup>)

A1 - Merge_GlobalNestedFn - sub-function called from A that merges gridded model results from global and nested simulations onto a 0.5x0.5 degree global grid

B - Calc_FractionalSourceContributions_GridFn - function that takes the PM<sub>2.5</sub> concentrations from each sensitivity simulation and calculates the fractional contributions (on a grid) from each emission source

C - Calc_AbsFrac_SourceContribution_RegionCountryListFn - calculates the absolute PM<sub>2.5</sub> concentrations from each emission source (ug/m<sup>-3</sup>) at the regional and national levels and then uses these to calculate the (population-weighted) absolute and fractional source contributions

C1 - Calc_AbsFrac_SourceContribution_UrbanListFn - sub-function called from C that also calculates the absolute PM<sub>2.5</sub> concentrations from each emission (ug/m<sup>-3</sup>) source and the population-weighted absolute and fractional source contribution at the sub-national level 

D - Calc_PMmort_SourceAttFn - uses the high resolution downscaled PM<sub>2.5</sub> exposure estimates (provided as input from Input_Data/DownscaledPM.zip) to calculate source-specific PM<sub>2.5</sub> attributable mortality estiamtes for 6 PM<sub>2.5</sub>-disease pairs and two neonatal disorders

D1 - GBD_MortalityCalculation_Tool_v1 - sub-function called from D that uses the relative risk curves from the Global Burden of Disease (GBD) and Global Exposure Mortality Model (GEMM), as well as national-level baseline mortality data from the GBD (from Input_Data/DiseaseBurden.zip), and high resolution PM<sub>2.5</sub> exposure estiamtes to calculate the PM<sub>2.5</sub> attributable disease burden. (Tool originally written and provided in MS Excel by Joseph V. Spadaro)

E - WriteFinalTableFn - function that writes the population-weighted fraction emission source contributions, annual mean PM<sub>2.5</sub> exposure level, and attributable PM<sub>2.5</sub> mortality estiamtes for 21 regions (+ global), 204 countries, and 200 sub-national areas (these data are reported in Supplemental Data Files 1 and 2 in McDuffie et al., 2021) 

The Main Analysis script can be run for the years 2017 or 2019, but requires GEOS-Chem model result data as input, in the correct format. 

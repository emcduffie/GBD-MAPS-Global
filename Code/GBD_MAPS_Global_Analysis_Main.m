%GBD_MAPS_Global_Analysis_Main
%Main analysis script for the GBD-MAPS GLOBAL Study
% Erin McDuffie - Last Updated: Nov. 8, 2020
%
%Description
% Main analysis script that follows three steps
% 1) Load and merge global and nested GEOS-Chem results files
% 2) Calculate the fractional gridded source contributions to PM2.5 mass
% 3) Calculate the regional, national, and urban-area frcational and
%       absolute contributions, as well as the associated
%       population-weighted PM2.5 mass and population
% 4) Calculate the disease burden associated with each source
% 5) Write results to final data file
% Details:
% Code can be run to calculate the sectoral contributions (w/ or w/out
% select fuel-sector sources (e.g., industry coal) OR can be run to
% calculate the total fuel contributions
%%%%%%%%%%%%
%INPUTS:
% Directory, specifies the location of raw the GEOS-Chem simulations results (monthly
%       average netCDF files)
% final_data_location_root (string), specify the root location of where to
%       put the output files
% load_resuts_data, (= 0 or 1), set to 0 if you have not yet run the script. 
%       Set to 1 if you have already run the script and would like to load 
%       merged results from 'outFile' rather than re-loading the GEOS-Chem
%       data (used to save computational time). Default set to 0.
% skip_frac, (= 0 or 1), set to 0 if you have not yet run the script. Set
%       to 1 if you have already run the script and would like to load
%       previously-calculated results (gridded fractional source
%       contributions,population-weighted results, PM2.5 concentrations, and
%       population) from the country_data_, region_data_, and
%       urban_data_outfiles (used to save computational time).Default set
%       to 0.
% use_fuel, (=0 or 1), Set to 0 to calculate sectoral results only. Set to
%       1 to calculate sectoral results with feul-specific cases (e.g.,
%       industry coal, residential biofuel). Default is set to 1. Ignored
%       if 'fuel_cases' set to 1.
% fuel_cases, (= 0 or 1), set to 0 to calculate sectoral contribtuions. Set
%       to 1 to calculate total fuel-specific cases (e.g., total coal,
%       biofuel, oil and gas). Default set to 0 (must run sectoral
%       calculation with fuel simualtions at least once before calculating
%       fuel-types).
% calc_urban, (=0 or 1), set to 0 to skip the calculation of urban-area
%       results. Set to 1 to calculate urban-area results. Default set to 1.
%       (used to save computational time)
% PM_product (string), specify the PM2.5 expsoure estimates to use. Can be
%       either 1x1km (36000 x 180000) or 10x10km (3600 x 1800) grids. Must
%       be in the form:  OR must run XXX to create downscaled exposure products. Options are 
% Year, (string), set to '2017' or '2019'. Default set to '2017'.
% grid_infile, (string), location of file that contains vectors of Lat/Lon
%       coordinates at different resolutions
% totalPM_file, (string), .mat files that contains the gridded PM2.5
%       mass from the sum of all sensitivity simulations. Created when code is
%       run with fuel_cases = 0
% pop_fileloc (string), location of gridded population data files (NOTE:
%       USER WILL NEED TO MAKE POPULATION FILES PRIOR TO RUNNING - SEE README.TXT
%       IN POPULATION_DATA FOLDER)
% mask_fileloc (string), location of files with gridded region, country,
%       and urban area masks (made using MakeGBDLandMasks_xMAPS)
% GBDdata_fileloc (string), location of the files with the relative risk
%       curves and baseline burden data for the disease burden calculation     
% DS_PM_fileloc (string), location of the files with the downscaled PM2.5
%       exposure estiamtes. (NOTE: USER WILL NEED TO MAKE POPULATION FILES 
%       PRIOR TO RUNNING - SEE README.TXT IN DOWNSCALED_PM FOLDER)
%%%%%%%%%%
% OUTPUTS:
% Four .csv files 
%   - GBD_MAPS_[outfilename]_Region_Results_Table.csv 
%       for each region - population-weighted annual average PM2.5 mass (ug/m3), fractional
%       contributions of each source, total excess number of deaths (from
%       GBD2019 CRFs and updated GEMM), fractional disease contributions to
%       total atttributable mortality (for GBD2019 CRFs and GEMM)
%   - GBD_MAPS_[outfilename]_Country_Results_Table.csv
%       for each country - same as above
%   - GBD_MAPS_[outfilename]_Urban_Results_Table.csv
%       for each urban area - population-weighted annual average PM2.5 mass
%       (ug/m3), fractional contributions of each source
%   - GBD_MAPS_[outfilename]_Final_Combined_Results_Table.csv (all combined, Data File 1 or 2 in main paper)
%
% Three .mat files (same as first three above, but in Matlab format)
%%%%%%%%%%
% DEPENDENCIES
%   A_Make_EmisSensSim_ResultsFileFn
%   B_Calc_FractionalSourceContributions_GridFn
%   C_Calc_AbsFrac_SourceContribution_RegionCountryListFn
%   D_Calc_PMmort_SourceAttFn
%   E_WriteFinalTableFn
%%%%%%%%%%%%%%%%
% USER INPUTS
Directory = 'GEOS-Chem_model_output_file_location';
final_data_location_root = 'final_output_location';
load_results_data=0;         %=0, re-calculate merged data, =1, load from outFile
skip_frac=0;                 %=0, re-calculate-fractional contributions, =1, already calculated
fuel_cases=0;                % calculate BIOFUEL & COAL & OILGAS (=1), or source sectors (=0)
calc_urban =1;               %re-calculate urban results? (1=yes, 0=no) 
PM_product = 'Downscaled_PM25_product';  %user-specified gridded PM2.5 exposure product
Year = '2017';               %year for the disease burden calculation
grid_infile = './Input_Data/Coordinates.mat';          %data file location for Lat/Lon coordinates
totalPM_file = 'gridded_PM_sum_file';                  %data file location for the sum of all sensitivity sims. (calc'd during 1st run)
pop_fileloc = './Input_Data/PopulationData/';          %location of population data files
mask_fileloc = './Input_Data/Masks/';                  %location of files with country, region, and urban gridded mask data
GBDdata_fileloc = './Input_Data/BurdenData/GBD2019/';  %location of the relevant disease burden input data
DS_PM_fileloc = './Input_Data/DownscaledPM/';          %location of the file that contains downscaled exposure estiamtes
%%%%%%%%%%%%%%%%
%PM_product =  %OPTIONS: C_DIMAQ_1km; C_DIMAQ_10km, C_DIMAQ_2019_1km; C_DIMAQ_2019_10km; C_Hammer_1km;

%%%Step 0. Create the base name for the output files, define sensitivity simulation cases for calculation
for doloop=1:1
if fuel_cases ==1; outfilename = sprintf('%s_%s','Fuel',Year);
elseif fuel_cases ==0; outfilename = sprintf('%s_%s','Sectors_wFuel',Year);
end
% specify output file names and locations
country_outfile = sprintf('%sGBD-MAPS_%s_Country_Results_Table.csv',final_data_location_root,outfilename);   %final data table location
region_outfile = sprintf('%sGBD-MAPS_%s_Region_Results_Table.csv',final_data_location_root,outfilename);
urban_outfile = sprintf('%sGBD-MAPS_%s_Urban_Results_Table.csv',final_data_location_root,outfilename);
country_data_outfile = sprintf('%sGBD-MAPS_%s_Country_Results_Data.mat',final_data_location_root,outfilename);                      
region_data_outfile = sprintf('%sGBD-MAPS_%s_Region_Results_Data.mat',final_data_location_root,outfilename);   
urban_data_outfile = sprintf('%sGBD-MAPS_%s_Urban_Results_Data.mat',final_data_location_root,outfilename);   
combined_final_outfile = sprintf('%sGBD-MAPS_%s_Final_Combined_Results_Table.csv',final_data_location_root,outfilename); 

outFile = sprintf('%sESimResults/EmissionSensitivityResults_%s.mat',final_data_location_root,outfilename);                    %location where model results will be stored (loaded with load_results_data =1)

%Define cases to calculate
if fuel_cases==0
    cases = {'base','noAGR','noENE','noENEcoal','noIND','noINDcoal','noNRTR','noROAD','noRCOR',...
        'noRCORcoal','noRCORbiofuel','noRCOC','noRCOO','noSLV','noWST','noSHP','noGFED',...
        'noGFEDagburn','noGFEDoburn','noWDUST','noAFCID','noNAT'};
elseif fuel_cases==1
    cases = {'base','noBIOFUEL','noCOAL','noOILGAS','noWDUST','noAFCID','noGFEDagburn','noGFEDoburn'}; %inlcudes extra non-fuel cases for figure purposes
end
end

%%% Step 1: Load and merge GEOS-Chem PM2.5 simulation results to a uniform 0.5x0.5 grid 
for doloop=1:1
if load_results_data==0                     % if 0, load the raw data files and merge the global and nested simulation results
    [C_merge,C_merge_aavg]=A_Make_EmisSensSim_ResultsFileFn(Directory,outFile,cases,grid_infile);
elseif load_results_data==1                 %if 0, load the results from a a previously saved version.
    load(outFile,'C_merge','C_merge_aavg')
end
end

%%% Step 2: Use model results to calulate fractional source contributions on a 0.5x0.5 grid 
[Frac_Ess_aavg]=B_Calc_FractionalSourceContributions_GridFn(cases,C_merge,totalPM_file);

%%% Step 3 - Calculate the gridded absolute contributions, then calcuate the average
%(population-weighted) absolute and fractional contributions, total
%population-weighted PM2.5 mass, and population, report by
%region/country/urban area
for doloop=1:1
if skip_frac==0
    %load the relevant gridded PM2.5 exposure level product
    if contains(Year,'2019'); pm=load(sprintf('%sGBDMAPS_DS_Results_2019.mat',DS_PM_fileloc),PM_product);
    else; pm=load(sprintf('%sGBDMAPS_DS_Results.mat',DS_PM_fileloc),PM_product);
    end
    pm = squeeze(pm.(char(sprintf('%s',PM_product))).PM25_ugm3(:,1,:,:));    %month x day x lat x lon
    pm = squeeze(nanmean(pm(:,:,:),1));
    %calculate the country and regional level products
    [AbsSource_PW_PM25_ugm3_RegionList,AbsSource_PW_PM25_ugm3_CountryList,...
        FracSource_PW_PM25_CountryList,FracSource_PW_PM25_RegionList,...
        PW_PM25_ugm3_RegionList,PW_PM25_ugm3_CountryList,Pop_CountryList,Pop_RegionList]...
     = C_Calc_AbsFrac_SourceContribution_RegionCountryListFn(...
        Frac_Ess_aavg,pm,Year,pop_fileloc,mask_fileloc,grid_infile);
    %calculate the urban-area products
    if calc_urban==1
        [AbsSource_PW_PM25_ugm3_UrbanList,FracSource_PW_PM25_UrbanList,PW_PM25_ugm3_UrbanList,Pop_UrbanList]...
        =C1_Calc_AbsFrac_SourceContribution_UrbanListFn(...
        Frac_Ess_aavg,pm,Year,pop_fileloc,mask_fileloc,grid_infile);
    else
        %if not re-calculating, make empty place holders to avoid erros in the remaining code
        AbsSource_PW_PM25_ugm3_UrbanList = NaN(200,1);
        FracSource_PW_PM25_UrbanList = NaN(200,1);
        PW_PM25_ugm3_UrbanList = NaN(200,1);
        Pop_UrbanList = cell(200,1);
    end
else
    % if skipping this re-calculation, load previously calculated data files
    load(country_data_outfile,'PW_PM25_ugm3_CountryList','Frac_Ess_aavg','FracSource_PW_PM25_CountryList','AbsSource_PW_PM25_ugm3_CountryList','Pop_CountryList');
    load(region_data_outfile,'PW_PM25_ugm3_RegionList','Frac_Ess_aavg','FracSource_PW_PM25_RegionList','AbsSource_PW_PM25_ugm3_RegionList','Pop_RegionList');
    load(urban_data_outfile,'PW_PM25_ugm3_UrbanList','Frac_Ess_aavg','FracSource_PW_PM25_UrbanList','AbsSource_PW_PM25_ugm3_UrbanList','Pop_UrbanList');
end
end

%%% Step 4: Calculate mortality associated with each source in each
%country/region - base should be equivalent to PM_product concentration (GBD2019 curve) 
%no mort data at the urban level
 [PMmort_CountryList,PMmort_RegionList]=D_Calc_PMmort_SourceAttFn(AbsSource_PW_PM25_ugm3_RegionList,AbsSource_PW_PM25_ugm3_CountryList,Year,mask_fileloc,GBDdata_fileloc);
 
%%% Step 5: Write country, region, and urban data to summary table
 [T_CountryResults, T_RegionResults,T_CityResults, T_CombinedFinal] =...
     E_WriteFinalTableFn(PMmort_CountryList,PW_PM25_ugm3_CountryList,...
     FracSource_PW_PM25_CountryList,PMmort_RegionList,...
     PW_PM25_ugm3_RegionList,FracSource_PW_PM25_RegionList,...
     PW_PM25_ugm3_UrbanList, FracSource_PW_PM25_UrbanList, country_outfile,...
     region_outfile, urban_outfile, combined_final_outfile, mask_fileloc);

%%% Step 6: Save final data and data tables as .mat and .csv files
save(country_data_outfile,'PW_PM25_ugm3_CountryList','Frac_Ess_aavg','FracSource_PW_PM25_CountryList','AbsSource_PW_PM25_ugm3_CountryList','PMmort_CountryList','Pop_CountryList');
save(region_data_outfile,'PW_PM25_ugm3_RegionList','Frac_Ess_aavg','FracSource_PW_PM25_RegionList','AbsSource_PW_PM25_ugm3_RegionList','PMmort_RegionList','Pop_RegionList');
save(urban_data_outfile,'PW_PM25_ugm3_UrbanList','Frac_Ess_aavg','FracSource_PW_PM25_UrbanList','AbsSource_PW_PM25_ugm3_UrbanList','Pop_UrbanList');

fprintf('Final Tables Saved to:\n%s\n%s\n%s\n',country_outfile,region_outfile, urban_outfile)
fprintf('Data Matrices Saved to:\n%s\n%s\n%s\n',country_data_outfile,region_data_outfile,urban_data_outfile)

clear grid interp load_results_data outFile PM_product skip_frac use_fuel cases Directory...
    cases doloop calc_urban combined_final_outfile country_data_outfile DS_PM_fileloc...
    final_data_loc_root fuel_cases GBD_fileloc grid_infile load_results_data mask_fileloc...
    outFile outfilename pm PM_product pop_fileloc region_data_outfile region_outfile skip_frac totalPM_file...
    urban_data_outfile urban_outfile Year

disp('GBD-MAPS ANALYSIS MAIN COMPLETE')
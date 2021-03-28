function [PMmort_CountryTable] = D1_GBD_MortalityCalculation_Tool_v1(countries,concentrations,HAPadjust,Year, GBDdata_fileloc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   GBD mortality calculation tool - calculate ambient PM2.5 premature
%   mortaility (combined genders) at the country level according to GBD2019 
%   concentration response functions (unadjusted for Household Air Pollutants (HAP))
%   and an updated version of the Global Exposure Mortality Model (GEMM).
%   Sources: [1] Pop: Global Health Data Exchange
%   (http://ghdx.healthdata.org/gbd-results-tool)
%   based on 'Burden Calculatr_IER 2017 model_v0.9' from Joe Spadaro
%   Mar 2020 - E. E. McDuffie
%   Updated: Nov. 8, 2020
%   Updates include new MR-BRT spline functions & GEMM & neonatal disorders
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%INPUTS:
% countries, (string vector), list of country names
% concentrations (numerical vector), list of annual average population-weighted PM2.5
%   concentrations (ug/m3) for each country
% HAPadjust, (=1 or 0), set to 1 to adjust GBD estimates for co-exposure to
%   household air pollution. Set to 0 to calculate the burden from both
%   outdoor and household air pollution.Default set to 1.
% Years, (string), exposure year, 2017 or 2019
% OUTPUTS:
% PMmort_CountryTable (table), with the following columns for each country:
%       column 1 = country names, column 2 = attributable premature mortality estiamtes (total, GBD2019), 
%       col 3 = attributable premature mortality estiamtes (total, updated GEMM),
%       col 4 = COPD (GBD2019), col 5 = COPD (GEMM), col 6 = DM (GBD2019),
%       col 7 = DM (GEMM), col 8 = LRI (GBD2019), col 9 = LRI (GEMM), col 10 = LC (GBD2019),
%       col 11 = LC (GEMM), col 12 = IHD (GBD2019), col 13 = IHD (GEMM),
%       col 14 = Stroke (GBD 2019), col 15 = Stroke (GEMM), col 16 = LBW (GBD2019)
%       col 17 = LBW (GEMM), col 18 = PTB (GBD2019), col 19 = PTB (GEMM)

%%% Step 0. Specify data paths, create country list, format output data placeholders
for doloop=1:1
% create country name list based on the input PM2.5 data
%outfile = '/misc/data6/emcduffie/ExtData/GBD/PMmortality_2017_CountryList.csv';
inputloc = GBDdata_fileloc;    %data from the 2019 GBD
lower =0;                       % set to one to calculate lower uncertainty bound (CRFs and baselines)
upper =0;                       % set to one to calculate upper uncertainty bound (CRFs and baselines)
lookup_base = 'MRBRT';
lookup_base_GEMM = 'GEMM';
%mort_base = sprintf('GBD19_%s_Baseline_Mortality',Year);
HAP_factors = 'CoExposure_Factors';
addpath(inputloc)
warning('off')
country_input_list = countries;
num_country = size(country_input_list,1);
country_PM25 = concentrations;

%if upper ==1
%    mort_base = sprintf('GBD19_%s_Baseline_Mortality_Upper',Year);
%elseif lower ==1
%    mort_base = sprintf('GBD19_%s_Baseline_Mortality_Lower',Year);
%else
    mort_base = sprintf('GBD19_%s_Baseline_Mortality',Year);
%end

%Set back to zero to use the mean CRF with baseline uncertainties
%upper=0;
%lower=0;
%
%Create output data matrices that will make up the final table
%column 1 = MRBRT (GBD2019), column 2 = GEMM
    T_country_list = cell(num_country,2);   
    T_total_PMmort = NaN(num_country,2);
    T_total_mortRate = NaN(num_country,2);
    T_AttributableFrac = NaN(num_country,2);
    T_COPDmort = NaN(num_country,2);
    T_DMmort = NaN(num_country,2);
    T_LRImort = NaN(num_country,2);
    T_LCmort = NaN(num_country,2);
    T_IHDmort = NaN(num_country,2);
    T_Strokemort = NaN(num_country,2);
    T_PTBmort = NaN(num_country,2);
    T_LBWmort = NaN(num_country,2);
end

%%% Step 1a. Load Relative Risk lookup tables for GBD2019 and updated GEMM
for doloop=1:1
% Load GBD 2019 Relative Risk lookup table data.
% Lookup table data are only function of PM2.5 (ug/m3) mass (annual average population-weighted) and age
RR_lookup.PM25_ugm3 = table2array(readtable(fullfile(inputloc,strcat(lookup_base,'_PM25.csv'))));
RR_lookup.COPD = table2array(readtable(fullfile(inputloc,strcat(lookup_base,'_COPD.csv'))));
RR_lookup.LC = table2array(readtable(fullfile(inputloc,strcat(lookup_base,'_LC.csv'))));
RR_lookup.LRI = table2array(readtable(fullfile(inputloc,strcat(lookup_base,'_LRI.csv'))));
RR_lookup.DM = table2array(readtable(fullfile(inputloc,strcat(lookup_base,'_DM.csv'))));
RR_lookup.IHD = table2array(readtable(fullfile(inputloc,strcat(lookup_base,'_IHD.csv'))));
RR_lookup.Stroke = table2array(readtable(fullfile(inputloc,strcat(lookup_base,'_Stroke.csv'))));
RR_lookup.PTB = table2array(readtable(fullfile(inputloc,strcat(lookup_base,'_PTB.csv'))));
RR_lookup.LBW = table2array(readtable(fullfile(inputloc,strcat(lookup_base,'_LBW.csv'))));

% Load GEMM (updated) Relative Risk lookup table data.
% Lookup table data are only function of PM2.5 (ug/m3) mass (annual average population-weighted) and age
RR_lookup_GEMM.PM25_ugm3 = table2array(readtable(fullfile(inputloc,strcat(lookup_base_GEMM,'_PM25.csv'))));
RR_lookup_GEMM.COPD = table2array(readtable(fullfile(inputloc,strcat(lookup_base_GEMM,'_COPD.csv'))));
RR_lookup_GEMM.LC = table2array(readtable(fullfile(inputloc,strcat(lookup_base_GEMM,'_LC.csv'))));
RR_lookup_GEMM.LRI = table2array(readtable(fullfile(inputloc,strcat(lookup_base_GEMM,'_LRI.csv'))));
RR_lookup_GEMM.DM = table2array(readtable(fullfile(inputloc,strcat(lookup_base_GEMM,'_DM.csv'))));
RR_lookup_GEMM.IHD = table2array(readtable(fullfile(inputloc,strcat(lookup_base_GEMM,'_IHD.csv'))));
RR_lookup_GEMM.Stroke = table2array(readtable(fullfile(inputloc,strcat(lookup_base_GEMM,'_Stroke.csv'))));
RR_lookup_GEMM.PTB = table2array(readtable(fullfile(inputloc,strcat(lookup_base_GEMM,'_PTB.csv'))));
RR_lookup_GEMM.LBW = table2array(readtable(fullfile(inputloc,strcat(lookup_base_GEMM,'_LBW.csv'))));
end

%%% Step 1b. Load disease- and age-specific 2017/2019 GBD19 baseline mortalities for each country
for doloop=1:1
% Load 'over25' average data for non-age specific diseases. 
% Load all age data for other variables.
% Use all ages for LRI (large fraction under 5yo)
opts=detectImportOptions(fullfile(inputloc,strcat(mort_base,'_COPD.csv')));
opts.SelectedVariableNames = {'x_Location'};
Baseline_lookup.country = table2cell(readtable(fullfile(inputloc,strcat(mort_base,'_COPD.csv')),opts)); %create baseline data country table
opts.SelectedVariableNames = {'over25'};
Baseline_lookup.COPD = table2array(readtable(fullfile(inputloc,strcat(mort_base,'_COPD.csv')), opts));  %load COPD baseline burden over 25 yo
opts=detectImportOptions(fullfile(inputloc,strcat(mort_base,'_LC.csv')));
opts.SelectedVariableNames = {'over25'};
Baseline_lookup.LC = table2array(readtable(fullfile(inputloc,strcat(mort_base,'_LC.csv')), opts));      %load LC baseline burden data over 25 yo
opts=detectImportOptions(fullfile(inputloc,strcat(mort_base,'_DM.csv')));
opts.SelectedVariableNames = {'over25'};
Baseline_lookup.DM = table2array(readtable(fullfile(inputloc,strcat(mort_base,'_DM.csv')), opts));      %load DM baseline burden data over 25 yo
opts=detectImportOptions(fullfile(inputloc,strcat(mort_base,'_LRI.csv')));
opts.SelectedVariableNames = {'Total'};
Baseline_lookup.LRI = table2array(readtable(fullfile(inputloc,strcat(mort_base,'_LRI.csv')), opts));    %load LRI baseline burden data for total (=under 5 and above 25 yo)
opts=detectImportOptions(fullfile(inputloc,strcat(mort_base,'_IHD.csv')));
opts.SelectedVariableNames = {'x25To29','x30To34','x35To39','x40To44','x45To49','x50To54','x55To59','x60To64','x65To69','x70To74','x75To79','x80To84','x85To89','x90To94','x95Plus'};
Baseline_lookup.IHD = table2array(readtable(fullfile(inputloc,strcat(mort_base,'_IHD.csv')), opts));    %load IHD baseline data from all age groups
opts=detectImportOptions(fullfile(inputloc,strcat(mort_base,'_Stroke.csv')));
opts.SelectedVariableNames = {'x25To29','x30To34','x35To39','x40To44','x45To49','x50To54','x55To59','x60To64','x65To69','x70To74','x75To79','x80To84','x85To89','x90To94','x95Plus'};
Baseline_lookup.Stroke = table2array(readtable(fullfile(inputloc,strcat(mort_base,'_Stroke.csv')), opts)); %load stroke baseline burden data from all ages
Baseline_lookup.Natural = NaN(num_country+1,1);                                                          % set natural deaths to NaN. %table2array(readtable(fullfile(inputloc,strcat(mort_base,'_Natural.csv')), opts)); not avialble for 2019 dataset
opts=detectImportOptions(fullfile(inputloc,strcat(mort_base,'_PTB.csv')));
varIn = sprintf('x%sPTBCounts',Year);
opts.SelectedVariableNames = {varIn};
Baseline_lookup.PTB = table2array(readtable(fullfile(inputloc,strcat(mort_base,'_PTB.csv')),opts));      %load baseline number of PTB incidents
opts=detectImportOptions(fullfile(inputloc,strcat(mort_base,'_LBW.csv')));
varIn = sprintf('x%sLBWCounts',Year);
opts.SelectedVariableNames = {varIn};
Baseline_lookup.LBW = table2array(readtable(fullfile(inputloc,strcat(mort_base,'_LBW.csv')),opts));     %load baseline number of LBW incidents    
end

%%%Step 2a. format country names between different input files
for doloop=1:1
%make sure the country names match those in the lookup tables. In some cases these are outdated
    Country_name_list = Baseline_lookup.country;
    Country_name_list=strrep(Country_name_list," ","_");
    Country_name_list=strrep(Country_name_list,"'","");
    Country_name_list=strrep(Country_name_list,".","");
    Country_name_list=strrep(Country_name_list,",","_");
    Country_name_list=strrep(Country_name_list,"-","_");
    Country_name_list(strcmp(Country_name_list,'Iran_(Islamic_Republic_of)'))='Iran';
    Country_name_list(strcmp(Country_name_list,'United_States_Virgin_Islands'))='Virgin_Islands_US';
    Country_name_list(strcmp(Country_name_list,'Russia'))='Russian_Federation';
    Country_name_list(strcmp(Country_name_list,'Czechia'))='Czech_Republic';
    Country_name_list(strcmp(Country_name_list,'North_Macedonia'))='Macedonia';
    Country_name_list(strcmp(Country_name_list,'Republic_of_Moldova'))='Moldova';
    Country_name_list(strcmp(Country_name_list,'Brunei_Darussalam'))='Brunei';
    Country_name_list(strcmp(Country_name_list,'Republic_of_Korea'))='South_Korea';
    Country_name_list(strcmp(Country_name_list,'United_States_of_America'))='United_States';
    Country_name_list(strcmp(Country_name_list,'Bolivia_(Plurinational_State_of)'))='Bolivia';
    Country_name_list(strcmp(Country_name_list,'Bahamas'))='The_Bahamas';
    Country_name_list(strcmp(Country_name_list,'Venezuela_(Bolivarian_Republic_of)'))='Venezuela';
    Country_name_list(strcmp(Country_name_list,'Syrian_Arab_Republic'))='Syria';
    Country_name_list(strcmp(Country_name_list,'Democratic_Peoples_Republic_of_Korea'))='North_Korea';
    Country_name_list(strcmp(Country_name_list,'Taiwan_(Province_of_China)'))='Taiwan';
    Country_name_list(strcmp(Country_name_list,'Micronesia_(Federated_States_of)'))='Federated_States_of_Micronesia';
    Country_name_list(strcmp(Country_name_list,'Lao_Peoples_Democratic_Republic'))='Laos';
    Country_name_list(strcmp(Country_name_list,'Viet_Nam'))='Vietnam';
    Country_name_list(strcmp(Country_name_list,'United_Republic_of_Tanzania'))='Tanzania';
    Country_name_list(strcmp(Country_name_list,'Eswatini'))='Swaziland';
    Country_name_list(strcmp(Country_name_list,'Cabo_Verde'))='Cape_Verde';
    Country_name_list(strcmp(Country_name_list,'Gambia'))='The_Gambia';
    
    %formulate the names from the neonatal files (just in case they are
    %different from the other diseases)
    opts=detectImportOptions(fullfile(inputloc,strcat(mort_base,'_PTB.csv')));
    opts.SelectedVariableNames = {'x_Location'};
    Country_neonatal_name_list = table2array(readtable(fullfile(inputloc,strcat(mort_base,'_PTB.csv')),opts));
    Country_neonatal_name_list=strrep(Country_neonatal_name_list," ","_");
    Country_neonatal_name_list=strrep(Country_neonatal_name_list,"'","");
    Country_neonatal_name_list=strrep(Country_neonatal_name_list,".","");
    Country_neonatal_name_list=strrep(Country_neonatal_name_list,",","_");
    Country_neonatal_name_list=strrep(Country_neonatal_name_list,"-","_");
    Country_neonatal_name_list(strcmp(Country_neonatal_name_list,'Iran_(Islamic_Republic_of)'))='Iran';
    Country_neonatal_name_list(strcmp(Country_neonatal_name_list,'United_States_Virgin_Islands'))='Virgin_Islands_US';
    Country_neonatal_name_list(strcmp(Country_neonatal_name_list,'Russia'))='Russian_Federation';
    Country_neonatal_name_list(strcmp(Country_neonatal_name_list,'Czechia'))='Czech_Republic';
    Country_neonatal_name_list(strcmp(Country_neonatal_name_list,'North_Macedonia'))='Macedonia';
    Country_neonatal_name_list(strcmp(Country_neonatal_name_list,'Republic_of_Moldova'))='Moldova';
    Country_neonatal_name_list(strcmp(Country_neonatal_name_list,'Brunei_Darussalam'))='Brunei';
    Country_neonatal_name_list(strcmp(Country_neonatal_name_list,'Republic_of_Korea'))='South_Korea';
    Country_neonatal_name_list(strcmp(Country_neonatal_name_list,'United_States_of_America'))='United_States';
    Country_neonatal_name_list(strcmp(Country_neonatal_name_list,'Bolivia_(Plurinational_State_of)'))='Bolivia';
    Country_neonatal_name_list(strcmp(Country_neonatal_name_list,'Bahamas'))='The_Bahamas';
    Country_neonatal_name_list(strcmp(Country_neonatal_name_list,'Venezuela_(Bolivarian_Republic_of)'))='Venezuela';
    Country_neonatal_name_list(strcmp(Country_neonatal_name_list,'Syrian_Arab_Republic'))='Syria';
    Country_neonatal_name_list(strcmp(Country_neonatal_name_list,'Democratic_Peoples_Republic_of_Korea'))='North_Korea';
    Country_neonatal_name_list(strcmp(Country_neonatal_name_list,'Taiwan_(Province_of_China)'))='Taiwan';
    Country_neonatal_name_list(strcmp(Country_neonatal_name_list,'Micronesia_(Federated_States_of)'))='Federated_States_of_Micronesia';
    Country_neonatal_name_list(strcmp(Country_neonatal_name_list,'Lao_Peoples_Democratic_Republic'))='Laos';
    Country_neonatal_name_list(strcmp(Country_neonatal_name_list,'Viet_Nam'))='Vietnam';
    Country_neonatal_name_list(strcmp(Country_neonatal_name_list,'United_Republic_of_Tanzania'))='Tanzania';
    Country_neonatal_name_list(strcmp(Country_neonatal_name_list,'Eswatini'))='Swaziland';
    Country_neonatal_name_list(strcmp(Country_neonatal_name_list,'Cabo_Verde'))='Cape_Verde';
    Country_neonatal_name_list(strcmp(Country_neonatal_name_list,'Gambia'))='The_Gambia';
fprintf('Calculating premature mortalities...\n')
end

%%% Step 2b. Do the calculation of attributable mortality for each country in the list
for doloop=1:1
%General Formula 
% PM2.5 Premature Excess Mortality = sum_diseases((1 - 1/RR_disease(PM2.5))*Baseline_country_disease)
% Diseases considered are: IHD, Stroke, COPD, LC, LRI, DM
% %RR's are calcualted from lookup tables (GBD2019 or GEMM)
% baselines from 2017 (or 2019) GBD
for icountry=1:num_country
    inputname = country_input_list{icountry};                                   %get country names
    fprintf('%d) %s\n', icountry,inputname)
    PM25_temp = round(country_PM25(strcmpi(country_input_list,inputname),1),1); %find the PM2.5 concentration for that country
    irow_PM = find(RR_lookup.PM25_ugm3-PM25_temp==0);                           %find the lookup table row that matches the given PM2.5 exposure level (GBD2019)
    irow2_PM = find(RR_lookup_GEMM.PM25_ugm3-PM25_temp==0);                     %find the lookup table row that matches the given PM2.5 exposure level (GEMM)
    if irow_PM ~= irow2_PM                                                      % the above sould match and throw an error if they don't. This indicates something has gone wrong in reading the burden files
        fprintf('NAMES DO NOT MATCH: %s\n',inputname)
    end
    irow_country = find(strcmpi(Country_name_list,inputname));                      %find the index of the row that corresponds to the given country
    irow_neonatal_country = find(strcmpi(Country_neonatal_name_list,inputname));    %find the index of the row that corresponds to the given country in the neonatal list
    if size(irow_country,1)==0
        fprintf('CHECK Country Name: %s\n',inputname)                               %throw and error if no country name is found. Indicates that a country in the GBD country list is not in the disease burden dataset (or the name strings don't match)
        continue;
    end
    if size(irow_PM,1)==0 
        irow_PM=1;                                                                  %if the population-weighted PM2.5 concentration does not appear in the RR lookup tables, set to zero (will occur for source categories with < 0 contributions)
        if ~(PM25_temp < 0)                                                         %only report error for NaN values (since negative values are okay). NaNs indicate an error in the previous calcualtion of pop.-weighted PM2.5 for that country
            fprintf('CHECK: %s has NaN PM2.5\n',inputname)
        end
    end
    
    %Calculate the attributable number of deaths (incidents) for the GBD2019 data
    %To calculate lower 95% CI, use 2nd column of data in the RR_lookup tables (set 1's to 2)
    %To calculate higher 95% CI, use 3nd column of data in the RR_lookup tables (set 1's to 3)
    if lower ==1
        COPD_mort_temp = (1-1/RR_lookup.COPD(irow_PM,2))* Baseline_lookup.COPD(irow_country);   
        LC_mort_temp = (1-1/RR_lookup.LC(irow_PM,2))* Baseline_lookup.LC(irow_country);
        DM_mort_temp = (1-1/RR_lookup.DM(irow_PM,2))* Baseline_lookup.DM(irow_country);
        LRI_mort_temp = (1-1/RR_lookup.LRI(irow_PM,2))* Baseline_lookup.LRI(irow_country);
        IHD_mort_temp = sum((1-1./RR_lookup.IHD(irow_PM,2:3:end)).* Baseline_lookup.IHD(irow_country,:),2);
        Stroke_mort_temp = sum((1-1./RR_lookup.Stroke(irow_PM,2:3:end)).* Baseline_lookup.Stroke(irow_country,:),2);
        Natural_mort_temp = NaN;%sum(Baseline_lookup.Natural(irow_country,:),2);
        PTB_mort_temp = (1-1/RR_lookup.PTB(irow_PM,2))*Baseline_lookup.PTB(irow_neonatal_country);
        LBW_mort_temp = (1-1/RR_lookup.LBW(irow_PM,2))*Baseline_lookup.LBW(irow_neonatal_country); 
    elseif upper ==1
        COPD_mort_temp = (1-1/RR_lookup.COPD(irow_PM,3))* Baseline_lookup.COPD(irow_country);   
        LC_mort_temp = (1-1/RR_lookup.LC(irow_PM,3))* Baseline_lookup.LC(irow_country);
        DM_mort_temp = (1-1/RR_lookup.DM(irow_PM,3))* Baseline_lookup.DM(irow_country);
        LRI_mort_temp = (1-1/RR_lookup.LRI(irow_PM,3))* Baseline_lookup.LRI(irow_country);
        IHD_mort_temp = sum((1-1./RR_lookup.IHD(irow_PM,3:3:end)).* Baseline_lookup.IHD(irow_country,:),2);
        Stroke_mort_temp = sum((1-1./RR_lookup.Stroke(irow_PM,3:3:end)).* Baseline_lookup.Stroke(irow_country,:),2);
        Natural_mort_temp = NaN;%sum(Baseline_lookup.Natural(irow_country,:),2);
        PTB_mort_temp = (1-1/RR_lookup.PTB(irow_PM,3))*Baseline_lookup.PTB(irow_neonatal_country);
        LBW_mort_temp = (1-1/RR_lookup.LBW(irow_PM,3))*Baseline_lookup.LBW(irow_neonatal_country);
    else
        COPD_mort_temp = (1-1/RR_lookup.COPD(irow_PM,1))* Baseline_lookup.COPD(irow_country);   
        LC_mort_temp = (1-1/RR_lookup.LC(irow_PM,1))* Baseline_lookup.LC(irow_country);
        DM_mort_temp = (1-1/RR_lookup.DM(irow_PM,1))* Baseline_lookup.DM(irow_country);
        LRI_mort_temp = (1-1/RR_lookup.LRI(irow_PM,1))* Baseline_lookup.LRI(irow_country);
        IHD_mort_temp = sum((1-1./RR_lookup.IHD(irow_PM,1:3:end)).* Baseline_lookup.IHD(irow_country,:),2);
        Stroke_mort_temp = sum((1-1./RR_lookup.Stroke(irow_PM,1:3:end)).* Baseline_lookup.Stroke(irow_country,:),2);
        Natural_mort_temp = sum(Baseline_lookup.Natural(irow_country,:),2);
        PTB_mort_temp = (1-1/RR_lookup.PTB(irow_PM,1))*Baseline_lookup.PTB(irow_neonatal_country);
        LBW_mort_temp = (1-1/RR_lookup.LBW(irow_PM,1))*Baseline_lookup.LBW(irow_neonatal_country);
    end
    
    %format table data
    T_country_list(icountry,1) = {inputname};                                                                           %country name
    T_total_PMmort(icountry,1) = COPD_mort_temp+LC_mort_temp+LRI_mort_temp+DM_mort_temp+IHD_mort_temp+Stroke_mort_temp; %Calculate the total number of deaths from 6 disease pairs
    T_total_mortRate(icountry,1) = T_total_PMmort(icountry,1)./1e5;                                                     %mortality rate
    T_AttributableFrac(icountry,1) = (T_total_PMmort(icountry,1)/Natural_mort_temp)*100;                                %percent = PM deaths / sum deaths
    T_COPDmort(icountry,1) = COPD_mort_temp;
    T_DMmort(icountry,1) = DM_mort_temp;
    T_LRImort(icountry,1) = LRI_mort_temp;
    T_LCmort(icountry,1) = LC_mort_temp;
    T_IHDmort(icountry,1) = IHD_mort_temp;
    T_Strokemort(icountry,1) = Stroke_mort_temp;
    T_PTBmort(icountry,1)= PTB_mort_temp;
    T_LBWmort(icountry,1) = LBW_mort_temp;
    
    %Calculate the attributable number of deaths (incidents) for the updated GEMM data
    %To calculate lower 95% CI, use 2nd column of data in the RR_lookup tables (set 1's to 2)
    %To calculate higher 95% CI, use 3nd column of data in the RR_lookup tables (set 1's to 3)
    if lower==1
        COPD_mort_temp = (1-1/RR_lookup_GEMM.COPD(irow_PM,2))* Baseline_lookup.COPD(irow_country);
        LC_mort_temp = (1-1/RR_lookup_GEMM.LC(irow_PM,2))* Baseline_lookup.LC(irow_country);
        DM_mort_temp = (1-1/RR_lookup_GEMM.DM(irow_PM,2))* Baseline_lookup.DM(irow_country);
        LRI_mort_temp = (1-1/RR_lookup_GEMM.LRI(irow_PM,2))* Baseline_lookup.LRI(irow_country);
        IHD_mort_temp = sum((1-1./RR_lookup_GEMM.IHD(irow_PM,2:3:end)).* Baseline_lookup.IHD(irow_country,:),2);
        Stroke_mort_temp = sum((1-1./RR_lookup_GEMM.Stroke(irow_PM,2:3:end)).* Baseline_lookup.Stroke(irow_country,:),2);
        Natural_mort_temp = NaN;%sum(Baseline_lookup.Natural(irow_country,:),2);
        PTB_mort_temp = (1-1/RR_lookup_GEMM.PTB(irow_PM,2))*Baseline_lookup.PTB(irow_neonatal_country);
        LBW_mort_temp = (1-1/RR_lookup_GEMM.LBW(irow_PM,2))*Baseline_lookup.LBW(irow_neonatal_country);
    elseif upper==1
        COPD_mort_temp = (1-1/RR_lookup_GEMM.COPD(irow_PM,3))* Baseline_lookup.COPD(irow_country);
        LC_mort_temp = (1-1/RR_lookup_GEMM.LC(irow_PM,3))* Baseline_lookup.LC(irow_country);
        DM_mort_temp = (1-1/RR_lookup_GEMM.DM(irow_PM,3))* Baseline_lookup.DM(irow_country);
        LRI_mort_temp = (1-1/RR_lookup_GEMM.LRI(irow_PM,3))* Baseline_lookup.LRI(irow_country);
        IHD_mort_temp = sum((1-1./RR_lookup_GEMM.IHD(irow_PM,3:3:end)).* Baseline_lookup.IHD(irow_country,:),2);
        Stroke_mort_temp = sum((1-1./RR_lookup_GEMM.Stroke(irow_PM,3:3:end)).* Baseline_lookup.Stroke(irow_country,:),2);
        Natural_mort_temp = NaN;%sum(Baseline_lookup.Natural(irow_country,:),2);
        PTB_mort_temp = (1-1/RR_lookup_GEMM.PTB(irow_PM,3))*Baseline_lookup.PTB(irow_neonatal_country);
        LBW_mort_temp = (1-1/RR_lookup_GEMM.LBW(irow_PM,3))*Baseline_lookup.LBW(irow_neonatal_country);
    else
        COPD_mort_temp = (1-1/RR_lookup_GEMM.COPD(irow_PM,1))* Baseline_lookup.COPD(irow_country);
        LC_mort_temp = (1-1/RR_lookup_GEMM.LC(irow_PM,1))* Baseline_lookup.LC(irow_country);
        DM_mort_temp = (1-1/RR_lookup_GEMM.DM(irow_PM,1))* Baseline_lookup.DM(irow_country);
        LRI_mort_temp = (1-1/RR_lookup_GEMM.LRI(irow_PM,1))* Baseline_lookup.LRI(irow_country);
        IHD_mort_temp = sum((1-1./RR_lookup_GEMM.IHD(irow_PM,1:3:end)).* Baseline_lookup.IHD(irow_country,:),2);
        Stroke_mort_temp = sum((1-1./RR_lookup_GEMM.Stroke(irow_PM,1:3:end)).* Baseline_lookup.Stroke(irow_country,:),2);
        Natural_mort_temp = sum(Baseline_lookup.Natural(irow_country,:),2);
        PTB_mort_temp = (1-1/RR_lookup_GEMM.PTB(irow_PM,1))*Baseline_lookup.PTB(irow_neonatal_country);
        LBW_mort_temp = (1-1/RR_lookup_GEMM.LBW(irow_PM,1))*Baseline_lookup.LBW(irow_neonatal_country);
    end
    %format table data
    T_country_list(icountry,2) = {inputname};
    T_total_PMmort(icountry,2) = COPD_mort_temp+LC_mort_temp+LRI_mort_temp+DM_mort_temp+IHD_mort_temp+Stroke_mort_temp; %total deaths from 6 disease pairs
    T_total_mortRate(icountry,2) = T_total_PMmort(icountry,2)./1e5;                                                     % mortality rate (#/100,000)
    T_AttributableFrac(icountry,2) = (T_total_PMmort(icountry,2)/Natural_mort_temp)*100;                                %percent  = PM deaths / sum deaths
    T_COPDmort(icountry,2) = COPD_mort_temp;
    T_DMmort(icountry,2) = DM_mort_temp;
    T_LRImort(icountry,2) = LRI_mort_temp;
    T_LCmort(icountry,2) = LC_mort_temp;
    T_IHDmort(icountry,2) = IHD_mort_temp;
    T_Strokemort(icountry,2) = Stroke_mort_temp;
    T_PTBmort(icountry,2)= PTB_mort_temp;
    T_LBWmort(icountry,2) = LBW_mort_temp;
end
end

%%% Step 3. Correct GBD2019 data to account for co-exposure to household air pollution
for doloop=1:1
% DO NOT APPLY TO GEMM
% details on the derivation of the co-exposure factors are described in the manuscript
if HAPadjust ==1
    selectyear = sprintf('x%s',Year);
    HAP_data = readtable(fullfile(inputloc,strcat(HAP_factors,'.csv')), 'ReadVariableNames',1);
    HAP_Sf = table2array(HAP_data(:,contains(HAP_data.Properties.VariableNames,selectyear)));           %assume country orders are the same as the 6 disease pairs
    HAP_countries = table2array(HAP_data(:,contains(HAP_data.Properties.VariableNames,'Location')));
    HAP_countries=strrep(HAP_countries," ","_");
    HAP_countries=strrep(HAP_countries,"'","");
    HAP_countries=strrep(HAP_countries,".","");
    HAP_countries=strrep(HAP_countries,",","_");
    HAP_countries=strrep(HAP_countries,"-","_");
    HAP_countries(strcmp(HAP_countries,'Iran_(Islamic_Republic_of)'))='Iran';
    HAP_countries(strcmp(HAP_countries,'United_States_Virgin_Islands'))='Virgin_Islands_US';
    HAP_countries(strcmp(HAP_countries,'Russia'))='Russian_Federation';
    HAP_countries(strcmp(HAP_countries,'Czechia'))='Czech_Republic';
    HAP_countries(strcmp(HAP_countries,'North_Macedonia'))='Macedonia';
    HAP_countries(strcmp(HAP_countries,'Republic_of_Moldova'))='Moldova';
    HAP_countries(strcmp(HAP_countries,'Brunei_Darussalam'))='Brunei';
    HAP_countries(strcmp(HAP_countries,'Republic_of_Korea'))='South_Korea';
    HAP_countries(strcmp(HAP_countries,'United_States_of_America'))='United_States';
    HAP_countries(strcmp(HAP_countries,'Bolivia_(Plurinational_State_of)'))='Bolivia';
    HAP_countries(strcmp(HAP_countries,'Bahamas'))='The_Bahamas';
    HAP_countries(strcmp(HAP_countries,'Venezuela_(Bolivarian_Republic_of)'))='Venezuela';
    HAP_countries(strcmp(HAP_countries,'Syrian_Arab_Republic'))='Syria';
    HAP_countries(strcmp(HAP_countries,'Democratic_Peoples_Republic_of_Korea'))='North_Korea';
    HAP_countries(strcmp(HAP_countries,'Taiwan_(Province_of_China)'))='Taiwan';
    HAP_countries(strcmp(HAP_countries,'Micronesia_(Federated_States_of)'))='Federated_States_of_Micronesia';
    HAP_countries(strcmp(HAP_countries,'Lao_Peoples_Democratic_Republic'))='Laos';
    HAP_countries(strcmp(HAP_countries,'Viet_Nam'))='Vietnam';
    HAP_countries(strcmp(HAP_countries,'United_Republic_of_Tanzania'))='Tanzania';
    HAP_countries(strcmp(HAP_countries,'Eswatini'))='Swaziland';
    HAP_countries(strcmp(HAP_countries,'Cabo_Verde'))='Cape_Verde';
    HAP_countries(strcmp(HAP_countries,'Gambia'))='The_Gambia';
    for icountry=1:num_country
        inputname = country_input_list{icountry}; %keep consistant with baseline data
        SF = HAP_Sf(strcmp(HAP_countries,inputname));
        fprintf('%d) %s: %d\n',icountry,inputname,SF)
        T_total_PMmort(icountry,1) = T_total_PMmort(icountry,1).*SF;            %remove deaths from indoor household pollution
        T_total_mortRate(icountry,1) = T_total_mortRate(icountry,1).*SF;        %remove deaths from indoor household pollution
        T_AttributableFrac(icountry,1) = T_AttributableFrac(icountry,1).*SF;    %remove deaths from indoor household pollution
        T_COPDmort(icountry,1) = T_COPDmort(icountry,1).*SF;                    %remove deaths from indoor household pollution
        T_DMmort(icountry,1) = T_DMmort(icountry,1).*SF;                        %remove deaths from indoor household pollution
        T_LRImort(icountry,1) = T_LRImort(icountry,1).*SF;                      %remove deaths from indoor household pollution
        T_LCmort(icountry,1) = T_LCmort(icountry,1).*SF;                        %remove deaths from indoor household pollution
        T_IHDmort(icountry,1) = T_IHDmort(icountry,1).*SF;                      %remove deaths from indoor household pollution
        T_Strokemort(icountry,1) = T_Strokemort(icountry,1).*SF;                %remove deaths from indoor household pollution
        T_PTBmort(icountry,1) = T_PTBmort(icountry,1).*SF;                      %remove incidents from neonatal diorders
        T_LBWmort(icountry,1) = T_LBWmort(icountry,1).*SF;                      %remove incidents from neonatal diorders
    end
end
end

%%% Step 4. Create output
for doloop=1:1
%make data table with country-level attributable mortality for GBD2019 and updated GEMM
T_data = table(T_country_list(:,1), T_total_PMmort(:,1),T_total_PMmort(:,2),T_total_mortRate(:,1),T_total_mortRate(:,2),...
    T_AttributableFrac(:,1), T_AttributableFrac(:,2),T_COPDmort(:,1),T_COPDmort(:,2),T_DMmort(:,1),T_DMmort(:,2),... 
    T_LRImort(:,1),T_LRImort(:,2), T_LCmort(:,1),T_LCmort(:,2),T_IHDmort(:,1),T_IHDmort(:,2),T_Strokemort(:,1),T_Strokemort(:,2),...
    T_PTBmort(:,1),T_PTBmort(:,2),T_LBWmort(:,1),T_LBWmort(:,2),...
    'VariableNames',{'Country_Name','Total_PMmort_MRBRT','Total_PMmort_GEMM','PMmort_rate_MRBRT','PMmort_rate_GEMM',...
    'Attributable_Fraction_MRBRT','Attributable_Fraction_GEMM','COPD_deaths_MRBRT','COPD_deaths_GEMM',...
    'DM_deaths_MRBRT','DM_deaths_GEMM','LRI_deaths_MRBRT','LRI_deaths_GEMM','LC_deaths_MRBRT','LC_GEMM',...
    'IHD_deaths_MRBRT','IHD_deaths_GEMM','Stroke_deaths_MRBRT','Stroke_deaths_GEMM','PTB_number_MRBRT','PTB_number_GEMM',...
    'LBW_number_MRBRT','LBW_number_GEMM'});
PMmort_CountryTable = T_data;
end

fprintf('COMPLETE\n')
end


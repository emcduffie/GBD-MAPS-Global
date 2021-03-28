function [T_out,T_out2,T_out3,T_out_combined]=...
    E_WriteFinalTableFn(PMmort_CountryList,PM_PW_CountryList,Frac_Country_List,...
    PMmort_RegionList,PM_PW_RegionList,Frac_Region_List,PM_PW_UrbanList,Frac_Urban_List,outfile,outfile2,outfile3, outfile4, mask_fileloc)
%E. McDuffie, last updated Nov. 8, 2020
%Write calculated data to three final tables and save to outfile locations
%
%INPUTS: 
% PMmort_CountryList (matrix), mortality data from all diseases (GBD2019 &
%       GEMM) for each country, for all sensitivity simulation cases
% PM_PW_CountryList (vector), list of annual average population-weighted
%       PM2.5 (ug/m3) for each country
% Frac_Country_List (matrix), list of population-weighted annual
%       average contributions from each source (%) in each country
% PMmort_RegionList (matrix), same as above for regions
% PM_PW_RegionList (vector), same as above for regions
% Frac_Region_List (matrix), same as above for regions
% PM_PW_UrbanList (vector), list of annual average population-weighted
%       PM2.5 (ug/m3) for each urban area
% Frac_Urban_List (matrix), list of population-weighted annual
%       average contributions from each source (%) in each urban area
% outfile (string), output location for country data
% outfile2 (string), output file location for region data
% outfile3 (string), output file location for urban area data
% outfile4 (string) output file location for combined (region, country, ubran) table
% mask_fileloc (string), file location of input mask data
%
%OUTPUTS:
% T_out (table), final data (PM2.5 annual average population-weighted
%       concentration, fractional source contributions, total attributable deaths
%       (GBD2019 & GEMM), fractional disease contributions) for each country
% T_out2 (table), same as above, for regions
% T_out3 (table), same as above, for urban areas
% T_out_combined (table), same as above, but for all region
doloop=1;

%%% Step 0. Load Country, Region, and Uran area names, format table output
%%% files, determine the names of sensitivity simulations
for doloop=1:1

load(sprintf('%sGBD_Country_Masks_0.10.mat', mask_fileloc)); %load just to get data names
load(sprintf('%sGBD_Region_Masks_0.10.mat', mask_fileloc));
load(sprintf('%sCity_Masks_0.10.mat', mask_fileloc));

%Names of sensitivity simulations
T_country_list = fieldnames(PMmort_CountryList);
numcountries = size(T_country_list,1);
Cases = fieldnames(PMmort_CountryList.(T_country_list{1}).MRBRT);
Cases = Cases(2:end);                       %don't include base case
numcases = size(Cases,1);
%format vectors to hold table data
T_total_PMmort = strings(numcountries,2);
T_COPDmort = strings(numcountries,2);
T_DMmort = strings(numcountries,2);
T_LRImort = strings(numcountries,2);
T_LCmort = strings(numcountries,2);
T_IHDmort = strings(numcountries,2);
T_Strokemort = strings(numcountries,2);
T_PTBmort = strings(numcountries,2);
T_LBWmort = strings(numcountries,2);
T_frac = strings(numcountries,numcases);
T_PW_PM = strings(numcountries,1);

%Var_shortnames
all_cases = {'AGR','ENE','ENEcoal','ENEother','IND','INDcoal','INDother','ROAD','NRTR','RCOR','RCORcoal','RCORbiofuel','RCORother','RCOC','RCOO','SLV',...
    'WST','SHP','GFEDagburn','GFEDoburn','AFCID','WDUST','BIOFUEL','COAL','OILGAS','NAT','OTHER'}';
all_longnames = {'Agriculture_Contribution_Percent','Energy_Contribution_Percent', 'Energy_Coal_Contribution_Percent',...
    'Energy_NonCoal_Contribution_Percent','Industry_Contribution_Percent','Industry_Coal_Contribution_Percent',...
    'Industry_NonCoal_Contribution_Percent','Road_Transport_Contribution_Percent','NonRoad_Transport_Contribution_Percent','Residential_Combustion_Contribution_Percent',...
    'Residential_Coal_Combustion_Contribution_Percent','Residential_Biofuel_Combustion_Contribution_Percent','Residential_Other_Combustion_Contribution_Percent',...
    'Commercial_Combustion_Contribution_Percent','Other_Combustion_Contribition_Percent', 'Solvent_Contribution_Percent','Waste_Contribution_Percent',...
    'International_Shipping_Contribution_Percent','Agricultural_Waste_Burning_Contribution_Percent', 'Other_Open_Fire_Contribution_Percent',...
    'AFCID_Dust_Contribution_Percent','Windblown_Dust_Contribution_Percent','Total_Biofuel_Contribution_Percent','Total_Coal_Contribution_Percent',...
    'Oil_and_Gas_Contribution_Percent','OtherNat_Contribution_Percent','Other_Contribution_Percent'}';
%make list of variable names based on the sensitivity simulation names in the input files
for icase=1:numcases
    idx = strcmp(all_cases,Cases{icase});
    if icase==1
        Var_longnames = all_longnames(idx);
    else
        Var_longnames = cat(1,Var_longnames,all_longnames(idx));
    end
end
end

%%% Step 1a. Format table data for each country (and calculate fraction disease contributions)
for icountry=1:numcountries
    countryname = char(T_country_list(icountry));
    % GBD2019 data
    T_total_PMmort(icountry,1) = sprintf('%1.0f',string(round(PMmort_CountryList.(countryname).MRBRT.Base.Total_PMmort/1)*1));      %round to the nearest 1
    total_mort_temp = PMmort_CountryList.(countryname).MRBRT.Base.Total_PMmort;
    T_COPDmort(icountry,1) = sprintf('%0.1f',string(PMmort_CountryList.(countryname).MRBRT.Base.PMmort_COPD*100/total_mort_temp));
    T_DMmort(icountry,1) = sprintf('%0.1f',string(PMmort_CountryList.(countryname).MRBRT.Base.PMmort_DM*100/total_mort_temp));
    T_LRImort(icountry,1) = sprintf('%0.1f',string(PMmort_CountryList.(countryname).MRBRT.Base.PMmort_LRI*100/total_mort_temp));
    T_LCmort(icountry,1) = sprintf('%0.1f',string(PMmort_CountryList.(countryname).MRBRT.Base.PMmort_LC*100/total_mort_temp));
    T_IHDmort(icountry,1) = sprintf('%0.1f',string(PMmort_CountryList.(countryname).MRBRT.Base.PMmort_IHD*100/total_mort_temp));
    T_Strokemort(icountry,1) = sprintf('%0.1f',string(PMmort_CountryList.(countryname).MRBRT.Base.PMmort_Stroke*100/total_mort_temp));
    T_PTBmort(icountry,1) = sprintf('%0.1f',string(PMmort_CountryList.(countryname).MRBRT.Base.PMmort_PTB));
    T_LBWmort(icountry,1) = sprintf('%0.1f',string(PMmort_CountryList.(countryname).MRBRT.Base.PMmort_LBW));
    T_total_PMmort(icountry,2) = sprintf('%1.0f', string(round(PMmort_CountryList.(countryname).GEMM.Base.Total_PMmort/1)*1));      %round to the nearest 500
    % Updated GEMM data
    total_mort_temp = PMmort_CountryList.(countryname).GEMM.Base.Total_PMmort;
    T_COPDmort(icountry,2) = sprintf('%0.1f',string(PMmort_CountryList.(countryname).GEMM.Base.PMmort_COPD*100/total_mort_temp));
    T_DMmort(icountry,2) = sprintf('%0.1f',string(PMmort_CountryList.(countryname).GEMM.Base.PMmort_DM*100/total_mort_temp));
    T_LRImort(icountry,2) = sprintf('%0.1f',string(PMmort_CountryList.(countryname).GEMM.Base.PMmort_LRI*100/total_mort_temp));
    T_LCmort(icountry,2) = sprintf('%0.1f',string(PMmort_CountryList.(countryname).GEMM.Base.PMmort_LC*100/total_mort_temp));
    T_IHDmort(icountry,2) = sprintf('%0.1f',string(PMmort_CountryList.(countryname).GEMM.Base.PMmort_IHD*100/total_mort_temp));
    T_Strokemort(icountry,2) = sprintf('%0.1f',string(PMmort_CountryList.(countryname).GEMM.Base.PMmort_Stroke*100/total_mort_temp));
    T_PTBmort(icountry,2) = sprintf('%0.1f',string(PMmort_CountryList.(countryname).GEMM.Base.PMmort_PTB));
    T_LBWmort(icountry,2) = sprintf('%0.1f',string(PMmort_CountryList.(countryname).GEMM.Base.PMmort_LBW));
    %sensitivity simulation, fractional source contribution results
    for icase=1:numcases
        T_frac(icountry,icase) = sprintf('%0.1f',string(Frac_Country_List.(countryname).(Cases{icase})*100));
    end
    %total annual average population-weighted PM2.5 mass (ug/3)
    T_PW_PM(icountry) = sprintf('%0.1f',string(PM_PW_CountryList.(countryname)));
end

%%% Step 1b. Save country table data
for doloop=1:1
Varnames = {'Country_Name','Population_Weighted_Annual_Average_PM25_ug_m3'};
Varnames = cat(2,Varnames,Var_longnames');
Varnames = cat(2,Varnames,'Total_Excess_Mortality_Deaths_MRBRT','Total_Excess_Mortality_Deaths_GEMM',...
    'COPD_Contribution_Percent_MRBRT','COPD_Contribution_Percent_GEMM','DM_Contribution_Percent_MRBRT','DM_Contribution_Percent_GEMM',...
    'LRI_Contribution_Percent_MRBRT','LRI_Contribution_Percent_GEMM','LC_Contribution_Percent_MRBRT','LC_Contribution_Percent_GEMM',...
    'IHD_Contribution_Percent_MRBRT','IHD_Contribution_Percent_GEMM','Stroke_Contribution_Percent_MRBRT','Stroke_Contribution_Percent_GEMM',...
    'PreTermBirths_Counts_MRBRT','PreTermBirths_Counts_GEMM','LowBirthWeight_Counts_MRBRT','LowBirthWeight_Counts_GEMM');
frac_vectors = num2cell(T_frac,1);
T_out = table(T_country_list,T_PW_PM,frac_vectors{1:end},(T_total_PMmort(:,1)),(T_total_PMmort(:,2)),...
        T_COPDmort(:,1), T_COPDmort(:,2),T_DMmort(:,1), T_DMmort(:,2),T_LRImort(:,1),T_LRImort(:,2),...
        T_LCmort(:,1),T_LCmort(:,2),T_IHDmort(:,1),T_IHDmort(:,2),T_Strokemort(:,1),T_Strokemort(:,2),...
        T_PTBmort(:,1),T_PTBmort(:,2),T_LBWmort(:,1),T_LBWmort(:,2),'VariableNames',Varnames);
writetable(T_out,outfile,'Delimiter',',')
end

%%% Step 2a. Format table data for each region (and calculate fraction disease contributions)
for doloop=1:1
%initialize matrices
T_region_list = fieldnames(PMmort_RegionList);
numregions = size(T_region_list,1);
T_total_PMmort = strings(numregions,2);
T_COPDmort = strings(numregions,2);
T_DMmort = strings(numregions,2);
T_LRImort = strings(numregions,2);
T_LCmort = strings(numregions,2);
T_IHDmort = strings(numregions,2);
T_Strokemort = strings(numregions,2);
T_PTBmort = strings(numregions,2);
T_LBWmort = strings(numregions,2);
T_frac = strings(numregions,numcases);
T_PW_PM = strings(numregions,1);

for iregion=1:numregions
    regionname = char(T_region_list(iregion));
    %GBD2019
    T_total_PMmort(iregion,1) = sprintf('%1.0f',string(round(PMmort_RegionList.(regionname).MRBRT.Base.Total_PMmort/1)*1));      %round to the nearest 1
    total_mort_temp = PMmort_RegionList.(regionname).MRBRT.Base.Total_PMmort;
    T_COPDmort(iregion,1) = sprintf('%0.1f',string(PMmort_RegionList.(regionname).MRBRT.Base.PMmort_COPD*100/total_mort_temp));
    T_DMmort(iregion,1) = sprintf('%0.1f',string(PMmort_RegionList.(regionname).MRBRT.Base.PMmort_DM*100/total_mort_temp));
    T_LRImort(iregion,1) = sprintf('%0.1f',string(PMmort_RegionList.(regionname).MRBRT.Base.PMmort_LRI*100/total_mort_temp));
    T_LCmort(iregion,1) = sprintf('%0.1f',string(PMmort_RegionList.(regionname).MRBRT.Base.PMmort_LC*100/total_mort_temp));
    T_IHDmort(iregion,1) = sprintf('%0.1f',string(PMmort_RegionList.(regionname).MRBRT.Base.PMmort_IHD*100/total_mort_temp));
    T_Strokemort(iregion,1) = sprintf('%0.1f',string(PMmort_RegionList.(regionname).MRBRT.Base.PMmort_Stroke*100/total_mort_temp));
    T_PTBmort(iregion,1) = sprintf('%0.1f',string(PMmort_RegionList.(regionname).MRBRT.Base.PMmort_PTB));
    T_LBWmort(iregion,1) = sprintf('%0.1f',string(PMmort_RegionList.(regionname).MRBRT.Base.PMmort_LBW));
    T_total_PMmort(iregion,2) = sprintf('%1.0f',string(round(PMmort_RegionList.(regionname).GEMM.Base.Total_PMmort/1)*1));
    %Updated GEMM
    total_mort_temp = PMmort_RegionList.(regionname).GEMM.Base.Total_PMmort;
    T_COPDmort(iregion,2) = sprintf('%0.1f',string(PMmort_RegionList.(regionname).GEMM.Base.PMmort_COPD*100/total_mort_temp));
    T_DMmort(iregion,2) = sprintf('%0.1f',string(PMmort_RegionList.(regionname).GEMM.Base.PMmort_DM*100/total_mort_temp));
    T_LRImort(iregion,2) = sprintf('%0.1f',string(PMmort_RegionList.(regionname).GEMM.Base.PMmort_LRI*100/total_mort_temp));
    T_LCmort(iregion,2) = sprintf('%0.1f',string(PMmort_RegionList.(regionname).GEMM.Base.PMmort_LC*100/total_mort_temp));
    T_IHDmort(iregion,2) = sprintf('%0.1f',string(PMmort_RegionList.(regionname).GEMM.Base.PMmort_IHD*100/total_mort_temp));
    T_Strokemort(iregion,2) = sprintf('%0.1f',string(PMmort_RegionList.(regionname).GEMM.Base.PMmort_Stroke*100/total_mort_temp));
    T_PTBmort(iregion,2) = sprintf('%0.1f',string(PMmort_RegionList.(regionname).GEMM.Base.PMmort_PTB));
    T_LBWmort(iregion,2) = sprintf('%0.1f',string(PMmort_RegionList.(regionname).GEMM.Base.PMmort_LBW));
    for icase=1:numcases
        T_frac(iregion,icase) = sprintf('%0.1f',string(Frac_Region_List.(regionname).(Cases{icase})*100));
    end
    T_PW_PM(iregion) = sprintf('%0.1f',string(PM_PW_RegionList.(regionname))); 
end
end

%%% Step 2b. Save region table data
for doloop=1:1
frac_vectors = num2cell(T_frac,1);
T_out2 = table(T_region_list,T_PW_PM,frac_vectors{1:end},(T_total_PMmort(:,1)),(T_total_PMmort(:,2)),...
        T_COPDmort(:,1), T_COPDmort(:,2),T_DMmort(:,1), T_DMmort(:,2),T_LRImort(:,1), T_LRImort(:,2),T_LCmort(:,1),...
        T_LCmort(:,2),T_IHDmort(:,1),T_IHDmort(:,2),T_Strokemort(:,1),T_Strokemort(:,2),...
        T_PTBmort(:,1),T_PTBmort(:,2),T_LBWmort(:,1),T_LBWmort(:,2),'VariableNames',Varnames);
writetable(T_out2,outfile2,'Delimiter',',')
end

%%% Step 3a. Format table data for each urban area (and calculate fraction disease contributions)
for doloop=1:1
%initialize matrices
T_urban_list = fieldnames(PM_PW_UrbanList);
numareas = size(T_urban_list,1);
T_total_PMmort = strings(numareas,2);
T_COPDmort = strings(numareas,2);
T_DMmort = strings(numareas,2);
T_LRImort = strings(numareas,2);
T_LCmort = strings(numareas,2);
T_IHDmort = strings(numareas,2);
T_Strokemort = strings(numareas,2);
T_PTBmort = strings(numareas,2);
T_LBWmort = strings(numareas,2);
T_frac = strings(numareas,numcases);
T_PW_PM = strings(numareas,1);

for iarea=1:numareas
    %no mortality data at the urban level
    cityname = char(T_urban_list(iarea));
    T_total_PMmort(iarea,1) = NaN;
    T_COPDmort(iarea,1) = NaN;
    T_DMmort(iarea,1) = NaN;
    T_LRImort(iarea,1) = NaN;
    T_LCmort(iarea,1) = NaN;
    T_IHDmort(iarea,1) = NaN;
    T_Strokemort(iarea,1) = NaN;
    T_PTBmort(iarea,1) = NaN;
    T_LBWmort(iarea,1) = NaN;
    T_total_PMmort(iarea,2) = NaN;
    T_COPDmort(iarea,2) = NaN;
    T_DMmort(iarea,2) = NaN;
    T_LRImort(iarea,2) = NaN;
    T_LCmort(iarea,2) = NaN;
    T_IHDmort(iarea,2) = NaN;
    T_Strokemort(iarea,2) = NaN;
    T_PTBmort(iarea,2) = NaN;
    T_LBWmort(iarea,2) = NaN;
    for icase=1:numcases
        T_frac(iarea,icase) = sprintf('%0.1f',string(Frac_Urban_List.(cityname).(Cases{icase})*100));
    end
    T_PW_PM(iarea) = sprintf('%0.1f',string(PM_PW_UrbanList.(cityname)));
end
end

%%% Step 3b. Save urban area table data
for doloop=1:1
frac_vectors = num2cell(T_frac,1);
T_out3 = table(T_urban_list,T_PW_PM,frac_vectors{1:end},(T_total_PMmort(:,1)),(T_total_PMmort(:,2)),...
        T_COPDmort(:,1), T_COPDmort(:,2),T_DMmort(:,1), T_DMmort(:,2),T_LRImort(:,1), T_LRImort(:,2),T_LCmort(:,1),...
        T_LCmort(:,2),T_IHDmort(:,1),T_IHDmort(:,2),T_Strokemort(:,1),T_Strokemort(:,2),...
        T_PTBmort(:,1),T_PTBmort(:,2),T_LBWmort(:,1),T_LBWmort(:,2),'VariableNames',Varnames);
writetable(T_out3,outfile3,'Delimiter',',')
end

%%% Step 4. Combine Country, Region, and Urban data into single table (Data Files 1 and 2 in McDuffie et al., Main Text)
for doloop=1:1
counter=1;
%get list of regions, countries, and urban areas
gbd_country_region_list = extractfield(sGBDCountries, 'Region');    %get list of countries that specify their region
gbd_country_region_list = strrep(gbd_country_region_list,' ','_');
gbd_country_region_list = strrep(gbd_country_region_list,'-','_');
gbd_country_list = extractfield(sGBDCountries, 'Name'); 
gbd_country_list = strrep(gbd_country_list,' ','_');
atlas_city_country_list = extractfield(sCityAtlas, 'Country');
atlas_city_list = extractfield(sCityAtlas,'City');
for iregion=1:numregions
    %for each region, extract the region data and find all the countries in that regions
    regionname = char(T_region_list(iregion));
    T_out_combined(counter,:) = T_out2(iregion,:);         %extract the data from the region table
    counter = counter+1;
        for icountry=1:numcountries
            %for each country, 
            country_name = char(T_country_list(icountry));
            if contains(country_name,'Virgin')
                indx = contains(gbd_country_list, 'Virgin');
            elseif contains(country_name,'Timor')
                indx = contains(gbd_country_list, 'Timor');
            elseif contains(country_name,'Ivoire')
                indx = contains(gbd_country_list, 'Ivoire');
            elseif contains(country_name,'Guinea_')
                indx = contains(gbd_country_list, 'Guinea-');
            else
                indx = strcmp(gbd_country_list, country_name);      %find the index of the given country
            end
            regiontmp = gbd_country_region_list(indx);              % find the corresponding region 
            if indx <1
                fprintf('%s, %s\n',char(regiontmp),gbd_country_list{contains(gbd_country_list, country_name(1:3))}) %print the names if the country can't be found in the list
            end
            if strcmp(regiontmp, regionname)
                T_out_combined(counter,:) = T_out(icountry, :);     %if the given country's region matches the current region, extract the data for that country and place in table
                counter = counter+1;
            else
                continue                                            %skip to next country if the given country is not in the current region
            end
            for iarea=1:numareas
                %for each urban area...
                urban_name = char(T_urban_list(iarea));
                indx = strcmp(atlas_city_list,urban_name);
                countrytmp = atlas_city_country_list(indx);
                countrytmp = strrep(countrytmp,' ','_');
                if strcmp(countrytmp,'Russia')
                    countrytmp = {'Russian_Federation'};
                elseif strcmp(countrytmp,'Congo,_DRC')
                    countrytmp = {'Democratic_Republic_of_the_Congo'};
                end
                if sum(strcmp(T_country_list, char(countrytmp{1}))) <1
                    fprintf('%s, %s\n',char(countrytmp),T_country_list{contains(T_country_list, char(countrytmp{1}(1:3)))}) %print names if the urban area name can't be found in list
                end
                if strcmp(countrytmp,country_name)
                    T_out_combined(counter,:) = T_out3(iarea,:); %if the given urban area is in the given country, then add to data table
                    counter = counter+1;
                else
                    continue                                     % skip and move on to next urban area if it is not in the given country
                end
            end 
        end
end
writetable(T_out_combined,outfile4,'Delimiter',',')             %save final data table
end

end
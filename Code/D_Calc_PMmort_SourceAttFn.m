function [PMmort_CountryList,PMmort_RegionList]=D_Calc_PMmort_SourceAttFn(AbsSource_PW_PM25_ugm3_RegionList,AbsSource_PW_PM25_ugm3_CountryList,Year,mask_fileloc,GBDdata_fileloc)
% E. McDuffie, last updated Nov. 8, 2020
% Use the GBD_Mortality calculation tool from Joe Spadaro to calculate the 
% cause-specific mortalities associated with each source of PM2.5 for each
% country/region
%
% INPUTS:
%   AbsSource_PW_PM25_ugm3_RegionList - region list of absolute average
%       population-weighte dPM2.5 source contributions (ug/m3)
%   AbsSource_PW_PM25_ugm3_CountryList - country list, same as above
%   Year - PM2.5 exposure year
%   mask_fileloc (string), file location of region and country gridded mask files
%   GBDdata_fileloc (string), file location of the GBD and GEMM relative
%   risk curves and the baseline burden data
% OUTPUTS:
%   PMmort_CountryList - country list of total attributable mortality from
%       PM2.5 according to GBD2019 and GEMM CRFs. Also report relative
%       contributions from 6 disease-pairs and number of incidents for LBW and
%       PTB
%   PMmort_RegionList - region list, same as above
% Dpendencies
%   GBD_MortalityCalculation_Tool_v1
%%%%%

    %%% Step 1. make lists of country and region names, format output files
    for doloop=1:1
    Regions = squeeze(fieldnames(AbsSource_PW_PM25_ugm3_RegionList));       % make list of region names
    Regions = unique(Regions,'stable');
    numregions = size(Regions,1);
    Countries = squeeze(fieldnames(AbsSource_PW_PM25_ugm3_CountryList));    %make list of country names
    Countries = unique(Countries, 'stable');
    numcountries = size(Countries,1);
    Cases = squeeze(fieldnames(AbsSource_PW_PM25_ugm3_CountryList.India));  %make list of sensitivity simulation names
    numcases = size(Cases,1);
    CountryList = cell(numcountries,1);
    PM25_ugm3 = NaN(numcountries,1);
    HAPadjust =1;                                                           % (=1, adjust for household air pollution co-exposure)
    clear PMmort_CountryList
    Var = {'Total_PMmort','PMmort_COPD','PMmort_DM','PMmort_LRI','PMmort_LC','PMmort_IHD','PMmort_Stroke','PMmort_PTB','PMmort_LBW'}; %disease names
    numvar=size(Var,2);
    end
    
    %%% Step 2. For each sensitivity simulation, calculate the burden associated with each disease 
    for doloop=1:1
    for icase=1:numcases
        %for each PM2.5 source...
        fprintf('Case: %s\n', Cases{icase})
        %make list of PM2.5 concenttations for each country in the country list
        for icountry=1:numcountries
            CountryList(icountry) = (Countries(icountry));
            PM25_ugm3(icountry) = AbsSource_PW_PM25_ugm3_CountryList.(Countries{icountry}).(Cases{icase});
        end
        %calculate the burden data for each country, using a MATLAB version of the tool from Joe Spadaro
        [PMmort_CountryTable] = D1_GBD_MortalityCalculation_Tool_v1(CountryList,PM25_ugm3,HAPadjust,Year,GBDdata_fileloc);
        for icountry=1:numcountries
            %for each country, reformat the results data from table to structure
            PMmort_CountryList.(Countries{icountry}).MRBRT.(Cases{icase}).Total_PMmort=table2array(PMmort_CountryTable(strcmp(table2array(PMmort_CountryTable(:,1)),Countries{icountry}),2));   %total deaths
            PMmort_CountryList.(Countries{icountry}).MRBRT.(Cases{icase}).PMmort_COPD=table2array(PMmort_CountryTable(strcmp(table2array(PMmort_CountryTable(:,1)),Countries{icountry}),8));    %COPD deaths
            PMmort_CountryList.(Countries{icountry}).MRBRT.(Cases{icase}).PMmort_DM=table2array(PMmort_CountryTable(strcmp(table2array(PMmort_CountryTable(:,1)),Countries{icountry}),10));     %DM deaths
            PMmort_CountryList.(Countries{icountry}).MRBRT.(Cases{icase}).PMmort_LRI=table2array(PMmort_CountryTable(strcmp(table2array(PMmort_CountryTable(:,1)),Countries{icountry}),12));    %LRI deaths
            PMmort_CountryList.(Countries{icountry}).MRBRT.(Cases{icase}).PMmort_LC=table2array(PMmort_CountryTable(strcmp(table2array(PMmort_CountryTable(:,1)),Countries{icountry}),14));     %LC deaths
            PMmort_CountryList.(Countries{icountry}).MRBRT.(Cases{icase}).PMmort_IHD=table2array(PMmort_CountryTable(strcmp(table2array(PMmort_CountryTable(:,1)),Countries{icountry}),16));    %IHD deaths
            PMmort_CountryList.(Countries{icountry}).MRBRT.(Cases{icase}).PMmort_Stroke=table2array(PMmort_CountryTable(strcmp(table2array(PMmort_CountryTable(:,1)),Countries{icountry}),18)); %Stroke deaths
            PMmort_CountryList.(Countries{icountry}).MRBRT.(Cases{icase}).PMmort_PTB=table2array(PMmort_CountryTable(strcmp(table2array(PMmort_CountryTable(:,1)),Countries{icountry}),20));    %PTB incidents
            PMmort_CountryList.(Countries{icountry}).MRBRT.(Cases{icase}).PMmort_LBW=table2array(PMmort_CountryTable(strcmp(table2array(PMmort_CountryTable(:,1)),Countries{icountry}),22));    %LBW incidents
            PMmort_CountryList.(Countries{icountry}).GEMM.(Cases{icase}).Total_PMmort=table2array(PMmort_CountryTable(strcmp(table2array(PMmort_CountryTable(:,1)),Countries{icountry}),3));    %total deaths
            PMmort_CountryList.(Countries{icountry}).GEMM.(Cases{icase}).PMmort_COPD=table2array(PMmort_CountryTable(strcmp(table2array(PMmort_CountryTable(:,1)),Countries{icountry}),9));     %COPD deaths
            PMmort_CountryList.(Countries{icountry}).GEMM.(Cases{icase}).PMmort_DM=table2array(PMmort_CountryTable(strcmp(table2array(PMmort_CountryTable(:,1)),Countries{icountry}),11));      %DM deaths
            PMmort_CountryList.(Countries{icountry}).GEMM.(Cases{icase}).PMmort_LRI=table2array(PMmort_CountryTable(strcmp(table2array(PMmort_CountryTable(:,1)),Countries{icountry}),13));     %LRI deaths
            PMmort_CountryList.(Countries{icountry}).GEMM.(Cases{icase}).PMmort_LC=table2array(PMmort_CountryTable(strcmp(table2array(PMmort_CountryTable(:,1)),Countries{icountry}),15));      %LC deaths
            PMmort_CountryList.(Countries{icountry}).GEMM.(Cases{icase}).PMmort_IHD=table2array(PMmort_CountryTable(strcmp(table2array(PMmort_CountryTable(:,1)),Countries{icountry}),17));     %IHD deaths
            PMmort_CountryList.(Countries{icountry}).GEMM.(Cases{icase}).PMmort_Stroke=table2array(PMmort_CountryTable(strcmp(table2array(PMmort_CountryTable(:,1)),Countries{icountry}),19));  %Stroke deaths
            PMmort_CountryList.(Countries{icountry}).GEMM.(Cases{icase}).PMmort_PTB=table2array(PMmort_CountryTable(strcmp(table2array(PMmort_CountryTable(:,1)),Countries{icountry}),21));     %PTB incidents
            PMmort_CountryList.(Countries{icountry}).GEMM.(Cases{icase}).PMmort_LBW=table2array(PMmort_CountryTable(strcmp(table2array(PMmort_CountryTable(:,1)),Countries{icountry}),23));     %LBW incidents
        end
    end
    end
    
    %%% Step 3. Calculate the regional data from a running sum of the country results (within each region)
    for doloop=1:1
    %load country names
    MaskLocName = sprintf('%sGBD_Country_Masks_0.10.mat',mask_fileloc);
    load(MaskLocName); %load just to get names
    clear PMmort_RegionList
    for iregion=1:numregions
        % for each region...
        counter=0;
        for icountry=1:numcountries
            %for each country...
            name_country = sGBDCountries(icountry).Name;
            name_country = strrep(name_country," ","_");            %replace spaces with underscores
            name_country = strrep(name_country,"-","_");            %replace dashes with underscores
            name_country = strrep(name_country,"'","");             %replace apostrophes with nothing (Cote d'Iviore)
            name_country = strrep(name_country,",","");             %replace commas with nothing (Virgin Islands, U.S.)
            name_country = char(strrep(name_country,".",""));       %replace periods with nothing (Virgin Islands, U.S.)
            country_region = sGBDCountries(icountry).Region;
            country_region = strrep(country_region," ","_");        %replace spaces with underscores
            country_region = strrep(country_region,"-","_");        %replace dashes with underscores
            if strcmp(Regions{iregion},country_region) || strcmp(Regions{iregion},'Global')
                counter=counter+1;
                %keep a running sum for each simulation case and disease variable
                for icase=1:numcases
                for ivar=1:numvar
                    if counter==1
                        PMmort_RegionList.(Regions{iregion}).MRBRT.(Cases{icase}).(Var{ivar}) = ...
                            PMmort_CountryList.(name_country).MRBRT.(Cases{icase}).(Var{ivar}) ;
                        PMmort_RegionList.(Regions{iregion}).GEMM.(Cases{icase}).(Var{ivar}) = ...
                            PMmort_CountryList.(name_country).GEMM.(Cases{icase}).(Var{ivar}) ;
                    else
                        PMmort_RegionList.(Regions{iregion}).MRBRT.(Cases{icase}).(Var{ivar}) = ...
                            PMmort_RegionList.(Regions{iregion}).MRBRT.(Cases{icase}).(Var{ivar})+PMmort_CountryList.(name_country).MRBRT.(Cases{icase}).(Var{ivar}) ;
                        PMmort_RegionList.(Regions{iregion}).GEMM.(Cases{icase}).(Var{ivar}) = ...
                            PMmort_RegionList.(Regions{iregion}).GEMM.(Cases{icase}).(Var{ivar})+PMmort_CountryList.(name_country).GEMM.(Cases{icase}).(Var{ivar}) ;
                    end
                end
                end
            end
        end
    end
    end
    

end
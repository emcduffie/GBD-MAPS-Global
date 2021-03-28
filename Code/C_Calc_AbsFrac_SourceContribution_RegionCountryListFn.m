function [AbsSource_PW_PM25_ugm3_RegionList,AbsSource_PW_PM25_ugm3_CountryList,FracSource_PW_PM25_CountryList,FracSource_PW_PM25_RegionList,PW_PM25_ugm3_RegionList,PW_PM25_ugm3_CountryList,Pop_CountryList,Pop_RegionList]...
    =C_Calc_AbsFrac_SourceContribution_RegionCountryListFn(Frac_Ess_aavg,pm_product,Year, pop_file,mask_fileloc, grid_infile)
% E. McDuffie, last updated Nov. 8, 2020
% calculates the absolute and frcational, average population-weighted
% source contributions to PM2.5 (for each world region and country)
%
% INPUTS:
%   Frac_Ess_aavg - 0.5x0.5 degree grids of fractional source contribtuions to total PM2.5
%   pm_product - specifed expsoure level product (1x1km or 10x10km resolution)
%   Year - expoure level year (2017 or 2019)
%   pop_file (string), population data file location
%   mask_fileloc (string), location of region and country gridded mask files
%   grid_infile (string), file with Lat/Lon coordinates at different resolutions
% OUTPUTS:
%   AbsSource_PW_PM25_ugm3_RegionList - region list of absolute average population-weighted source
%       contributions (ug/m3)
%   AbsSource_PW_PM25_ugm3_CountryList - country list, same as above
%   FracSource_PW_PM25_RegionList - region list of fractional average
%       population-weighted source contributions (%)
%   FracSource_PW_PM25_CountryList - country list, same as above
%   PW_PM25_ugm3_RegionList - region list of population weighted average
%       PM2.5 concentrations (ug/m3)
%   PW_PM25_ugm3_CountryList - country list, same as above
%   Pop_CountryList - region list of population
%   Pop_RegionList - country list of population
% Dependencies
%   ReGrid_xMAPS
%%%%%%


%%% Step 0. Define PM product size, load data, make list of country and
%%% region names
for doloop=1:1

%load the downscaled data and calculate annual average
%pm_product, annual average

%determine resolution of exposure product
if size(pm_product,2)==18000; high_res=1;
else; high_res=0;
end

%load the appropriate region and country masks for the given resolution
if high_res==1
    PopDataLoc = sprintf('%s%s-0.01.h5',pop_file,Year);  %0.01 resolution   
    load(sprintf('%sGBD_Country_Masks_0.10.mat',mask_fileloc));                                   %load just to get data names
    load(sprintf('%sGBD_Region_Masks_0.10.mat', mask_fileloc));                                    %load just to get data names
    grid=0;
elseif high_res==0
    PopDataLoc = sprintf('%s%s-0.10.h5',pop_file,Year);  %0.1 resolution
    load(sprintf('%sGBD_Country_Masks_0.10.mat',mask_fileloc));                                   %load just to get data names
    load(sprintf('%sGBD_Region_Masks_0.10.mat', mask_fileloc));                                    %load just to get data names
    grid =1;
end
disp('Loading population data...')
WorldPop = h5read(PopDataLoc,'/WorldPop')';
disp('Load Complete')
WorldPop(isnan(WorldPop))=0;

%grid_infile = '/data6/emcduffie/Matlab_Scripts/Coordinates.mat';   %file that defines Lat and Lon at different resolutions
load(grid_infile, 'Lat05','Lon05','Lat001','Lon001','Lon01','Lat01');  %load the appropriate lat and lon from the coordinates file

interp=0;                                               %don't interpolate when re-gridding
nest='na';                                              %placeholder variable
Regions = squeeze(extractfield(sGBDRegions,'name'));    %make list of GBD region names
Regions = unique(Regions,'stable');
numregions = size(Regions,2);
Countries = squeeze(extractfield(sGBDCountries,'Name'));%make list of GBD country names
Countries = unique(Countries,'stable');
numcountries = size(Countries,2);
Cases = squeeze(fieldnames(Frac_Ess_aavg));             % make a list of sensitivity simulation names
numcases = size(Cases,1);
end

%%% Step 1. Calculate the regional population-weighted PM2.5 and absolute
%%% and fractional contributions for each source
for doloop=1:1
for iregion=1:numregions
    %for each region...
    nametmp = strrep(Regions{iregion}," ","_"); %replace spaces with underscores in names
    nametmp = strrep(nametmp,"-","_");          %replace dashes with underscores in names
    if high_res==1
        FileName = sprintf('%sGBD_Region_Masks_%s_0.01.mat',mask_fileloc,nametmp);
        load(FileName);
        masktmp = sRegionMask.Mask; clear sRegionMask;
    elseif high_res==0
        masktmp = double(sGBDRegions(iregion).Mask);    %lonxlat mask (double)
    end
    masktmp(masktmp>0)=1;                               %ensure that all values are zero or 1
    GBD_pm = pm_product;                                % find PM for each region
    GBD_pm(masktmp==0)=0;                               % set PM data to zero outside of given region
    fractmp = Frac_Ess_aavg.Base;                       % gridded base PM2.5 at 0.05x0.05 (set to 0 in previous function)
    [fractmp2,~,~] = ReGridFn_xMAPS(grid,nest,interp,fractmp,Lat05,Lon05); %regrid to desired res (0.01 or 0.1). fractmp is 1 or 0 for base case so it is being used as mask here
    Pop_RegionList.(char(nametmp)) = nansum(nansum(WorldPop.*masktmp.*fractmp2,1));   %calculate the total population in that region (inlcudes all grid boxes, whether or not there is PM2.5 data)
    PW_PM25_ugm3_RegionList.(char(nametmp)) = nansum(nansum(WorldPop.*GBD_pm.*fractmp2,1))/Pop_RegionList.(char(nametmp));   % average population-weighted PM2.5 concentration in the region
    AbsSource_PW_PM25_ugm3_RegionList.(char(nametmp)).Base = nansum(nansum(WorldPop.*GBD_pm.*fractmp2,1))/Pop_RegionList.(char(nametmp)); %population-weighted absolute source contribution in region (= PW_PM25_ugm3 for base)
    FracSource_PW_PM25_RegionList.(char(nametmp)).Base = AbsSource_PW_PM25_ugm3_RegionList.(char(nametmp)).Base/PW_PM25_ugm3_RegionList.(char(nametmp)); %population-weighted average fractional source contribution in region (=1 in base)
    for icase=1:numcases
        % for all cases other than base...
        fractmp = Frac_Ess_aavg.(Cases{icase});                                 %gridded fractional contribution at 0.05x0.05 resolution
        [fractmp2,~,~] = ReGridFn_xMAPS(grid,nest,interp,fractmp,Lat05,Lon05);  %regrid to desired res. (0.01 or 0.1)
        fractmp2(masktmp==0)=0;                                                 %set all values outside mask to NaN
        AbsSource_PW_PM25_ugm3_RegionList.(char(nametmp)).(Cases{icase}) = nansum(nansum(WorldPop.*GBD_pm.*fractmp2,1))/Pop_RegionList.(char(nametmp)); %population-weighted absolute source contribution
        FracSource_PW_PM25_RegionList.(char(nametmp)).(Cases{icase}) = ...
            AbsSource_PW_PM25_ugm3_RegionList.(char(nametmp)).(Cases{icase})/PW_PM25_ugm3_RegionList.(char(nametmp));                                   %population-weighted fractional source contribution for region
    end
    fprintf('COMPLETE: %s\n', nametmp)
end
end

%%% Step 1b. Calculate the population-weighted PM2.5 mass and absolute and
%%% fractional source contributions for the global region. Add results to
%%% region matrix
for doloop =1:1
%calculate the global data, then add results to regional data vectors
GBD_pm = pm_product;                     % gridded annual average PM2.5
fractmp = Frac_Ess_aavg.Base;            % gridded fractional contributions at 0.05x0.05 (=1 for base)
[fractmp2,~,~] = ReGridFn_xMAPS(grid,nest,interp,fractmp,Lat05,Lon05); %regrid to desired res (0.01 or 0.1). fractmp is 1 or 0 for base case so it is being used as mask here
Pop_RegionList.(char('Global')) = nansum(nansum(WorldPop.*fractmp2,1));   %calculate the total population in that region (inlcudes all grid boxes, whether or not there is PM2.5 data)
PW_PM25_ugm3_RegionList.(char('Global')) = nansum(nansum(WorldPop.*GBD_pm.*fractmp2,1))/Pop_RegionList.(char('Global'));   %calculate the average population-weighted PM2.5 concentration in that region
AbsSource_PW_PM25_ugm3_RegionList.(char('Global')).Base = nansum(nansum(WorldPop.*GBD_pm.*fractmp2,1))/Pop_RegionList.(char('Global')); %absolute average population-weighted absolute source contribution
FracSource_PW_PM25_RegionList.(char('Global')).Base = AbsSource_PW_PM25_ugm3_RegionList.(char('Global')).Base/PW_PM25_ugm3_RegionList.(char('Global')); %populatio-weighted average fractional source contribution
for icase=1:numcases
    %for all cases other than the base...
	fractmp = Frac_Ess_aavg.(Cases{icase});                          %gridded at 0.05x0.05
	[fractmp2,~,~] = ReGridFn_xMAPS(grid,nest,interp,fractmp,Lat05,Lon05); %regrid to desired res (0.01 or 0.1)
	AbsSource_PW_PM25_ugm3_RegionList.(char('Global')).(Cases{icase}) = nansum(nansum(WorldPop.*GBD_pm.*fractmp2,1))/Pop_RegionList.(char('Global')); %absolute average population-weighted absolute source contribution
	FracSource_PW_PM25_RegionList.(char('Global')).(Cases{icase}) = ...
        AbsSource_PW_PM25_ugm3_RegionList.(char('Global')).(Cases{icase})/PW_PM25_ugm3_RegionList.(char('Global'));                                 %populatio-weighted average fractional source contribution
end
fprintf('COMPLETE: %s\n', 'Global')
end

%%% Step 2. Calculate the population-weighted PM2.5 mass and absolute and
%%% fractional source contributions for each country
for doloop=1:1
for icountry=1:numcountries
    %load gridded mask data
    nametmp = strrep(Countries{icountry}," ","_");  %replace spaces with underscores
    nametmp = strrep(nametmp,"-","_");              %replace dashes with underscores
    nametmp = strrep(nametmp,"'","");               %replace apostrophes with nothing (Cote d'Iviore)
    nametmp = strrep(nametmp,",","");               %replace commas with nothing (Virgin Islands, U.S.)
    nametmp = strrep(nametmp,".","");               %replace periods with nothing (Virgin Islands, U.S.)
    if high_res==1
        FileName = sprintf('%sGBD_Country_Masks_%s_0.01.mat',mask_fileloc,nametmp);
        load(FileName);
        masktmp = single(sCountryMask.Mask); clear sCountryMask
    elseif high_res==0
        if (strcmp(nametmp,'Marshall_Islands') || strcmp(nametmp,'Bermuda'))  %too small to be captured in 0.1x0.1 mask; set using 0.01x0.01 masks
             FileName = sprintf('%sGBD_Country_Masks_%s_0.01.mat',mask_fileloc,nametmp);
             load(FileName);
             masktmp = single(sCountryMask.Mask); clear sCountryMask;                     %0.01x0.1 mask
             [masktmp,~,~] = ReGridFn_xMAPS(grid,nest,interp,masktmp,Lat001,Lon001);      %regrid to 0.1x0.1
        else
            masktmp = double(sGBDCountries(icountry).Mask);                               %lonxlat mask (double)
        end
    end
    masktmp(masktmp>0)=1;                                                           %ensure that all values are zero or one
    GBD_pm = pm_product;                                                            % find PM for each country
    GBD_pm(masktmp==0)=0;                                                           % set PM2.5 data to zero outside of region
    fractmp = Frac_Ess_aavg.Base;                                                   %gridded at 0.05x0.05
    [fractmp2,~,~] = ReGridFn_xMAPS(grid,nest,interp,fractmp,Lat05,Lon05); %regrid to desired res (0.01 or 0.1) fractmp2 is 0 or 1 for base case
    Pop_CountryList.(char(nametmp)) = nansum(nansum(WorldPop.*masktmp.*fractmp2,1));   %calculate the total population in that region (inlcudes all grid boxes, whether or not there is PM2.5 data)
    PW_PM25_ugm3_CountryList.(char(nametmp)) = nansum(nansum(WorldPop.*GBD_pm.*fractmp2,1))/Pop_CountryList.(char(nametmp));   %calculate the average population-weighted PM2.5 concentration in that region
    AbsSource_PW_PM25_ugm3_CountryList.(char(nametmp)).Base = nansum(nansum(WorldPop.*GBD_pm.*fractmp2,1))/Pop_CountryList.(char(nametmp)); %abs pw avg PM2.5 source = country pop*frac PM2.5 soure*GBD PM/country pop
    FracSource_PW_PM25_CountryList.(char(nametmp)).Base = AbsSource_PW_PM25_ugm3_CountryList.(char(nametmp)).Base/PW_PM25_ugm3_CountryList.(char(nametmp)); %frac = abs pw average source / pw total PM2.5
    for icase=1:numcases
        %for all cases other than base...
        fractmp = Frac_Ess_aavg.(Cases{icase});
        [fractmp2,~,~] = ReGridFn_xMAPS(grid,nest,interp,fractmp,Lat05,Lon05);  %regrid to 0.1x0.1 or to 0.01x0.01
        fractmp2(masktmp==0)=NaN;
        AbsSource_PW_PM25_ugm3_CountryList.(char(nametmp)).(Cases{icase}) = nansum(nansum(WorldPop.*GBD_pm.*fractmp2,1))/Pop_CountryList.(char(nametmp)); %abs pw avg PM2.5 source = country pop*frac PM2.5 soure*GBD PM/country pop
        FracSource_PW_PM25_CountryList.(char(nametmp)).(Cases{icase}) = ...
            AbsSource_PW_PM25_ugm3_CountryList.(char(nametmp)).(Cases{icase})/PW_PM25_ugm3_CountryList.(char(nametmp)); %frac = abs pw average source / pw total PM2.5
    end
    fprintf('COMPLETE: %s\n', nametmp)
end
end

end
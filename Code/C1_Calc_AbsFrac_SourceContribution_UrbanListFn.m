function [AbsSource_PW_PM25_ugm3_UrbanList,FracSource_PW_PM25_UrbanList,PW_PM25_ugm3_UrbanList,Pop_UrbanList]...
    =C1_Calc_AbsFrac_SourceContribution_UrbanListFn(Frac_Ess_aavg,pm_product,Year, pop_file, mask_fileloc,grid_infile)
% E. McDuffie, last updated Nov. 8, 2020
% calculates the absolute and frcational, average population-weighted
% source contributions to PM2.5 (for each urban area)
%
% INPUTS:
%   Frac_Ess_aavg - 0.5x0.5 degree grids of fractional source contribtuions to total PM2.5
%   pm_product - specifed expsoure level product (1x1km or 10x10km resolution)
%   Year - expoure level year (2017 or 2019)
%   pop_file (string), population data file location
%   mask_fileloc (string), location of region and country gridded mask files
%   grid_infile (string), file with Lat/Lon coordinates at different resolutions
% OUTPUTS:
%   AbsSource_PW_PM25_ugm3_UrbanList - urban area list of absolute average population-weighted source
%       contributions (ug/m3)
%   FracSource_PW_PM25_UrbanList - urban area list of fractional average
%       population-weighted source contributions (%)
%   PW_PM25_ugm3_UrbanList - uran area list of population weighted average
%       PM2.5 concentrations (ug/m3)
%   Pop_UrbanList - urban area list of population
% Dependencies:
%   ReGrid_xMAPS
%%%%%%

%%% Step 0. Determine PM data matrix size, load masks and population, make
%%% urban area name lists
for doloop=1:1
%load population data
if size(pm_product,2)==18000; high_res=1;
else; high_res=0;
end
%load the appropriate region and country masks for the given resolution
if high_res==1
    PopDataLoc = sprintf('%s%s-0.01.h5',pop_file,Year);  %0.01 resolution   
    load(sprintf('%sCity_Masks_0.10.mat',mask_fileloc));                                   %load just to get data names
    grid=0;
elseif high_res==0
    PopDataLoc = sprintf('%s%s-0.10.h5',pop_file,Year);  %0.1 resolution
    load(sprintf('%sCity_Masks_0.10.mat',mask_fileloc));                                   
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
UrbanAreas = squeeze(extractfield(sCityAtlas,'City'));  %make list of urban area names (identified by nearest city)
UrbanAreas = unique(UrbanAreas,'stable');
numareas = size(UrbanAreas,2);
Cases = squeeze(fieldnames(Frac_Ess_aavg));             %make list of sensitivity simulation names
numcases = size(Cases,1);
end

%%% Step 1. Calculate population-weighted PM2.5 for each region, as well as
%%% the absolute and fractional source contributions
for doloop=1:1
for iarea=1:numareas
    nametmp = strrep(UrbanAreas{iarea}," ","_"); %replace spaces with underscores
    nametmp = strrep(nametmp,"-","_");           %replace dashes with underscores
    nametmp = strrep(nametmp,"'","");          %replace apostrophes with nothing (Cote d'Iviore)
    nametmp = strrep(nametmp,",","");          %replace commas with nothing (Virgin Islands, U.S.)
    nametmp = strrep(nametmp,".","");          %replace periods with nothing (Virgin Islands, U.S.)
    if high_res==1
        FileName = sprintf('%sCity_Masks_%s_0.01.mat',mask_fileloc,nametmp);
        load(FileName);
        masktmp = single(sCityMask); clear sCityMask
        fprintf('%s (%d)\n',nametmp,iarea)
    elseif high_res==0
        masktmp = double(sCityAtlas(iarea).Mask);                               %lonxlat mask (double)
    end
    masktmp(masktmp>0)=1;                                                       %ensure that all values are zero or one
    GBD_pm = pm_product;                                                        % find PM for each country
    GBD_pm(masktmp==0)=0;
    fractmp = Frac_Ess_aavg.Base;                                               %gridded at 0.05x0.05
    [fractmp2,~,~] = ReGridFn_xMAPS(grid,nest,interp,fractmp,Lat05,Lon05); %regrid to desired res (0.01 or 0.1) fractmp2 is 0 or 1 for base case
    Pop_UrbanList.(char(nametmp)) = nansum(nansum(WorldPop.*masktmp.*fractmp2,1));   %calculate the total population in that region (inlcudes all grid boxes, whether or not there is PM2.5 data)
    PW_PM25_ugm3_UrbanList.(char(nametmp)) = nansum(nansum(WorldPop.*GBD_pm.*fractmp2,1))/Pop_UrbanList.(char(nametmp));   %calculate the average population-weighted PM2.5 concentration in that region
    AbsSource_PW_PM25_ugm3_UrbanList.(char(nametmp)).Base = nansum(nansum(WorldPop.*GBD_pm.*fractmp2,1))/Pop_UrbanList.(char(nametmp)); %abs pw avg PM2.5 source = country pop*frac PM2.5 soure*GBD PM/country pop
    FracSource_PW_PM25_UrbanList.(char(nametmp)).Base = AbsSource_PW_PM25_ugm3_UrbanList.(char(nametmp)).Base/PW_PM25_ugm3_UrbanList.(char(nametmp)); %frac = abs pw average source / pw total PM2.5
    for icase=1:numcases
        %for all cases other than base
        fractmp = Frac_Ess_aavg.(Cases{icase});
        [fractmp2,~,~] = ReGridFn_xMAPS(grid,nest,interp,fractmp,Lat05,Lon05); %regrid to 0.1x0.1 or to 0.01x0.01
        fractmp2(masktmp==0)=NaN;
        AbsSource_PW_PM25_ugm3_UrbanList.(char(nametmp)).(Cases{icase}) = nansum(nansum(WorldPop.*GBD_pm.*fractmp2,1))/Pop_UrbanList.(char(nametmp)); %abs pw avg PM2.5 source = country pop*frac PM2.5 soure*GBD PM/country pop
        FracSource_PW_PM25_UrbanList.(char(nametmp)).(Cases{icase}) = ...
            AbsSource_PW_PM25_ugm3_UrbanList.(char(nametmp)).(Cases{icase})/PW_PM25_ugm3_UrbanList.(char(nametmp));                                   %frac = abs pw average source / pw total PM2.5
    end
    fprintf('COMPLETE: %s\n', nametmp)
end
end

end
function [WHO_PM25,WHO_LAT,WHO_LON,WHO_WHORegion,WHO_GBDRegion] = LoadWHOData_XMAPS(File,year)
% Function to pull extra data (for the given year) from original WHO PM files 
%- E.E. McDuffie - (11/8/18)

% Filters for points that measured PM2.5, does not include sites that
%   converted PM10 measruements to PM2.5, filters for > 75% reported data coverage
% Data are reported with UNKNOWN hygroscopic growth factors
% *** Assume PM2.5 data are at 35% RH and in volumetric units (no confirmation available) ***
%INPUTS:
%   Input data file location
%   Year of choice
%OUTPUTS: 
%   Vector of WHOPM2.5 data (annual, ug m-3), Lat, Lon, Region name, corresponding GBD
%   region for each measurement site
% INPUT DATA SOURCE: https://www.who.int/airpollution/data/cities/en/.

opts = detectImportOptions(File);
opts.SelectedVariableNames ={'Year','PM25','Latitude','Longitude','WHORegion','GBDRegion','GBDSuperRegion','PM25Conv','WebLink','PM25PercCoverage'};
Data = readtable(File,opts);
WHO_Year = table2array(Data(:,1));
WHO_PM25 = table2array(Data(:,2));
WHO_LAT = table2array(Data(:,3));
WHO_LON = table2array(Data(:,4));
WHORegion = table2cell(Data(:,5));
GBDRegion = table2cell(Data(:,6));
WHO_PM25Conv = table2array(Data(:,8));
WHO_PM25PerCov = str2double(table2array(Data(:,10)));
WHO_PM25PerCov(isnan(WHO_PM25PerCov)) =0;

%Filter data
WHO_LAT=WHO_LAT(WHO_Year==year & WHO_PM25Conv ==0 & WHO_PM25PerCov >= 0.75);
WHO_LON=WHO_LON(WHO_Year==year & WHO_PM25Conv ==0 & WHO_PM25PerCov >= 0.75);
WHO_PM25=WHO_PM25(WHO_Year==year & WHO_PM25Conv ==0 & WHO_PM25PerCov >= 0.75);
WHO_GBDRegion=GBDRegion(WHO_Year==year & WHO_PM25Conv ==0 & WHO_PM25PerCov >= 0.75);
WHO_WHORegion=WHORegion(WHO_Year==year & WHO_PM25Conv ==0 & WHO_PM25PerCov >= 0.75);

WHO_PM25 = str2double(WHO_PM25);

%    %%% *** CODE FOR CONVERSION FROM LOCAL TO STANDARD *** %%%
%     %correct data from volumetric to standard conditions using MERRA2 T and P
%     FileLoc = '/misc/data10/ctm/GEOS_0.5x0.625_NA/MERRA2/';
%     FileName = sprintf('%s%s/01/MERRA2.%s0101.A1.05x0625.NA.nc4',FileLoc,YearStr,YearStr);
%     MERRA_Lat = ncread(FileName,'lat');
%     MERRA_Lon = ncread(FileName,'lon');
%     res = [0.5,0.625];
%     for i=1:size(WHO_PM25,1)
%         if isnan(WHO_PM25(i))==0
%             iLat = find(abs(MERRA_Lat-WHO_LAT(i)) <= res(1)/2,1,'first');%find closest value of cLat for each NAPS site
%             iLon = find(abs(MERRA_Lon-WHO_LON(i)) <= res(2)/2,1,'first');
%             if isempty(iLat) || isempty(iLon)
%                if isempty(iLat) || isempty(iLon) %if point is not in nested region, find within the larger region
%                   FileLoc = '/misc/data10/ctm/GEOS_2x2.5/MERRA2/';
%                   FileName = sprintf('%s%s/01/MERRA2.%s0101.A1.2x25.nc4',FileLoc,YearStr,YearStr);
%                   MERRA_Lat = ncread(FileName,'lat');
%                   MERRA_Lon = ncread(FileName,'lon');
%                   res = [2,2.5];
%                   iLat = find(abs(MERRA_Lat - WHO_LAT(i)) <= res(1)/2,1,'first');%find closest value of cLat for each NAPS site
%                   iLon = find(abs(MERRA_Lon-WHO_LON(i)) <= res(2)/2,1,'first');
%                else
%                   FileLoc = '/misc/data10/ctm/GEOS_0.5x0.625_NA/MERRA2/';%if ok, make sure the res and fole loc are set correctly
%                   res = [0.5,0.625];
%               end
%             end
%             %Date = char(IMPROVE_Date(i));
%             %monthstr = Date(5:6);
%             %daystr = Date(7:8);
%             daysinmonth = [31,28,31,30,31,30,31,31,30,31,30,31];
%             for mn=1:12
%             if mn<=9;monthstr=sprintf('0%d',mn);else;monthstr=string(mn);end
%             for iday =1:daysinmonth(mn)
%             if iday<=9;daystr=sprintf('0%d',iday);else;daystr=string(iday);end
%             if res(1) <1
%                 FileName = sprintf('%s%s/%s/MERRA2.%s%s%s.A1.05x0625.NA.nc4',FileLoc,YearStr,monthstr,YearStr,monthstr,daystr);
%                 TmpT = ncread(FileName,'T2M');
%                 FileName = sprintf('%s%s/%s/MERRA2.%s%s%s.I3.05x0625.NA.nc4',FileLoc,YearStr,monthstr,YearStr,monthstr,daystr);
%                 TmpP = ncread(FileName,'PS');
%             else
%                 FileName = sprintf('%s%s/%s/MERRA2.%s%s%s.A1.2x25.nc4',FileLoc,YearStr,monthstr,YearStr,monthstr,daystr);
%                 TmpT = ncread(FileName,'T2M');
%                 FileName = sprintf('%s%s/%s/MERRA2.%s%s%s.I3.2x25.nc4',FileLoc,YearStr,monthstr,YearStr,monthstr,daystr);
%                 TmpP = ncread(FileName,'PS');
%             end
%             
%             MERRA2_T_d(iday) = squeeze(nanmean(TmpT(iLon,iLat,:),3)); %get daily average temp at site
%             MERRA2_P_d(iday) = squeeze(nanmean(TmpP(iLon,iLat,:),3));
%             end
%             MERRA2_T_mn(mn) = squeeze(nanmean(MERRA_T_d(:))); %get monthly average temp at site
%             MERRA2_P_mn(mn) = squeeze(nanmean(MERRA_P_d(:)));
%             end
%             MERRA2_T_yr = squeeze(nanmean(MERRA_T_mn(:))); %get yearly average temp at site
%             MERRA2_P_yr = squeeze(nanmean(MERRA_P_mn(:)));
%             if rem(i,100)==1
%                 disp(i)
%             end
%             WHO_PM25(i) = WHO_PM25(i)./((MERRA2_P_yr/101325)*(298/MERRA2_T_yr));
%         end
%     end
%    %%% END CONVERSION CODE %%%

GBDRegionDefs = {'Asia Pacific, High Income','Asia, Central','Asia, East','Asia, South','Asia, Southeast','Oceania'...
    'Australasia','Caribbean','Europe, Central','Europe, Eastern','Europe, Western','North America, High Income',...
    'Latin America, Andean',...
    'Latin America, Central','Latin America, Southern','Latin America, Tropical','North Africa / Middle East',...
    'Sub-Saharan Africa, Central','Sub-Saharan Africa, East','Sub-Saharan Africa, Southern','Sub-Saharan Africa, West','Unknown'};
WHO_GBDRegion_ColorWave = NaN((size(WHO_LAT,1)),1);
for i=1:size(GBDRegionDefs,2)
    WHO_GBDRegion_ColorWave(strcmp(WHO_GBDRegion,GBDRegionDefs{i}))=i;
end
end
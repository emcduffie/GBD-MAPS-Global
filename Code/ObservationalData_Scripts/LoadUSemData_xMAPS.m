function [USem_PM25,USem_LAT,USem_LON] = LoadUSemData_XMAPS(File,year)
% Function to extract the data (for the given year) from original US embassy PM files
% E.E. McDuffie - (9/28/19) - 
% Data are reported with UNKNOWN hygroscopic growth factors. US growth factors
%    are frequently reported at 35% RH, Europe is typically 50% RH, unknown for the others
% *** Assume PM2.5 data are at 35% RH and in volumetric units (no confirmation available) ***
%Filters for data coverage < 75%
%INPUTS:
%   Input data file location
%   desired year
%OUTPUTS
%   Vector of the annual US embassy PM2.5 concentrations (ug m-3), lat, and
%   lon at each site
%INPUT DATA SOURCE: https://www.airnow.gov/international/us-embassies-and-consulates. 

opts = detectImportOptions(File);
opts.SelectedVariableNames ={'year','mean_conc','latitude','longitude','US_Embassy_site','percent_missing'};
Data = readtable(File,opts);
USem_Year = table2array(Data(:,1));
USem_PM25 = table2array(Data(:,2));
USem_LAT = table2array(Data(:,3));
USem_LON = table2array(Data(:,4));
USem_percentmissing = table2array(Data(:,6));
USem_percentmissing(isnan(USem_percentmissing))=1;

USem_LAT=USem_LAT(USem_Year==year & USem_percentmissing < 0.25);
USem_LON=USem_LON(USem_Year==year & USem_percentmissing < 0.25);
USem_PM25=USem_PM25(USem_Year==year & USem_percentmissing < 0.25);

end
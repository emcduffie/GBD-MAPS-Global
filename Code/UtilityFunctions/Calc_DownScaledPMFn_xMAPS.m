function [C_DS] = Calc_DownScaledPMFn_xMAPS(year)
%%% Function Calc_DownscaledPMFn()
% Calculated downscaled PM2.5 gridded products using 10x10km DIMAQ and 1x1km Hammer et al., (2020) data. 
% EXAMPLES:
%[C_DIMAQ_1km]=Calc_DownScaledPMFn_xMAPS(2017);

%Description 
%%% Function to calculate the downscaled PM2.5 gridded product at ~1x1 km resolution. 
%%% This method uses the DIMAQ data (at 10x10km) as the average of the
%%% underlying 10x10 box subgrid, whose values relative to the average are taken from the Hammer et al., (2020) data.
%%% Written by E.E. McDuffie, 3/28/2019

%References:
% HAMMER data described in :
%   M. S. Hammer, A. van Donkelaar, C. Li, A. Lyapustin, A. M. Sayer, N. C. Hsu,
%   R. C. Levy, M. Garay, O. Kalashnikova, R. A. Kahn, M. Brauer, J. S. Apte, D. K. Henze,
%   L. Zhang, Q. Zhang, B. Ford, J. R. Pierce, R. V. Martin, Global Estimates and Long-Term
%   Trends of Fine Particulate Matter Concentrations (1998-2018). Environ. Sci. Technol.,(2020).
%   and avialable here: 
%       https://sites.wustl.edu/acag/datasets/surface-pm2-5/
%      
% DIMAQ data described in:
%   G. Shaddick, M. L. Thomas, H. Amini, D. Broday, A. Cohen, J. Frostad, A. Green,
%   S. Gumy, Y. Liu, R. V. Martin, A. Pruss-Ustun, D. Simpson, A. van Donkelaar, 
%   M. Brauer, Data Integration for the Assessment of Population Exposure to Ambient 
%   Air Pollution for Global Burden of Disease Assessment. Environ. Sci. Technol. 52, 9069-9078 (2018).
%   and avialable here:
%       https://www.who.int/airpollution/data/modelled-estimates/en/
%%% INPUTS:
% Year (numerical), 2017 or 2019
%%% OUTPUTS:
% C_DS.PM25_ugm3 = downscaled mass of all components and PM2.5
%   (1x1x36000x18000x1) (structure follows month x day x lat x lon x height
%   format. Since data are annual averages of surface data, there is no
%   month, day, or height information in the dataset)
%%% Dependencies
% ReGridFn_xMAPS - function to regrid datasets to different resolution
%%%%%%%

%Load raw HAMMER or DIMAQ Data (for 2017 or 2019)
% Data Locations
if year==2017
    HammerFile ='/path/to/Hammer/data/ACAG_PM25_GWR_V4GL03_201701_201712_0p01.nc' ;     
    DIMAQDataLoc = '/path/to/DIMAQ data/';                                             
    DIMAQFile = sprintf('%sDIMAQ_file_name_%s',DIMAQDataLoc,string(year));
    latpos = 2;lonpos=1;medianpos=14;
elseif year==2019
    HammerFile = '/path/to/hammer/data/ACAG_PM25_GWR_V4GL03_201901_201912_0p01.nc';     
    DIMAQDataLoc = '/path/to/DIMAQ/data/';
    DIMAQFile = sprintf('%sDIMAQ_file_name',DIMAQDataLoc);
    latpos = 3;lonpos=2;medianpos=15;
end
grid_infile = '/path/to/scripts/Coordinates.mat';                           %file that defines Lat and Lon at different resolutions
load(grid_infile,'Lat001','Lon001');
load(grid_infile,'Lat01','Lon01');


%%% Step 1. Load DIMAQ and HAMMER data, DOWN-SCALE AND BIAS-CORRECT
for doloop=1:1
    %Load DIMAQ data
	T = readtable(DIMAQFile);                                               %read .csv file (10x10 km resolution)
	SatLat=table2array(T(:,latpos));                                        %get vector of lat data
	SatLon=table2array(T(:,lonpos));                                        %get vector of lon data
	SatPM25tmp=table2array(T(:,medianpos));                                 %get vector of median PM2.5 concentration
    SatPM25 = NaN(size(Lon01,1),size(Lat01,1));                             %make palceholder matrix
    for ibox=1:size(SatLat,1)                                               % reformat the DIMAQ data so it's on the same grid as the model
        iLat = find(abs(Lat01-SatLat(ibox))<=0.05,1,'first');
        iLon = find(abs(Lon01-SatLon(ibox))<=0.05,1,'first');
        SatPM25(iLon,iLat) = SatPM25tmp(ibox);
    end
    %Load Hammer Data
	HammerPM25 = ncread(HammerFile,'PM25')';                                %read netCDF files (1x1km resolution)
	SatLat001 = ncread(HammerFile,'LAT');                                   %get vector of lat data (lon data are the same as Lon001)
	HammerFullPM25 = NaN(size(Lon001,1),size(Lat001,1));                    %make placeholder matrix
 	iLat1 = find(abs(Lat001-SatLat001(1))<=0.005,1,'first');                %find the index of the first lat value
  	iLat2 = find(abs(Lat001-SatLat001(end))<=0.005,1,'first');              %find the index of the last lat value
	HammerFullPM25(:,iLat1:iLat2) = HammerPM25(:,:);                        %insert Hammer data into full grid
	HammerFullPM25(HammerFullPM25<0)=0;                                     %there are negative values in some of the products, set to zero instead
	SatSpatialFrac = NaN(size(Lon001,1),size(Lat001,1));                    %make placeholder for a matrix where each value will represent the fraction of that grid cell relative to the surrounding 10x10 grid cells
end

%%% Step 2. Use the HAMMER data to calcualte the fractional contributions of 
%%%each point in the 1x1km grid relative to the average of the larger 10x10 grid box
for doloop=1:1
    for i=1:size(Lat01,1)-1
        if i==1
            latstart = 1;
            latrange(1)=latstart;                                               %find the index of the first box (in 10x10 km grid)
        else
            latstart = find(round(abs(Lat001-Lat01(i-1)),3)<=0.005,1,'first');  %find center of the previous 10x10km box (in corresponding 1x1km grid)
        end
        latend = find(round(abs(Lat001-Lat01(i+1)),3)<=0.005,1,'first');        % find center of next 10x10km box (in the corresponding 1x1km grid)
        distance = round(size(Lat001,1)/size(Lat01,1))./2;                      % find number of grid boxes (in 1x1km grid) along one edge of the 10x10km grid box divided by 2 (=5)
        if ~(i==1)
            latrange(1)=latstart+distance+1;                                    % set the bottom edge point in the 1x1km grid as the center lat value + the distance from the center to the edge of the larger 10x10km box
        end
        latrange(2)=latend-distance;                                            % set top edge point in the 1x1km grid as the next center lat value - the distance from the center to the edge of the 10x10km box
        fprintf('%6f\n',latrange(end))                                          % print lat value to keep track of the script status
        for j=1:size(Lon01,1)-1                                                 % for each latitude value, repeat the above process for each longitude value
            if j==1
                lonstart=1;
                lonrange(1)=lonstart;
            else
                lonstart = find(round(abs(Lon001-Lon01(j-1)),3)<=0.005,1,'first');
            end
            lonend = find(round(abs(Lon001-Lon01(j+1)),3)<=0.005,1,'first');
            distance = round(size(Lon001,1)/size(Lon01,1))./2;
            if Lon01(j) < 0                                                     % if we are in negative lon space...
                if ~(j==1)
                    lonrange(1)=lonstart+distance+1;   
                end
                    lonrange(2)=lonend-distance;
            elseif Lon01(j) ==0
                lonrange(1)=lonstart+distance+1;                                %if doing the equator, extend the box by one 
                lonrange(2)=lonend-distance+1;
            else
                lonrange(1)=lonstart+distance+1;                                %if we are in positive lon space...
                lonrange(2)=lonend-distance;
            end
            Sattmp = HammerFullPM25(lonrange(1):lonrange(end),latrange(1):latrange(end));   %extract 10x10 matrix subgrid from Hammer
            if nansum(Sattmp(:))>0                                                          %if there is Hammer data in the matrix, then calulate fractional values
                SatSpatialFrac(lonrange(1):lonrange(end),latrange(1):latrange(end)) = Sattmp./nanmean(Sattmp(:)); %fractional contributions of each grid box, relative to the mean of the 10x10 subset
            else 
                SatSpatialFrac(lonrange(1):lonrange(end),latrange(1):latrange(end)) = NaN;  %if there are entirely nans in sat data, set to nan
            end
        end
    end
	SatSpatialFrac(isnan(SatSpatialFrac))=1;    %***SET TO ONE SO THAT DIMAQ WILL NOT BE CHANGED IF NO HAMMER DATA PRESENT (important for small islands)
end

%%%Step. 3
%Downscale DIMAQ usingf fractional contributions to the average from HAMMER
%Each value in the 1x1km DIMAQ grid represent the average of the larger 10x10 km grid. 
%These values are multiplied by the fraction of each grid cell relative to
%larger average, calculated from the Hammer data
for doloop=1:1
    fprintf('PM25 Scalar Calculation Complete\n');
	%save('/misc/data6/emcduffie/Downscaling_Results2.mat', 'SatSpatialFrac','PM25Scalar','-v7.3');
	%interp=0;
    nest_region = 'na'; %placeholder (doesn't do anything)
	clear HammerPM25 SatLat001 SatPM25tmp T SatLat...
            SatPM25 PM25_ugsm3_tmp DataTmp outTmp Lat001 Lon001 SatLon
    [outTmp,~,~] = ReGridFn_xMAPS(Grid,nest_region,interp,SatPM25,Lat01,Lon01);     %re-grid DIMAQ down to 1km (without interpolating)
    outTmp=outTmp.*SatSpatialFrac;                                                   
    outTmp = single(outTmp);
    C_DS.PM25_ugm3(1,1,:,:,1) = outTmp;
end

fprintf('%s complete\n','PM2.5')
      
end
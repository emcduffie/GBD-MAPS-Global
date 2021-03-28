function [C_merge,C_merge_aavg]=A1_Merge_GlobalNestedFn(Directory,cases,speclist,interp, grid_infile)
%Merge_GlobalNestedFn
%Erin McDuffie, Oct. 16, 2020

%Merge the global GEOS-Chem simulation results with nested results
% In thefinal grid, global simulation results are replaced with the associated nested results
% OUTPUT: 
%   C_merge.{case}.{spec} at 0.5 x 0.5 degree resolution (monthly results)
%   C_merge_aavg.{case}.{spec} at 0.5 x 0.5 degree resolution (annual
%   average results)
% INPUT (examples):
%   Directory = '/data12/emcduffie/GCv12runs/rundirs/newCEDS/compute1_directories/';
%   cases = {'base'};   %can also set to 'all'
%   interp=0;   %interpolate when re-gridding (zero =no)
%   speclist = {'PM25','NIT','NH4','BC','SO4','DST1','DST2','SALA','POA','SOA','DST'};
% Dependencies:
% functions:
%   - LoadData_xMAPS
%   - ReGridFn_xMAPS
%%%%%

clear C_merge C_merge_aavg                             % clear any previous results
Loc1 = Directory;                                      % Specify directory
S=load(grid_infile, 'Lat2','Lon25');                   % load vectors of relevant grids            
Lat = S.Lat2;Lon=S.Lon25;   
S=load(grid_infile, 'Lat05_as','Lon0625_as');
Lat_as = S.Lat05_as;Lon_as=S.Lon0625_as;
S=load(grid_infile, 'Lat05_eu','Lon0625_eu');
Lat_eu = S.Lat05_eu;Lon_eu=S.Lon0625_eu;
S=load(grid_infile, 'Lat05_na','Lon0625_na');
Lat_na = S.Lat05_na;Lon_na=S.Lon0625_na;
%Lat_global = -90:0.5:90;
%Lon_global = -180:0.625:180;

%Specify simulation cases and chemical compounds to load
if strcmp(cases{1}, 'all')
    varname = {'base','noAFCID','noAGR','noBIOFUEL','noCOAL','noENE','noENEcoal','noGFED',...
    'noGFEDagburn','noGFEDoburn','noIND','noINDcoal','noNAT','noNRTR','noRCOC',...
    'noRCOO'};
else
    varname = cases;
end
numvar = size(varname,2);
if strcmp(speclist,'all')
    speclist = {'PM25','NIT','NH4','BC','SO4','DST1','DST2','SALA','POA','SOA','DST'};
end
numspec = size(speclist,2);


% for each case, make output structure (c.{case}_PM25_ugm3) with global and nested data merged onto global 0.5x0.5 grid
for ivar=1:numvar
    % for each sensitivity simulation case...
    %specify directory names
    global_directory = sprintf('%s%s_%s/base/',Loc1,'merra2_2x25_tropchem',varname{ivar});
    as_directory = sprintf('%s%s_%s/nest_as_base/',Loc1,'merra2_2x25_tropchem',varname{ivar});
    eu_directory = sprintf('%s%s_%s/nest_eu_base/',Loc1,'merra2_2x25_tropchem',varname{ivar});
    na_directory = sprintf('%s%s_%s/nest_na_base/',Loc1,'merra2_2x25_tropchem',varname{ivar});
    
    %load global and nested data structures (monthly netCDF files from GEOS-Chem)
    [C_global,~,~] = LoadData_xMAPS(1,'na',global_directory);
    [C_as,~,~] = LoadData_xMAPS(2,'as',as_directory);
    [C_eu,~,~] = LoadData_xMAPS(2,'eu',eu_directory);
    [C_na,~,~] = LoadData_xMAPS(2,'na',na_directory);
    
    for ispec=1:numspec
        %for each chemical compound...
        C_merge.(char(varname{ivar})).(sprintf('%s_ugm3',speclist{ispec})) = zeros(12,1,720,360,1);     %make output file
        C_merge_aavg.(char(varname{ivar})).(sprintf('%s_ugm3',speclist{ispec})) = zeros(720,360);       %make output file
        %for each month, 
        %1) extract the nested data, regrid to 0.5x0.5, 
        %2) extract and re-grid (0.5x0.5) global data, 
        %3) Combine all nested data into one matrix
        %4) In all places where no nested data, replace w/ corresponding global value
        for i=1:12
            intmp = squeeze(C_as.(sprintf('%s_ugm3',speclist{ispec}))(i,1,:,:,1));         %extra lonxlat data
            [outData1,~,~] = ReGridFn_xMAPS(2,'na',interp,intmp,Lat_as,Lon_as);
            outData1(isnan(outData1))=0;
            intmp = squeeze(C_eu.(sprintf('%s_ugm3',speclist{ispec}))(i,1,:,:,1));         %extra lonxlat data
            [outData2,~,~] = ReGridFn_xMAPS(2,'na',interp,intmp,Lat_eu,Lon_eu);
            outData2(isnan(outData2))=0;
            intmp = squeeze(C_na.(sprintf('%s_ugm3',speclist{ispec}))(i,1,:,:,1));         %extra lonxlat data
            [outData3,~,~] = ReGridFn_xMAPS(2,'na',interp,intmp,Lat_na,Lon_na);
            outData3(isnan(outData3))=0;
            intmp = squeeze(C_global.(sprintf('%s_ugm3',speclist{ispec}))(i,1,:,:,1));
            [outData4,~,~] = ReGridFn_xMAPS(2,'na',1,intmp,Lat,Lon);
            outData_total = outData1+outData2+outData3;             %data are zero other than nested region so can sum regions
            x = (outData_total==0);                                 %locations where nested data are zero
            outData_total(x)=outData4(x);                           %in all places where there is no nested data, replaced with global data
            C_mergetmp.(sprintf('%s_%s_ugm3',char(varname{ivar}),speclist{ispec}))(i,1,:,:,1)=outData_total; %place into final data structure
        end
        C_merge_aavg.(char(varname{ivar})).(sprintf('%s_ugm3',speclist{ispec}))(:,:)=...
            nanmean(C_mergetmp.(sprintf('%s_%s_ugm3',char(varname{ivar}),speclist{ispec}))(:,1,:,:,1),1); %calculate the annual average from monthly results
    end
    fprintf('COMPLETE...%s\n',string(varname{1}))
    
end

%Reformat the final monthy data (C_merge.(case).(compounds))
for ivar=1:numvar
    if contains(varname{ivar},'base') %save the base files for all compounds for later use
        %for ispec = 1:numspec
        %    C_merge_base.(sprintf('%s_ugm3',speclist{ispec})) = C_mergetmp.(sprintf('base_%s_ugm3',speclist{ispec}));
            C_merge.(char(varname{ivar})).(sprintf('%s_ugm3',speclist{ispec}))= C_mergetmp.(sprintf('%s_%s_ugm3',char(varname{ivar}),speclist{ispec}));
        %end
        %save('/misc/data6/emcduffie/ESimResults/GBD_MAPS/Global_2017Merge_base_ugm3.mat','C_merge_base','-v7.3');
    else
        C_merge.(char(varname{ivar})).(sprintf('%s_ugm3',speclist{ispec}))= C_mergetmp.(sprintf('%s_%s_ugm3',char(varname{ivar}),speclist{ispec}));
    end
end

fprintf('COMPLETE\n')
end
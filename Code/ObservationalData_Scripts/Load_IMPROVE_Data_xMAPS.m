function [IMPROVE_MonAvg_ugm3,IMPROVE_MonAvg_StdDev,IMPROVE_MonAvg_SiteInfo]=Load_IMPROVE_Data_xMAPS(DataFile,Year,dayfilter)
%Load IMRPOVE Data from raw data text file - raw data are reported at local
%conditions 
%E. McDuffie, Oct. 16, 2020
%INPUTS:
%   Data File - IMPROVE input data file location
%   Year - desired data year
%   dayfilter - minimum days of measurement data per month
%OUTPUTS:
%   IMPROVE_MonAvg_ugm3 is a matrix of the monthly averages for each PM2.5 component at
%   each site. 
%   IMPROVE_MonAvg_StdDev - same format as above, but StdDev
%   IMPROVE_MonAvg_SiteInfo  - lat, lon, elevation of each site

% Data are reported at local conditions - 35% RH for PM, dry for others
% Data have been filtered to remove point below LOD
%INPUT DATA SOURCE: http://views.cira.colostate.edu/fed/QueryWizard/, 
%   selecting the ‘IMPROVE Aerosol (RHR II (New Equation))’ data set
%   All sites were requested for dates between January and December 2017. 
%   Downloaded parameters include data values, site code, site location, 
%   measurement uncertainty, detection limit, and status flag. 
%   Data are saved as colon delimited .txt files in the normalized format

    YearStr = string(Year);

%Load Data
   opts = detectImportOptions(DataFile);
   opts.SelectedVariableNames = {'SiteCode','Date','ParamCode','Value','Method','Unc','MDL','StatusFlag','Latitude','Longitude','Elevation'};
   opts = setvartype(opts,{'Date','StatusFlag','Method'},'char');
   Data = readtable(DataFile,opts);
   Tmp_SiteCode = table2array(Data(:,1));
   Tmp_Date = table2cell(Data(:,2));
   Tmp_VarCode = table2array(Data(:,3));
   Tmp_Value = table2array(Data(:,4));
   Tmp_Method = table2cell(Data(:,5));
   Tmp_DL = table2array(Data(:,7));
   Tmp_SF = table2cell(Data(:,8));
   Tmp_Lat = table2array(Data(:,9));
   Tmp_Lon = table2array(Data(:,10));
   Tmp_Elev = table2array(Data(:,11));
   
   
   %Filter for the desired year and for data measured by one of the IMPROVE Module's, or a quantity calculated by CIRA (i.e. OC, EC, SS, Soil)
   IMPROVE_SiteCode =Tmp_SiteCode(contains(Tmp_Date,YearStr) & (contains(Tmp_Method, 'Module') | contains(Tmp_Method,'CIRA')));
   IMPROVE_Date =Tmp_Date(contains(Tmp_Date,YearStr) & (contains(Tmp_Method, 'Module') | contains(Tmp_Method,'CIRA')));
   IMPROVE_VarCode =Tmp_VarCode(contains(Tmp_Date,YearStr) & (contains(Tmp_Method, 'Module') | contains(Tmp_Method,'CIRA')));
   IMPROVE_Value =Tmp_Value(contains(Tmp_Date,YearStr) & (contains(Tmp_Method, 'Module') | contains(Tmp_Method,'CIRA')));
   IMPROVE_DL =Tmp_DL(contains(Tmp_Date,YearStr) & (contains(Tmp_Method, 'Module') | contains(Tmp_Method,'CIRA')));
   IMPROVE_SF = Tmp_SF(contains(Tmp_Date,YearStr) & (contains(Tmp_Method, 'Module') | contains(Tmp_Method,'CIRA')));
   IMPROVE_Lat = Tmp_Lat(contains(Tmp_Date,YearStr) & (contains(Tmp_Method, 'Module') | contains(Tmp_Method,'CIRA')));
   IMPROVE_Lon = Tmp_Lon(contains(Tmp_Date,YearStr) & (contains(Tmp_Method, 'Module') | contains(Tmp_Method,'CIRA')));
   IMPROVE_Elev = Tmp_Elev(contains(Tmp_Date,YearStr) & (contains(Tmp_Method, 'Module') | contains(Tmp_Method,'CIRA')));
   
   %filter for Values above detection limits & StatusFlag = V0, valid value
   IMPROVE_SiteCode =IMPROVE_SiteCode(IMPROVE_Value > IMPROVE_DL & strcmp(IMPROVE_SF, 'V0')); %IMPROVE_Value > IMPROVE_DL &
   IMPROVE_Date =IMPROVE_Date(IMPROVE_Value > IMPROVE_DL & strcmp(IMPROVE_SF, 'V0'));
   IMPROVE_VarCode =IMPROVE_VarCode(IMPROVE_Value > IMPROVE_DL & strcmp(IMPROVE_SF, 'V0'));
   IMPROVE_Lat =IMPROVE_Lat(IMPROVE_Value > IMPROVE_DL & strcmp(IMPROVE_SF, 'V0'));
   IMPROVE_Lon =IMPROVE_Lon(IMPROVE_Value > IMPROVE_DL & strcmp(IMPROVE_SF, 'V0'));
   IMPROVE_Elev =IMPROVE_Elev(IMPROVE_Value > IMPROVE_DL & strcmp(IMPROVE_SF, 'V0'));
   IMPROVE_Value =IMPROVE_Value(IMPROVE_Value > IMPROVE_DL & strcmp(IMPROVE_SF, 'V0')); %filter value at end to make sure logic wave is still correct   
   
   %Make monthly averages for each site (variable avg: site x month x variable )
   MonStr = {strcat(YearStr,'01'), strcat(YearStr,'02'), strcat(YearStr,'03'), strcat(YearStr,'04'), ...
       strcat(YearStr,'05'), strcat(YearStr,'06'), strcat(YearStr,'07'), strcat(YearStr,'08'), ...
       strcat(YearStr,'09'), strcat(YearStr,'10'), strcat(YearStr,'11'), strcat(YearStr,'12')};
   MonStr = cellstr(MonStr)';
   SiteStr = unique(IMPROVE_SiteCode,'stable');
   VarStr = {'MF','RCFM','NO3f','SO4f','Sf','ECf','OCf','CHLf','ALf','SIF','CAf','FEf','TIf,''NH4f'};
   numVar = size(VarStr,2);
   numsites = size(SiteStr,1);
   IMPROVE_MonAvg_SiteInfo=NaN(numsites,3); %make placeholder for SIte Info wave (site X 3 (lat, lon, elev))
   IMPROVE_MonAvg_ugm3=NaN(numsites,12,numVar);
   IMPROVE_MonAvg_StdDev=NaN(numsites,12,numVar);


   for i=1:numsites
       %Make waves of site lat, lon, elevation     
       tmp = strcmp(IMPROVE_SiteCode, SiteStr(i)); %find all points for each particular site
       IMPROVE_MonAvg_SiteInfo(i,1) = squeeze(unique(IMPROVE_Lat(tmp),'stable')); %should give single value if all lats, lons, and elevations are reported the same
       IMPROVE_MonAvg_SiteInfo(i,2) = squeeze(unique(IMPROVE_Lon(tmp),'stable'));
       IMPROVE_MonAvg_SiteInfo(i,3) = squeeze(unique(IMPROVE_Elev(tmp),'stable'));
       if size(unique(IMPROVE_Lat(tmp),'stable')) > 1 | size(unique(IMPROVE_Lon(tmp),'stable')) > 1 | size(unique(IMPROVE_Elev(tmp),'stable')) > 1
           disp('!!!!!!!!FRM LOAD ERROR!!!!!!!!!')
       end
       for j=1:12
        filtertmp = strcmp(IMPROVE_SiteCode, SiteStr(i)) & contains(IMPROVE_Date, MonStr(j));
        for k=1:numVar
            %filter for values for this site, at this month, for each
            %desired variable and set the MonAvg array accordingly
            tmp = IMPROVE_Value(filtertmp & strcmp(IMPROVE_VarCode, VarStr(k)));
            tmp(tmp<0) = NaN;
            if contains(MonStr(j),'02') 
                dayfilter2=8;
                else;dayfilter2=dayfilter;
            end
            if sum(~isnan(tmp),1) >= dayfilter2
                IMPROVE_MonAvg_ugm3(i,j,k) = nanmean(tmp);   %take average of all month data from each site for each variable
                IMPROVE_MonAvg_StdDev(i,j,k) = std(tmp,'omitnan');
            else
                IMPROVE_MonAvg_ugm3(i,j,k) = NaN;
                IMPROVE_MonAvg_StdDev(i,j,k) = NaN;
            end
        end
       end
   end
   savefile = sprintf('/misc/data6/emcduffie/ExtData/IMPROVE_Data/IMPROVE_MonthlyAvgs_%s.mat',string(Year));
   save(savefile','IMPROVE_MonAvg_ugsm3','IMPROVE_MonAvg_SiteInfo');
end
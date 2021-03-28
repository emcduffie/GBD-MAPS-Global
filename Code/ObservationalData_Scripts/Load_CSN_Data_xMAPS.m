function [CSN_MonAvg_ugsm3, CSN_MonAvg_SiteInfo] = Load_CSN_Data_xMAPS(DataFile,year,dayfilter)
%Function to load CSN data. Raw data are reported at local conditions 
%E. McDuffie - Oct 16, 2020
%INPUTS:
%   Data File - CSN input data file location
%   Year - desired data year
%   dayfilter - minimum days of measurement data per month
%OUTPUTS:
%   CSN_MonAvg_ugm3 is a matrix of the monthly average components at each site. 
%   CSN_MonAvg_SiteInfo  - lat, lon, elevation of each site
% Data have been filtered for values above DL.
%   Data were additionally filtered for the following parameter codes and methods:
%   PM2.5 (CSNIFM; gravimetric), BC (ECf_NIOSH, 88313, 88316, 88357, 88380; TOR),
%   OC (OCf_NIOSH, 88355, 88370; TOR), Chloride (CHLf, 88203; IC). 
% Data are reported at local conditions - 35% RH for PM, dry for others
%INPUT DATA SOURCE: http://views.cira.colostate.edu/fed/QueryWizard/, 
%   selecting the ‘EPA PM2.5 Speciation (CSN) – Daily’ data set. All 
%   sites were requested for dates between January and December 2017. 
%   Downloaded parameters include data values, parameter, method, and 
%   site codes, site location, measurement uncertainty, detection limit, 
%   and status flag. Data are saved as colon delimited .txt files in the normalized format

   YearStr = string(year);
%Load Data
   opts = detectImportOptions(DataFile);
   opts.SelectedVariableNames = {'SiteCode','Date','ParamCode','Value','Method','Unc','MDL','StatusFlag','Latitude','Longitude','Elevation'};%{'Var2','Var4','Var5','Var6','Var7','Var8','Var9','Var10','Var11','Var12','Var13'};%
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
   
   
   %Filter for the desired year
   CSN_SiteCode =Tmp_SiteCode(contains(Tmp_Date,YearStr));
   CSN_Date =Tmp_Date(contains(Tmp_Date,YearStr));
   CSN_VarCode =Tmp_VarCode(contains(Tmp_Date,YearStr));
   CSN_Value =Tmp_Value(contains(Tmp_Date,YearStr));
   CSN_DL =Tmp_DL(contains(Tmp_Date,YearStr));
   CSN_SF = Tmp_SF(contains(Tmp_Date,YearStr));
   CSN_Method = Tmp_Method(contains(Tmp_Date,YearStr));
   CSN_Lat = Tmp_Lat(contains(Tmp_Date,YearStr));
   CSN_Lon = Tmp_Lon(contains(Tmp_Date,YearStr));
   CSN_Elev = Tmp_Elev(contains(Tmp_Date,YearStr));
   
   %filter for Values above detection limits & StatusFlag = V0, valid value
   CSN_SiteCode =CSN_SiteCode(CSN_Value > CSN_DL & strcmp(CSN_SF, 'V0')); 
   CSN_Date =CSN_Date(CSN_Value > CSN_DL & strcmp(CSN_SF, 'V0'));
   CSN_VarCode =CSN_VarCode(CSN_Value > CSN_DL & strcmp(CSN_SF, 'V0'));
   CSN_Method =CSN_Method(CSN_Value > CSN_DL & strcmp(CSN_SF, 'V0'));
   CSN_Lat =CSN_Lat(CSN_Value > CSN_DL & strcmp(CSN_SF, 'V0'));
   CSN_Lon =CSN_Lon(CSN_Value > CSN_DL & strcmp(CSN_SF, 'V0'));
   CSN_Elev =CSN_Elev(CSN_Value > CSN_DL & strcmp(CSN_SF, 'V0'));
   CSN_Value =CSN_Value(CSN_Value > CSN_DL & strcmp(CSN_SF, 'V0'));     %filter value at end to make sure logic wave is still correct
   
   %Make monthly averages for each site (variable avg: site x month x variable )
   MonStr = {strcat(YearStr,'01'), strcat(YearStr,'02'), strcat(YearStr,'03'), strcat(YearStr,'04'), ...
       strcat(YearStr,'05'), strcat(YearStr,'06'), strcat(YearStr,'07'), strcat(YearStr,'08'), ...
       strcat(YearStr,'09'), strcat(YearStr,'10'), strcat(YearStr,'11'), strcat(YearStr,'12')};
   MonStr = cellstr(MonStr)';
   SiteStr = unique(CSN_SiteCode,'stable');
   VarStr = {'CSNIFM','NO3f','NH4f','SO4f','Sf','ECf','OCf','CHLf','ALf','SIf','CAf','FEf','TIf','MGf','NAf'};
   numVar = size(VarStr,2);
   % manually set the logic waves for the methods filters to apply for each
   % variable at each site below.s
   MethodsFilter = cell(numVar,1);
   for ivar=1:numVar
       switch true
           case strcmp(char(VarStr(ivar)), 'CSNIFM') %extract gravimetric, 24-hour and 1-hr methods, not from IMPROVE network %%can change this to inlcude beta-attneuation at a later date
               MethodsFilter{ivar,1} = contains(CSN_Method, 'Grav','IgnoreCase',true) & ~contains(CSN_Method, 'Module');% & contains(CSN_Method,'24 hour','IgnoreCase',true);
           case strcmp(char(VarStr(ivar)), 'NO3f') | strcmp(char(VarStr(ivar)), 'NH4f') | strcmp(char(VarStr(ivar)), 'SO4f') |...
                   strcmp(char(VarStr(ivar)), 'Sf') | strcmp(char(VarStr(ivar)), 'ALf') | strcmp(char(VarStr(ivar)), 'SIf') |...
                   strcmp(char(VarStr(ivar)), 'CAf') | strcmp(char(VarStr(ivar)), 'FEf') | strcmp(char(VarStr(ivar)), 'TIf')...
                   | strcmp(char(VarStr(ivar)), 'MGf') | strcmp(char(VarStr(ivar)), 'NAf') % extract 24 hour methods, not from IMPROVE network
               MethodsFilter{ivar,1} =  ~contains(CSN_Method,'Module');% & (contains(CSN_Method, '24 hour','IgnoreCase',true);
           case strcmp(char(VarStr(ivar)), 'ECf')
               CSN_VarCode(strcmp(CSN_VarCode,'ECf_NIOSH')) = {'ECf'};
               CSN_VarCode(strcmp(CSN_VarCode,'88313')) = {'ECf'};
               CSN_VarCode(strcmp(CSN_VarCode,'88316')) = {'ECf'};
               CSN_VarCode(strcmp(CSN_VarCode,'88357')) = {'ECf'};
               CSN_VarCode(strcmp(CSN_VarCode,'88380')) = {'ECf'};
               %Only 24-hour, TOR data, not from IMPROVE
               MethodsFilter{ivar,1} = (~contains(CSN_Method,'TOT') & (contains(CSN_Method,'SASS') | contains(CSN_Method, 'TOR')) & ~contains(CSN_Method,'Sunset'));%contains(CSN_Method, '24 hour','IgnoreCase',true) &
               %If I want the TOT, 24 hour data:
               %MethodsFilter{ivar,1} = (contains(CSN_Method, '24 hour','IgnoreCase',true) & contains(CSN_Method,'TOT'));
               %If I want 1-hour NIOSH data (TOT):
               %MethodsFilter{ivar,1} = (contains(CSN_Method, '1 hour','IgnoreCase',true) & contains(CSN_Method,'TOT'));
           case strcmp(char(VarStr(ivar)), 'OCf')
               CSN_VarCode(strcmp(CSN_VarCode,'88355')) = {'OCf'};
               CSN_VarCode(strcmp(CSN_VarCode,'88370')) = {'OCf'};
               CSN_VarCode(strcmp(CSN_VarCode,'OCf_NIOSH')) = {'OCf'};
               %Only 24-hour, TOR data, not from IMPROVE
               MethodsFilter{ivar,1} = (~contains(CSN_Method,'TOT') & (contains(CSN_Method,'SASS') | contains(CSN_Method, 'TOR')) & ~contains(CSN_Method,'Sunset'));%contains(CSN_Method, '24 hour','IgnoreCase',true) & 
               %If I want the TOT, 24 hour data:
               %MethodsFilter{ivar,1} = (contains(CSN_Method, '24 hour','IgnoreCase',true) & contains(CSN_Method,'TOT'));
               %If I want 1-hour NIOSH data (TOT):
               %MethodsFilter{ivar,1} = (contains(CSN_Method, '1 hour','IgnoreCase',true) & contains(CSN_Method,'TOT'));
           case strcmp(char(VarStr(ivar)), 'CHLf')
               CSN_VarCode(strcmp(CSN_VarCode,'88203')) = {'CHLf'};% VarStr(ivar) = {'88203'}; %rename data variable code to CHLf
               MethodsFilter{ivar,1} = (~contains(CSN_Method,'Module'));%contains(CSN_Method, '24 hour','IgnoreCase',true) &
       end
   end
   
   
   numsites = size(SiteStr,1);
   CSN_MonAvg_SiteInfo=NaN(numsites,3);                         %make placeholder for SIte Info wave (site X 3 (lat, lon, elev))
   CSN_MonAvg_ugsm3=NaN(numsites,12,numVar);
   for i=1:numsites
       tmp = (CSN_SiteCode == SiteStr(i));                      %find all points for each particular site
       CSN_MonAvg_SiteInfo(i,1) = unique(CSN_Lat(tmp),'stable');%should give single value if all lats, lons, and elevations are reported the same
       CSN_MonAvg_SiteInfo(i,2) = unique(CSN_Lon(tmp),'stable');
       CSN_MonAvg_SiteInfo(i,3) = unique(CSN_Elev(tmp),'stable');
       if size(unique(CSN_Lat(tmp),'stable')) > 1 | size(unique(CSN_Lon(tmp),'stable')) > 1 | size(unique(CSN_Elev(tmp),'stable')) > 1
           disp('!!!!!!!!CSN LOAD ERROR!!!!!!!!!')
       end
       for j=1:12
        filtertmp = CSN_SiteCode==SiteStr(i) & contains(CSN_Date, MonStr(j));
        for k=1:numVar
            %filter for values for this site, at this month, for each
            %desired variable and set the MonAvg array accordingly
            Methodstmp = MethodsFilter{k,1};
            tmp = CSN_Value(filtertmp & strcmp(CSN_VarCode, VarStr(k)) & Methodstmp);
            tmp(tmp<0) = NaN;
            if contains(MonStr(j),'02') 
                dayfilter2=8;
                else;dayfilter2=dayfilter;
            end
            if sum(~isnan(tmp),1) >= dayfilter2                          %only record data if >= 4 measurements per month
                CSN_MonAvg_ugsm3(i,j,k) = nanmean(tmp);         %take average of all month data from each site for each variable
            else
                CSN_MonAvg_ugsm3(i,j,k) = NaN;
            end
        end
       end
   end
   savefile = sprintf('/misc/data6/emcduffie/ExtData/CSN_Data/CSN_MonthlyAvgs_%s.mat', YearStr);
   save(savefile,'CSN_MonAvg_ugsm3','CSN_MonAvg_SiteInfo');
end
function [FRM_MonAvg_ugm3,FRM_MonAvg_SiteInfo] = Load_FRM_Data_xMAPS(DataFile,Year,dayfilter)
%Function to load FRM data. Raw data are reported at local conditions 
% E. McDuffie - Oct 16, 2020
%INPUTS:
%   Data File - FRM input data file location
%   Year - desired data year
%   dayfilter - minimum days of measurement data per month
%OUTPUTS:
%   FRM_MonAvg_ugm3 is a matrix of the monthly average PM2.5 at each site. 
%   FRM_MonAvg_SiteInfo  - lat, lon, elevation of each site

% Data are reported at local conditions - 35% RH for PM, dry for others
% Data have been filtered to remove point below LOD and for gravimetric
% measurement only
%INPUT DATA SOURCE: http://views.cira.colostate.edu/fed/QueryWizard/, 
%   selecting the ‘EPA PM2.5 Mass FRM (88101) – Daily’ data set. The 
%   EPA parameter code 88101 specifies PM2.5 measured following the FRM 
%   protocol. All sites were requested for dates between January and December 
%   2017. Downloaded parameters include data values, site codes, site location, 
%   measurement uncertainty, detection limit, and status flag. Data are saved 
%   as colon delimited .txt files in the normalized format

    YearStr = string(Year);
    
%Load Data
   opts = detectImportOptions(DataFile);
   opts.SelectedVariableNames = {'SiteCode','Date','ParamCode','Value','Unc','MDL','StatusFlag','Latitude','Longitude','Elevation','Method'};
   opts = setvartype(opts,{'Date','StatusFlag','Method'},'char');
   Data = readtable(DataFile,opts);
   Tmp_SiteCode = table2array(Data(:,1));
   Tmp_Date = table2cell(Data(:,2));
   Tmp_VarCode = table2array(Data(:,3));
   Tmp_Value = table2array(Data(:,4));
   Tmp_DL = table2array(Data(:,6));
   Tmp_SF = table2cell(Data(:,7));
   Tmp_Lat = table2cell(Data(:,8));
   Tmp_Lon = table2cell(Data(:,9));
   Tmp_Elev = table2cell(Data(:,10));
   Tmp_Method = table2cell(Data(:,11));
   
   %Filter for the desired year & Method %%%Currently filtering only for
   %Gravimetric methods, can change this to include beta attenuation
   FRM_SiteCode =Tmp_SiteCode(contains(Tmp_Date,YearStr) & (contains(Tmp_Method, 'Grav','IgnoreCase',true) & contains(Tmp_Method, '24')));
   FRM_Date =Tmp_Date(contains(Tmp_Date,YearStr) & (contains(Tmp_Method, 'Grav','IgnoreCase',true) & contains(Tmp_Method, '24')));
   FRM_VarCode =Tmp_VarCode(contains(Tmp_Date,YearStr) & (contains(Tmp_Method, 'Grav','IgnoreCase',true) & contains(Tmp_Method, '24')));
   FRM_Value =Tmp_Value(contains(Tmp_Date,YearStr) & (contains(Tmp_Method, 'Grav','IgnoreCase',true) & contains(Tmp_Method, '24')));
   FRM_DL =Tmp_DL(contains(Tmp_Date,YearStr) & (contains(Tmp_Method, 'Grav','IgnoreCase',true) & contains(Tmp_Method, '24')));
   FRM_SF = Tmp_SF(contains(Tmp_Date,YearStr) & (contains(Tmp_Method, 'Grav','IgnoreCase',true) & contains(Tmp_Method, '24')));
   FRM_Lat = Tmp_Lat(contains(Tmp_Date,YearStr) & (contains(Tmp_Method, 'Grav','IgnoreCase',true) & contains(Tmp_Method, '24')));
   FRM_Lon = Tmp_Lon(contains(Tmp_Date,YearStr) & (contains(Tmp_Method, 'Grav','IgnoreCase',true) & contains(Tmp_Method, '24')));
   FRM_Elev = Tmp_Elev(contains(Tmp_Date,YearStr) & (contains(Tmp_Method, 'Grav','IgnoreCase',true) & contains(Tmp_Method, '24')));
   
   %filter for Values above detection limits & StatusFlag = V0, valid value
   FRM_SiteCode =FRM_SiteCode(FRM_Value > FRM_DL & strcmp(FRM_SF, 'V0'));
   FRM_Date =FRM_Date(FRM_Value > FRM_DL & strcmp(FRM_SF, 'V0'));
   FRM_VarCode =FRM_VarCode(FRM_Value > FRM_DL & strcmp(FRM_SF, 'V0'));
   FRM_Lat =FRM_Lat(FRM_Value > FRM_DL & strcmp(FRM_SF, 'V0'));
   FRM_Lon =FRM_Lon(FRM_Value > FRM_DL & strcmp(FRM_SF, 'V0'));
   FRM_Elev =FRM_Elev(FRM_Value > FRM_DL & strcmp(FRM_SF, 'V0'));
   FRM_Value =FRM_Value(FRM_Value > FRM_DL & strcmp(FRM_SF, 'V0'));

   
   MetaFRM_SiteCode = unique(FRM_SiteCode,'stable');                        %pull out all the individual site codes

   %Make monthly averages for each site (variable avg: site x month x variable )
   MonStr = {strcat(YearStr,'01'), strcat(YearStr,'02'), strcat(YearStr,'03'), strcat(YearStr,'04'), ...
       strcat(YearStr,'05'), strcat(YearStr,'06'), strcat(YearStr,'07'), strcat(YearStr,'08'), ...
       strcat(YearStr,'09'), strcat(YearStr,'10'), strcat(YearStr,'11'), strcat(YearStr,'12')};
   MonStr = cellstr(MonStr)';
   SiteStr = MetaFRM_SiteCode;
   VarStr = {'MF'};
   numVar = size(VarStr,2);
   numsites = size(SiteStr,1);
   FRM_MonAvg_SiteInfo = NaN(numsites,3);                                   %make wave of sites X 3 (lat, lon, elev)
   FRM_MonAvg_ugm3=NaN(numsites,12,numVar);
   
   for i=1:numsites
       SiteStr = MetaFRM_SiteCode;
       tmp = (FRM_SiteCode == SiteStr(i));                                  %find all points for each particular site
       FRM_MonAvg_SiteInfo(i,1) = unique(cell2mat(FRM_Lat(tmp)),'stable');  %should give single value if all lats, lons, and elevations are reported the same
       FRM_MonAvg_SiteInfo(i,2) = unique(cell2mat(FRM_Lon(tmp)),'stable');
       FRM_MonAvg_SiteInfo(i,3) = unique(cell2mat(FRM_Elev(tmp)),'stable');
       if size(unique(cell2mat(FRM_Lat(tmp)),'stable')) > 1
           disp('!!!!!!!!FRM LOAD ERROR!!!!!!!!!')
       end
       for j=1:12
        filtertmp = (FRM_SiteCode == SiteStr(i)) & contains(FRM_Date, MonStr(j));
        for k=1:numVar
            %filter for values for this site, at this month, for each
            %desired variable and set the MonAvg array accordingly
            tmp = FRM_Value(filtertmp & strcmp(FRM_VarCode, VarStr(k)));
            tmp(tmp<0) = NaN;
            if contains(MonStr(j),'02') 
                dayfilter2=8;
                else;dayfilter2=dayfilter;
            end
            if sum(~isnan(tmp),1) >= dayfilter2
                FRM_MonAvg_ugm3(i,j,k) = nanmean(tmp);                     %take average of all month data from each site for each variable
            else
                FRM_MonAvg_ugm3(i,j,k) = NaN;
            end
        end
       end
   end
   savefile = sprintf('/misc/data6/emcduffie/ExtData/FRM_Data/FRM_MonthlyAvgs_%s.mat',YearStr);
   save(savefile,'FRM_MonAvg_ugsm3','FRM_MonAvg_SiteInfo');
end
function [SPARTAN_MonAvg_ugsm3, SPARTAN_MonAvg_SiteInfo] = Load_SPARTAN_Data_xMAPS(DataFileLoc,year)
% loads data from SPARTAN files, downloaded from the Website
% data are reported at ambinet conditions
% E. McDuffie - Oct. 16, 2020
%INPUTS:
%   Data File - SPARTAN input data file location
%   Year - desired data year
%OUTPUTS:
%   SPARTAN_MonAvg_ugm3 is a matrix of the monthly average components at each site. 
%   SPARTAN_MonAvg_SiteInfo  - lat, lon, elevation of each site
% PM2.5 is analyzed at 35% RH
%INPUT DATA SOURCE: https://www.spartan- network.org/data

%Additional notes:
%   Data are filtered to include the following parameter codes: 
%   PM2.5 (28101), SO42- (28401), BC (28201), OM (28306), 
%   trace metals (28801, 28804, 28805, 28902, 28904, 28908), Dust (28305), and Sea Salt (28304). 
%   Due to known loss of ammonium and nitrate on Teflon filters, 
%   these data have been excluded from this analysis. 
    
    YearStr = string(year);
    listing = dir(DataFileLoc);
    numFiles = size(listing,1);
    FileNames = cell(numFiles,1);
   for i=1:numFiles
       FileNames(i,1) = cellstr(listing(i).name);   %makes cell array of the file names in the directory
   end
   FileNames = FileNames(contains(FileNames, '_PM25_speciation') | contains(FileNames, '_RCFM'));    %filter file names for those with raw data from the desired year
   numFiles = size(FileNames,1);
   warning('off','all')
%Load Data
   for i=1:numFiles
       filename = strcat(DataFileLoc,char(FileNames(i)));
       Data = readtable(filename);
       VarNames = Data.Properties.VariableNames; %table2cell(Data(1,:));
       Description_Keys = {'Site','Lat','Lon','Elev','Start_Year','Start_Month','Start_Day','Start_hour',...
           'End_Year','End_Month','End_Day','End_hour','Hours_sampled','Parameter_Code','Value','Method'};
       startpos = zeros(1,numel(Description_Keys));
       for j=1:numel(Description_Keys)
        [~,ipos] = max(contains(VarNames,Description_Keys{j}));           % find the lines that gives the site latitude, longitude, and elevations
        startpos(j)=ipos;
       end
       SiteCode = table2array(Data(2:end,startpos(1)));
       SiteLat = table2array(Data(2:end,startpos(2)));SiteLat(contains(SiteCode,'Ambient local'))=[];
       SiteLon = table2array(Data(2:end,startpos(3)));SiteLon(contains(SiteCode,'Ambient local'))=[];
       SiteElev = table2array(Data(2:end,startpos(4)));SiteElev(contains(SiteCode,'Ambient local'))=[];
       sYear = table2array(Data(2:end,startpos(5)));sYear(contains(SiteCode,'Ambient local'))=[];
       sMonth = table2array(Data(2:end,startpos(6)));sMonth(contains(SiteCode,'Ambient local'))=[];
       sDay = table2array(Data(2:end,startpos(7)));sDay(contains(SiteCode,'Ambient local'))=[];
       eYear = table2array(Data(2:end,startpos(9)));eYear(contains(SiteCode,'Ambient local'))=[];
       eMonth = table2array(Data(2:end,startpos(10)));eMonth(contains(SiteCode,'Ambient local'))=[];
       eDay = table2array(Data(2:end,startpos(11)));eDay(contains(SiteCode,'Ambient local'))=[];
       hr_sampled = table2array(Data(2:end,startpos(13)));hr_sampled(contains(SiteCode,'Ambient local'))=[];
       Param_Code = table2array(Data(2:end,startpos(14)));Param_Code(contains(SiteCode,'Ambient local'))=[];
       Value = table2array(Data(2:end,startpos(15)));Value(contains(SiteCode,'Ambient local'))=[];
       Method = table2array(Data(2:end,startpos(16)));Method(contains(SiteCode,'Ambient local'))=[];
       SiteCode(contains(SiteCode,'Ambient local'))=[];
       
   if i==1 
       SiteTmp = SiteCode;
       LatTmp = SiteLat;
       LonTmp = SiteLon;
       ElevTmp = SiteElev;
       sYearTmp = sYear;
       sMonthTmp = sMonth;
       sDayTmp = sDay;
       eYearTmp = eYear;
       eMonthTmp = eMonth;
       eDayTmp = eDay;
       hr_sampledTmp = hr_sampled;
       Param_CodeTmp = Param_Code;
       ValueTmp = Value;
       MethodTmp = Method;
   else
       SiteTmp = cat(1,SiteTmp,SiteCode);
       LatTmp = cat(1,LatTmp,SiteLat);
       LonTmp = cat(1,LonTmp,SiteLon);
       ElevTmp = cat(1,ElevTmp,SiteElev);
       sYearTmp = cat(1,sYearTmp,sYear);
       sMonthTmp = cat(1,sMonthTmp,sMonth);
       sDayTmp = cat(1,sDayTmp,sDay);
       eYearTmp = cat(1,eYearTmp,eYear);
       eMonthTmp = cat(1,eMonthTmp,eMonth);
       eDayTmp = cat(1,eDayTmp,eDay);
       hr_sampledTmp = cat(1,hr_sampledTmp,hr_sampled);
       Param_CodeTmp = cat(1,Param_CodeTmp,Param_Code);
       ValueTmp = cat(1,ValueTmp,Value);
       MethodTmp = cat(1,MethodTmp,Method);
   end
   fprintf('%s: added.\n',filename);    
   end
   
   NaNLoc = sYearTmp~=year | hr_sampledTmp <24; %find locations not measured for 24 hours and not in the given year
   SiteTmp(NaNLoc)= [];
   LatTmp(NaNLoc)= [];
   LonTmp(NaNLoc)= [];
   ElevTmp(NaNLoc)= [];
   sYearTmp(NaNLoc)= [];
   sMonthTmp(NaNLoc)= [];
   sDayTmp(NaNLoc)= [];
   eYearTmp(NaNLoc)= [];
   eMonthTmp(NaNLoc)= [];
   eDayTmp(NaNLoc)= [];
   Param_CodeTmp(NaNLoc)= [];
   ValueTmp(NaNLoc)= [];
   
   %Calc mid time
   Date = NaN(size(sYearTmp,1),1);
   for i=1:size(sYearTmp,1)
       midyear = eYearTmp(i)-sYearTmp(i);
       if midyear ~=0   %skip data point if spans multiple years
           midmonth = 0;midday=0;
       else
       midmonth = eMonthTmp(i)-sMonthTmp(i);
       if midmonth ~=0  % if the data spans two months, assign to month with most values
           if sMonthTmp(i)==1 || sMonthTmp(i)==3 || sMonthTmp(i)==5|| sMonthTmp(i)==7 ||...
                   sMonthTmp(i)==8||sMonthTmp(i)==10||sMonthTmp(i)==12
               dnum1 = 31-sMonthTmp(i)/2;
           elseif sMonthTmp(i)==2
               dnum1 = 28-sMonthTmp(i)/2;
           else
               dnum1 = 30-sMonthTmp(i);
           end
           dnum2 = eMonthTmp(i);
           if dnum1 > dnum2
               midmonth = sMonthTmp(i);
               midday = sDayTmp(i);
           else
               midmonth = eMonthTmp(i);
               midday = eDayTmp(i);
           end
       else
           midmonth = round(sMonthTmp(i));
           midday = round(sDayTmp(i)+((eDayTmp(i)-sDayTmp(i))/2));
       end
       end
       if midmonth<10;midmonth = strcat('0',string(midmonth));end
       if midday<10;midday = strcat('0',string(midday));end
       Date(i) = strcat(string(sYearTmp(i)),string(midmonth),string(midday));
   end
   
   MonStr = {strcat(YearStr,'01'), strcat(YearStr,'02'), strcat(YearStr,'03'), strcat(YearStr,'04'), ...
       strcat(YearStr,'05'), strcat(YearStr,'06'), strcat(YearStr,'07'), strcat(YearStr,'08'), ...
       strcat(YearStr,'09'), strcat(YearStr,'10'), strcat(YearStr,'11'), strcat(YearStr,'12')};
   MonStr = cellstr(MonStr)';
   SiteStr = unique(SiteTmp,'stable');
   VarStr = {'28101','28401','28402','28802','28201','28306','28801','28804','28805','28902','28904','28908','28305','28304'}; %'PM, SO4, NIT, AM, BC,OC, Na, Mg,Ca,Al,Ti,Fe,Calc_Dust,Calc_SS
   %PM2.5, SO4, NO3, NH4, BC, Residual Matter (reconstructed), Na, Mg, Ca, Al, Ti, Fe, Fine soil (reconstructed), sea salt (reconstructed)
   numVar = size(VarStr,2);
   numsites = size(SiteStr,1);
   SPARTAN_MonAvg_SiteInfo=repmat({NaN},[numsites,4]);                         %make placeholder for SIte Info wave (site X 3 (lat, lon, elev))
 
   SPARTAN_MonAvg_ugsm3=NaN(numsites,12,numVar);
   for i=1:numsites
       tmp = contains(SiteTmp,SiteStr{i});                      %find all points for each particular site
       SPARTAN_MonAvg_SiteInfo(i,1) = num2cell(unique(LatTmp(tmp),'stable'));%should give single value if all lats, lons, and elevations are reported the same
       SPARTAN_MonAvg_SiteInfo(i,2) = num2cell(unique(LonTmp(tmp),'stable'));
       SPARTAN_MonAvg_SiteInfo(i,3) = num2cell(unique(ElevTmp(tmp),'stable'));
       SPARTAN_MonAvg_SiteInfo(i,4) = unique(SiteTmp(tmp),'stable');
       if size(unique(LatTmp(tmp),'stable')) > 1 | size(unique(LonTmp(tmp),'stable')) > 1 | size(unique(ElevTmp(tmp),'stable')) > 1
           disp('!!!!!!!!SPARTAN LOAD ERROR!!!!!!!!!')
       end
       for j=1:12
        filtertmp = contains(SiteTmp,SiteStr{i}) & contains(string(Date), MonStr(j));
        for k=1:numVar
            %filter for values for this site, at this month, for each
            %desired variable and set the MonAvg array accordingly
            tmp = ValueTmp(filtertmp & strcmp(string(Param_CodeTmp), VarStr(k)));% & Methodstmp);
            tmp(tmp<0) = NaN;
            if sum(~isnan(tmp),1) >= 2                          %only record data if >= 2 measurements per month
                SPARTAN_MonAvg_ugsm3(i,j,k) = nanmean(tmp);         %take average of all month data from each site for each variable
            else
                SPARTAN_MonAvg_ugsm3(i,j,k) = NaN;
            end
        end
       end
   end
   savefile = sprintf('/misc/data6/emcduffie/ExtData/SPARTAN_Data/SPARTAN_MonthlyAvgs2_%s.mat', YearStr);
   save(savefile,'SPARTAN_MonAvg_ugsm3','SPARTAN_MonAvg_SiteInfo');
   
end
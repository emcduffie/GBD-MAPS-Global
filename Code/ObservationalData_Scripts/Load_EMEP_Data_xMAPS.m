function [EMEP_MonAvg_ugsm3, EMEP_MonAvg_SiteInfo] = Load_EMEP_Data_xMAPS(DataFileLoc,year)
 %load EMEP Data from data file source
% E. McDuffie - Oct. 16, 2020
%INPUTS:
%   Data File - EMEP input data file location
%   Year - desired data year
%OUTPUTS:
%   EMEP_MonAvg_ugm3 is a matrix of the monthly average components at each site. 
%   EMEP_MonAvg_SiteInfo  - lat, lon, elevation of each site
% Methods selsected for each component are provided below
 % ***Assume PM2.5 data are reported at 35% RH
 % ***Species are reported as dry concentrations
 % ***Assume all data are reported at local conditions (no confirmation available)
%INPUT DATA SOURCE: http://ebas.nilu.no, 
%   selecting all available countries, instrument types, and components associated with the ‘pm25’ matrix

%Additional Notes:
%   PM2.5 measurements have been filtered to exclude beta-attenuation methods.
%   Analysis methods for organic and black carbon are not frequently reported.
%   Differences in OC and BC concentrations derived from different methods has 
%   been previously reported for US networks. It is unknown whether the EMEP data
%   are reported at ambient conditions or whether they have been corrected using
%   standard temperature and pressure. Here, we assume data are reported in 
%   local conditions for consistency in units with national ambient air quality standards. 
  
	YearStr = string(year);
	DataFileLoc = sprintf('%sEMEP_%d/', DataFileLoc,year);
	listing = dir(DataFileLoc);
	numFiles = size(listing,1);
	FileNames = cell(numFiles,1);
    for i=1:numFiles
        FileNames(i,1) = cellstr(listing(i).name);   %makes cell array of the file names in the directory
    end
	FileNames = FileNames(contains(FileNames, '.nas'));    %filter file names for nas files
	numFiles = size(FileNames,1);
	EMEPSite = repmat({''},numFiles*3000,1);
	EMEPSiteLat = repmat(NaN,numFiles*3000,1);
	EMEPSiteLon = EMEPSiteLat;
    EMEPSiteElev = EMEPSiteLat;
    EMEPDate = EMEPSiteLat;
    EMEPVarCode = EMEPSite;
    EMEPUnitCode = EMEPSite;
    EMEPValue = double(EMEPSiteLat);
    EMEPVarTxt = EMEPSite;
    EMEPInstType = EMEPSite;
    EMEPInstName = EMEPSite;
    EMEPMethodRef = EMEPSite;
   
    AllDates =[];
    AllDates = [AllDates Juln2Grg((year*1000+1):Grg2Juln(year*10000+1231))]'; %make 1-year timewave by day
    counter = 1;
   
    for i=1:numFiles
        filename = strcat(DataFileLoc,char(FileNames(i)));
        disp(filename)
        V = ICARTTreader(filename);
        numvar = size(V.varnames,1);
        Resolution = V.header(contains(V.header,'Resolution Code:','IgnoreCase',true));
        if ~isempty(Resolution)
            Resolution = char(Resolution);
            Resolution = strtrim(Resolution(17:end));
        else
            disp('***UNKNOWN Data Resolution: ')
            disp(filename)
            Resolution = 'NaN';
        end
       
        if strcmp(Resolution,'1d') || strcmp(Resolution, '15h')
            offset = 1;%+round(V.starttime(1)); %figure out what day this is closest
            numdays = size(V.starttime,1);
            for j=1:numdays                  % for each day and variable, record the site info, date, variable info and value
                itime = round(V.starttime(j))+offset; 
                if itime <= numdays && itime > 0       %skip the last value if >days in year and skip the first if from the prior year
                    for k=3:numvar                      %var1 and var2 are start and endtime
                        vartmp = char(V.varnames(k));
                        EMEPSite(counter) = V.SiteInfo{1};
                        EMEPSiteLat(counter) = V.SiteInfo{2};
                        EMEPSiteLon(counter) = V.SiteInfo{3};
                        EMEPSiteElev(counter) = V.SiteInfo{4};
                        EMEPDate(counter) = AllDates(itime);
                        EMEPVarCode(counter) = V.varnames(k);  
                        EMEPUnitCode(counter) = V.units(k);
                        EMEPValue(counter) = V.(vartmp)(j);     %get value at each time point
                        EMEPVarTxt(counter) = V.vartxt(k);
                        EMEPInstType(counter) = V.insttype(1);
                        EMEPInstName(counter) = V.instname(1);
                        EMEPMethodRef(counter) = V.methodref(1);
                        counter = counter+1;
                    end
                end
            end
        elseif strcmp(Resolution,'1h')
            offset = 1;
            numpnts = size(V.starttime,1);
            idx1=1;
            while idx1 < numpnts
                itime = round(V.starttime(idx1))+offset;
                istop = floor(V.endtime(idx1))+offset;
                idx2 = idx1;
                while istop ==itime
                    idx2=idx2+1;
                    if idx2 < numpnts                           %check whether we've reached the end of the file
                        istop = floor(V.endtime(idx2))+offset;   %find the point number where the day changes
                    else
                        istop = idx2;                            %set to the end of the file if we've reached the end
                    end
                end
                if itime < size(AllDates,1) && itime > 0
                    for k=3:numvar                          %var1 and var2 are start and endtimes
                        vartmp = char(V.varnames(k));
                        EMEPSite(counter) = V.SiteInfo{1};
                        EMEPSiteLat(counter) = V.SiteInfo{2};
                        EMEPSiteLon(counter) = V.SiteInfo{3};
                        EMEPSiteElev(counter) = V.SiteInfo{4};
                        EMEPDate(counter) = AllDates(itime);
                        EMEPVarCode(counter) = V.varnames(k);  
                        EMEPUnitCode(counter) = V.units(k);
                        EMEPValue(counter) = nanmean(V.(vartmp)(idx1:idx2));     %get value at each time point (average over a day)
                        EMEPVarTxt(counter) = V.vartxt(k);
                        EMEPInstType(counter) = V.insttype(1);
                        EMEPInstName(counter) = V.instname(1);
                        EMEPMethodRef(counter) = V.methodref(1);
                        counter = counter+1;
                    end
                end
                idx1=idx2+1;
            end
        elseif strcmp(Resolution,'6d') || strcmp(Resolution, '4d') || strcmp(Resolution, '1w') || strcmp(Resolution, '1mo') || strcmp(Resolution, '2w') ||strcmp(Resolution, '2d')
            offset = 1;%+round(V.starttime(1));
            numpnts = size(V.starttime,1);
            for j=1:numpnts
                if V.starttime(j) < 0
                    itime = round(round(V.starttime(j)) + (V.endtime(j) - V.starttime(j))/2)+offset; %find the middle date of the week
                    if round(V.endtime(j) - V.starttime(j)) ==1
                        itime = round(V.starttime(j))+offset;
                    end
                else
                    itime = round(round(V.starttime(j)) + (V.endtime(j) - V.starttime(j))/2)+offset; %find the middle date of the week
                    if round((V.endtime(j) - V.starttime(j))) ==1        %check to see whether the data duration is actually 1 day
                        itime = round(V.starttime(j))+offset;
                    end
                end
                if itime < size(AllDates,1) && itime > 0
                    for k=3:numvar                          %var1 and var2 are start and endtimes
                        vartmp = char(V.varnames(k));
                        EMEPSite(counter) = V.SiteInfo{1};
                        EMEPSiteLat(counter) = V.SiteInfo{2};
                        EMEPSiteLon(counter) = V.SiteInfo{3};
                        EMEPSiteElev(counter) = V.SiteInfo{4};
                        EMEPDate(counter) = AllDates(itime);
                        EMEPVarCode(counter) = V.varnames(k);  
                        EMEPUnitCode(counter) = V.units(k);
                        EMEPValue(counter) = V.(vartmp)(j);     %get value at each time point
                        EMEPVarTxt(counter) = V.vartxt(k);
                        EMEPInstType(counter) = V.insttype(1);
                        EMEPInstName(counter) = V.instname(1);
                        EMEPMethodRef(counter) = V.methodref(1);
                        counter = counter+1;
                    end
                end
            end
        else
            disp('**Alternative Resolution')
            disp(filename)
        end
    end                        %load all Files
    EMEPSite = EMEPSite(1:counter-1);             %reformat to remove extra nans
    EMEPSiteLat = EMEPSiteLat(1:counter-1);
    EMEPSiteLon = EMEPSiteLon(1:counter-1);
    EMEPSiteElev = EMEPSiteElev(1:counter-1);
    EMEPDate = EMEPDate(1:counter-1);
    EMEPVarCode = EMEPVarCode(1:counter-1);
    EMEPUnitCode = EMEPUnitCode(1:counter-1);
    EMEPValue = EMEPValue(1:counter-1);
    EMEPVarTxt = EMEPVarTxt(1:counter-1);
    EMEPInstType = EMEPInstType(1:counter-1);
    EMEPInstName = EMEPInstName(1:counter-1);
    EMEPMethodRef = EMEPMethodRef(1:counter-1);

    S.Site = EMEPSite(~contains(EMEPVarCode, 'numflag'));                          %make structure with all the data, eliminate 'data flags'
	S.Elev = single(EMEPSiteElev(~contains(EMEPVarCode, 'numflag')));
	S.Lat = single(EMEPSiteLat(~contains(EMEPVarCode, 'numflag')));
	S.Lon = single(EMEPSiteLon(~contains(EMEPVarCode, 'numflag')));
    S.Date = string(EMEPDate(~contains(EMEPVarCode, 'numflag')));
    S.Value = single(EMEPValue(~contains(EMEPVarCode, 'numflag')));
    S.Unit = EMEPUnitCode(~contains(EMEPVarCode, 'numflag'));
    S.VarTxt = EMEPVarTxt(~contains(EMEPVarCode, 'numflag'));
    S.InstName = EMEPInstName(~contains(EMEPVarCode, 'numflag'));
	S.InstType = EMEPInstType(~contains(EMEPVarCode, 'numflag'));
	S.Method = EMEPMethodRef(~contains(EMEPVarCode, 'numflag'));
    S.VarCode = EMEPVarCode(~contains(EMEPVarCode, 'numflag'));
    
    clear listing AllDates counter filename FileNames i idx1 idx2 istop itime...
            i j numdays numFiles numpnts numvar numVar offset Resolution SiteStr...
            Units V Varnames VarStr vartmp k EMEPSite EMEPSiteElev EMEPSiteLat...
        EMEPSiteLon EMEPDate EMEPValue EMEPVarCode EMEPUnitCode EMEPVarTxt...
        EMEPInstName EMEPInstType EMEPMethodRef

    S.MDL = S.Lat;
    for i=1:size(S.Lat,1)
        check = strcmp(S.VarTxt(i),'Detection limit'); %check if MDL is reported
        if check~=0
            if contains(S.Unit(i),'ug')
                S.MDL(i) = str2double(extractBetween(S.VarTxt(i),'Detection limit=','u'));
            elseif contains(S.Unit(i), 'ng')
                S.MDL(i) = str2double(extractBetween(S.VarTxt(i),'Detection limit=','n')); 
            else
                if contains(S.VarCode(i), 'numflag') %ignore the 'numflag' data
                    S.MDL(i) = 0;
                else
                    disp('**UNKNOWN MDL Unit')
                    disp(i)                              %report any other cases
                    disp(S.Unit(i))
                end
            end
        else
            S.MDL(i) = 0;
        end
    end
   
    %filter for data below the DL's (if given)
    S.Value(S.Value < S.MDL) = NaN;
    %convert from ug of N and S to NO3, NH4, and SO4
    for i=1:size(S.Lat)
        if contains(S.Unit(i),'N/')
            if strcmp(S.VarCode(i),'nitrate')
                S.Value(i) = S.Value(i).*(62/14); 
            elseif strcmp(S.VarCode(i),'ammonium')
                S.Value(i) = S.Value(i).*(18/14);
            end
            if contains(S.Unit(i),'ng')
                S.Value(i) = S.Value(i) ./ 1e3;
            end
        elseif contains(S.Unit(i),'S/')
            S.Value(i) = S.Value(i).*(96/32);
            if contains(S.Unit(i),'ng')
                S.Value(i) = S.Value(i)./1e3;
            end
        elseif contains(S.Unit(i),'ng')
            S.Value(i) = S.Value(i)./1e3;
        elseif ~contains(S.Unit(i),'C/') && ~contains(S.Unit(i),'ug/m3')&& ~contains(S.Unit(i),'hPa')&& ~contains(S.Unit(i),'K')
            disp('**Unknown Variable Unit')
            disp(S.Unit(i))
        end
    end

    %Make monthly averages for each site (variable avg: site x month x variable )
    MonStr = {strcat(YearStr,'01'), strcat(YearStr,'02'), strcat(YearStr,'03'), strcat(YearStr,'04'), ...
       strcat(YearStr,'05'), strcat(YearStr,'06'), strcat(YearStr,'07'), strcat(YearStr,'08'), ...
       strcat(YearStr,'09'), strcat(YearStr,'10'), strcat(YearStr,'11'), strcat(YearStr,'12')};
    MonStr = cellstr(MonStr)';
    SiteStr = unique(S.Site,'stable');
    VarStr = {'pm25_mass','nitrate','ammonium','sulphate_corrected','sulphate_total','elemental_carbon','equivalent_black_carbon','organic_carbon','total_carbon','chloride','aluminium','calcium','iron','titanium', 'magnesium', 'sodium'};
    numVar = size(VarStr,2);

    %set the logic waves for the methods filters to apply for each
    % variable at each site below.s
    MethodsFilter = cell(numVar,1);
    for ivar=1:numVar
       %var = VarStr(ivar);
       switch true
            %char(VarStr(ivar))
            case strcmp(char(VarStr(ivar)),'pm25_mass') %extract gravimetric, 24-hour methods,  %%can change this to inlcude beta-attneuation at a later date, *all VarTxt is OK
                MethodsFilter{ivar,1} = (~contains(S.InstType, 'b-atten','IgnoreCase',true) & ~contains(S.InstType, 'beta', 'IgnoreCase',true));% & contains(EMEPMethodRef,'24 hour','IgnoreCase',true);
            %extract all measurements (they all seem to use IC), but could
            %filter for 'IC' in S.Method, S.VarCode, S.VarTxt just to be
            %sure - unsure about filter2pk and filter3pk methods
            case strcmp(char(VarStr(ivar)), 'nitrate') | strcmp(char(VarStr(ivar)), 'ammonium') | strcmp(char(VarStr(ivar)), 'sulphate_corrected')...
                    |strcmp(char(VarStr(ivar)),'sulphate_total') | strcmp(char(VarStr(ivar)),'chloride')  %all nitrate vartxt OK
                MethodsFilter{ivar,1} = logical(repmat(1,[size(S.Value),1]));%(contains(S.InstType,'IC','IgnoreCase',true) | contains(S.InstName,'IC','IgnoreCase',true) | contains(S.Method,'IC','IgnoreCase',true));% ones(size(S.Value,1));
            case strcmp(char(VarStr(ivar)),'elemental_carbon') %dont include data that are given as uncertainties - have no idea what the methods are
                MethodsFilter{ivar,1} = ~contains(S.VarTxt,'Uncertainty','IgnoreCase',true);
            case strcmp(char(VarStr(ivar)),'equivalent_black_carbon') %dont include these for now
                MethodsFilter{ivar,1} = logical(repmat(1,[size(S.Value),1]));
            case strcmp(char(VarStr(ivar)),'organic_carbon') %no idea what method
                MethodsFilter{ivar,1} = ~(contains(S.VarTxt,'Uncertainty','IgnoreCase',true)| contains(S.VarTxt,'Fraction','IgnoreCase',true)|contains(S.VarTxt,'Artifact','IgnoreCase',true));   
            case strcmp(char(VarStr(ivar)),'total_carbon') %no idea what methods
                MethodsFilter{ivar,1} = ~contains(S.VarTxt,'Uncertainty','IgnoreCase',true);
            case strcmp(char(VarStr(ivar)),'aluminium') | strcmp(char(VarStr(ivar)),'calcium') | strcmp(char(VarStr(ivar)),'iron') | strcmp(char(VarStr(ivar)),'titanium')|...
                    strcmp(char(VarStr(ivar)),'magnesium') | strcmp(char(VarStr(ivar)),'sodium') %these measruements seem to be ICPMS, unclear of all though
                MethodsFilter{ivar,1} = logical(repmat(1,[size(S.Value),1]));
       end
    end
    numsites = size(SiteStr,1);
    EMEP_MonAvg_SiteInfo=repmat({NaN},[numsites,4]); %make placeholder for SIte Info wave (site X 3 (lat, lon, elev))
    EMEP_MonAvg_ugsm3 = repmat(NaN,[numsites,12,numVar]);
 
    for i=1:numsites
        %Make waves of site lat, lon, elevation
        tmp = strcmp(S.Site,SiteStr(i)); %find all points for each particular site
        EMEP_MonAvg_SiteInfo(i,1) = num2cell(unique(S.Lat(tmp),'stable')); %should give single value if all lats, lons, and elevations are reported the same
        EMEP_MonAvg_SiteInfo(i,2) = num2cell(unique(S.Lon(tmp),'stable'));
        EMEP_MonAvg_SiteInfo(i,3) = num2cell(unique(S.Elev(tmp),'stable'));
        EMEP_MonAvg_SiteInfo(i,4) = unique(S.Site(tmp),'stable');
        if size(unique(S.Lat(tmp),'stable')) > 1 | size(unique(S.Lon(tmp),'stable')) > 1 | size(unique(S.Elev(tmp),'stable')) > 1| size(unique(S.Site(tmp),'stable')) > 1
           disp('!!!!!!!!EMEP LOAD ERROR!!!!!!!!!')
        end
        for j=1:12
        filtertmp = contains(S.Site,SiteStr(i)) & contains(S.Date, MonStr(j));
        for k=1:numVar
            %filter for values for this site, at this month, for each
            %desired variable and set the MonAvg array accordingly
            Methodstmp = MethodsFilter{k,1};
            tmp = S.Value(filtertmp & strcmp(S.VarCode, VarStr(k)) & Methodstmp);
            tmp(tmp<0) = NaN;
           % if contains(MonStr(j),'02') 
           %     dayfilter2=8;
           %     else;dayfilter2=dayfilter;
           % end
            if sum(~isnan(tmp),1) >= 4%dayfilter                              %only record data if >= 4 measurements per month
                EMEP_MonAvg_ugsm3(i,j,k) = nanmean(tmp);   %take average of all month data from each site for each variable
            else
                EMEP_MonAvg_ugsm3(i,j,k) = NaN;
            end
        end
       end
    end
   savestr = sprintf('/misc/data6/emcduffie/ExtData/EMEP_Data/EMEP/EMEP_MonthlyAvgs_%s.mat',YearStr);
   save(savestr,'EMEP_MonAvg_ugsm3','EMEP_MonAvg_SiteInfo');
   
end
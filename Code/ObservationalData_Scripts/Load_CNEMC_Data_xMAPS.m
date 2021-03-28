function [CNEMC_MonAvg_ugsm3, CNEMC_MonAvg_SiteInfo] = Load_CNEMC_Data_xMAPS(DataFileLoc,year)
%Function to load CNEMC data. Assume data are reported at local conditions 
% E. McDuffie - Oct. 16, 2020
%INPUTS:
%   Data File - CNEMC input data file location
%   Year - desired data year
%OUTPUTS:
%   CNEMC_MonAvg_ugm3 is a matrix of the monthly average PM2.5 each site. 
%   CNEMC_MonAvg_SiteInfo  - lat, lon, elevation of each site
% *** Assume PM2.5 data is at 35% RH and in volumetric units (no confirmation available) ***
%INPUT DATA SOURCE: http://www.cnemc.cn/en/ 


    YearStr = string(year);
   
    tAllDates = [];
    tAllDates = [tAllDates Juln2Grg((str2double(YearStr)*1000+1):Grg2Juln(str2double(YearStr)*10000+1231))];

    SiteInfoFile = [DataFileLoc 'China_Sites.csv'];
    SiteID = [];
    Latitude = [];
    Longitude = [];
    
    fid = fopen(SiteInfoFile);
    header = fgetl(fid);
    while(~feof(fid))
        fline = fgetl(fid);
        cspot = strfind(fline,',');
        SiteID = [SiteID; str2double(fline(1:4))];
        Longitude = [Longitude; str2double(fline(cspot(3)+1:cspot(4)-1))];
        Latitude = [Latitude; str2double(fline(cspot(4)+1:end))];
    end
    fclose(fid);
    
    numsites = size(SiteID,1);
    CNEMC_MonAvg_SiteInfo = NaN(numsites,3);  %make wave of sites X 3 (lat, lon, elev)
    CNEMC_MonAvg_ugsm3 = NaN(numsites,12); %monthly data
    CNEMC_MonAvg_SiteInfo(:,1) = SiteID;
    CNEMC_MonAvg_SiteInfo(:,2) = Latitude;
    CNEMC_MonAvg_SiteInfo(:,3) = Longitude;
    
    tPM25Data = NaN(numsites,numel(tAllDates),24);
    tdir = dir([DataFileLoc sprintf('China_%d',str2double(YearStr)) '*']);
    numdir = numel(tdir);
    for j = 2:numdir %**hard coded to avoid zip files, remove in the future
        if tdir(j).isdir == 1
            tfiles = dir([DataFileLoc tdir(j).name '/*.csv']);
            for k = 1:numel(tfiles)
                tdatespot = find(tAllDates == str2double(tfiles(k).name(end-11:end-4)));
                fid = fopen([DataFileLoc tdir(j).name '/' tfiles(k).name]);
                header = fgetl(fid);
                cspot = strfind(header,',');
                if ~isempty(cspot)
                    typespot = 2;
                    hrspot = 1;
                    tSiteID = header(cspot(typespot+1)+1:end);
                    tSiteID(strfind(tSiteID,'A')) = '';
                    tSiteID = str2num(tSiteID)';

                    tSiteIDspot = NaN(numel(tSiteID),1);
                    for x = 1:numel(tSiteID)
                        t = find(SiteID == tSiteID(x));
                        if ~isempty(t)
                            tSiteIDspot(x) = t;
                        end
                    end
                    tSiteIDspoti = find(~isnan(tSiteIDspot));

                    while(~feof(fid))
                        fline = fgetl(fid);
                        cspot = strfind(fline,',');
                        if strcmp(fline(cspot(typespot)+1:cspot(typespot+1)-1),'PM2.5')
                            tHr = str2double(fline(cspot(hrspot)+1:cspot(hrspot+1)-1));
                            if tHr == 0
                                tHr = 24;
                            end
                            tPM25 = NaN(numel(tSiteID),1);
                            for x = (typespot+1):numel(cspot)
                                if (x+1 <= numel(cspot)) && cspot(x)+1 ~= cspot(x+1)
                                    tPM25(x-typespot) = str2double(fline(cspot(x)+1:cspot(x+1)-1));
                                elseif (x == numel(cspot))
                                    t = str2double(fline(cspot(x)+1:end));
                                    if ~isempty(t)
                                        tPM25(x-typespot) = t;
                                    end
                                end
                            end
                            tPM25Data(tSiteIDspot(tSiteIDspoti),tdatespot,tHr) = tPM25(tSiteIDspoti);
                        end
                    end
                end
                
                fclose(fid);
                fprintf('%s complete\n',[DataFileLoc tdir(j).name '/' tfiles(k).name]);
            end
            %break
        end
        MonStr = {strcat(YearStr,'01'), strcat(YearStr,'02'), strcat(YearStr,'03'), strcat(YearStr,'04'), ...
        strcat(YearStr,'05'), strcat(YearStr,'06'), strcat(YearStr,'07'), strcat(YearStr,'08'), ...
        strcat(YearStr,'09'), strcat(YearStr,'10'), strcat(YearStr,'11'), strcat(YearStr,'12')};
    
        %disp('test')
        for i=1:numsites
            for j=1:12
            filtertmp = contains(string(tAllDates)', MonStr{j});
                %filter for values for this site, at this month, for each
                %desired variable and set the MonAvg array accordingly
                tmp = tPM25Data(i,filtertmp,:);
                tmp(tmp<0) = NaN;
                daysum = nansum(~isnan(tmp(1,:,:)),3); %sum number of hourly data points
                tmp2=NaN(size(daysum,2),1);
                for m=1:size(daysum,2)
                if daysum(m) >=24                           %if at least 24 data points...(must have all daily values)
                    tmp2(m) = squeeze(nanmean(tmp(1,m,1:24),3))';
                else
                    tmp2(m) = NaN;
                end
                end
                if sum(~isnan(tmp2),1) >= 20              % if at least 20 sampling days during month
                    CNEMC_MonAvg_ugsm3(i,j) = nanmean(tmp2);   %take average of all month data from each site for each variable
                else
                    CNEMC_MonAvg_ugsm3(i,j) = NaN;
                end
            end
        end
    end
   savefile = sprintf('/misc/data6/emcduffie/ExtData/CNEMC_Data/CNEMC_MonthlyAvgs_%s.mat',YearStr);
   save(savefile,'CNEMC_MonAvg_ugsm3','CNEMC_MonAvg_SiteInfo');
end
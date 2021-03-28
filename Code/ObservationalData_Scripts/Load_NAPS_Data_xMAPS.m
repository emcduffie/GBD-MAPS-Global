function [NAPS_MonAvg_ugsm3, NAPS_MonAvg_SiteInfo] = Load_NAPS_Data_xMAPS(DataFileLoc,year)
%Load NAPS Data. Raw data are reported at local conditions 
% E. McDuffie - Oct. 16, 2020
%INPUTS:
%   Data File - NAPS input data file location
%   Year - desired data year
%OUTPUTS:
%   NAPS_MonAvg_ugm3 is a matrix of the monthly average components at each site. 
%   NAPS_MonAvg_SiteInfo  - lat, lon, elevation of each site
% Data have been filtered for values above DL.
% Data are reported at local conditions - 35% RH for PM, dry for others
%INPUT DATA SOURCE: 
%   Hourly and Integrated files: http://data.ec.gc.ca/data/air/monitor/national-air-pollution-surveillance-naps-program/Data-Donnees/2017/?lang=en

%Additional Notes:
%   In 2017, hourly measurements are reported from a variety of instruments...
%   including the Tapered Element Oscillating Microbalance (TEOM), Scientific
%   Synchronized Hybrid Ambient Real-time Particulate (SHARP) model 5030, 
%   and Met-One Beta-Attenuation Monitors (BAM). TEOM instruments heat the 
%   sample line to 40C to reduce particle bound water, but result in a loss of
%   volatile species, particularly during cold seasons, leading to the phase-out 
%   of these instruments in some provinces (40). Therefore, the TEOM data from the
%   NAPS network (Method codes 706 and 760) have been corrected to account for 
%   differences between these measurements and FRM methods during warm and cold 
%   seasons across Canada, based on recommendations from Health Canada in 2013. 
%   SHARP data (Method code 184) have similarly been corrected. BAM monitors also %
%   have uncertainties associated with an over-measurement of aerosol water relative 
%   to FRM methods at high ambient RH, though this uncertainty has been reduced in 
%   NAPS data with the installation of heaters that control sample flow RH at 
%   35% and BAM data (Method code 170) are not corrected in this work.

%%%************************************
   Hourly =1; %Load hourly NAPS data? (these need met data for the conversion and a transformation factor (which has uncertainties))
   ApplyTransformations = 1; %make corrections to hourly data based on email from A. van Donkelaar (1 = yes)
           
   YearStr = string(year);
   listing = dir(DataFileLoc);
   numFiles = size(listing,1);
   FileNames = cell(numFiles,1);
   for i=1:numFiles
       FileNames(i,1) = cellstr(listing(i).name);   %makes cell array of the file names in the directory
   end
   FileNames = FileNames(contains(FileNames, '_PM25_') & contains(FileNames, YearStr));    %filter file names for those with raw data from the desired year
   numFiles = size(FileNames,1);
   orgcounter = 1;
   xrfcounter =1;
   iccounter = 1;
   warning('off','all')
   NAPS_MonAvg_SiteInfo = NaN(numFiles,4);
%Load Data
   for i=1:numFiles
       filename = strcat(DataFileLoc,char(FileNames(i)));
       disp(filename)
       MetaData = readtable(filename,'Sheet','Station Info');
       Descriptions = table2cell(MetaData(:,1));
       [~,iLat] = max(contains(Descriptions,'Lat'));           % find the lines that gives the site latitude, longitude, and elevations
       [~,iLon] = max(contains(Descriptions,'Lon'));
       [~,iElev] = max(contains(Descriptions,'Elev'));
       [~,iID] = max(contains(Descriptions,'NAPS ID'));
       SiteLat = table2array(MetaData(iLat,2));
       SiteLon = table2array(MetaData(iLon,2));
       SiteElev = table2array(MetaData(iElev,2));
       SiteID = table2array(MetaData(iID,2));
       NAPS_MonAvg_SiteInfo(i,1) = str2double(SiteLat); %should give single value if all lats, lons, and elevations are reported the same
       NAPS_MonAvg_SiteInfo(i,2) = str2double(SiteLon);
       NAPS_MonAvg_SiteInfo(i,3) = str2double(SiteElev);
       NAPS_MonAvg_SiteInfo(i,4) = str2double(SiteID);
       LoadOrganics = max(contains(table2cell(MetaData(:,2)),'OCEC & OC Artifacts')); %find out if organic analysis is included in data
       LoadXRF = max(contains(table2cell(MetaData(:,2)),'Elements_ED-XRF'));
       LoadIC = max(contains(table2cell(MetaData(:,2)),'Ions-Spec_IC'));
       
       %PM2.5 Data Load
       Data = readtable(filename,'Sheet','PM2.5','ReadVariableNames',false); %load data from PM2.5 sheet
       [~,irow] = max(contains(table2cell(Data(:,1)),'Sampler'));           %find the sampler header line   
       SamplerLine = table2cell(Data(irow,:));                            %record all the sampler headers
       SamplerLine(cellfun(@isempty,SamplerLine)) = {'NaN'};          % convert all non-character elements into characters
       Samplers = SamplerLine(~contains(SamplerLine,'NaN') & ~contains(SamplerLine,'Sampler')); %find all sample types
       Samplers = unique(Samplers,'stable');
       numSamplers = size(Samplers,2);                         %get the number of samplers
       [~,irow] = max(contains(table2cell(Data(:,1)),'NAPS Site ID'));   %find the data header line
       HeaderLine = table2cell(Data(irow,:));                                %record all the data headers
       HeaderLine(cellfun(@isempty,HeaderLine)) = {'NaN'};

       for j=1:numSamplers
            [~,icol] = max(contains(HeaderLine,'NAPS Site ID'));
            SiteTmp = table2array(Data(irow+1:end,icol));      %extract all the site code data
            [~,icol] = max(contains(HeaderLine,'Sampling Date'));
            DateTmp = table2cell(Data(irow+1:end,icol));       %extract all sample times
            [~,icol] = max(contains(HeaderLine,'Sample Type'));
            TypeTmp = table2cell(Data(irow+1:end,icol));        % extract all sample types
            [~,icol] = max(strcmp(HeaderLine, 'PM2.5') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
            DataTmp = table2array(Data(irow+1:end,icol));
            [~,icol] = max(contains(HeaderLine, 'PM2.5-MDL') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
            MDLTmp = table2array(Data(irow+1:end,icol));
            [~,icol] = max(contains(HeaderLine, 'PM2.5-Vflag') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
            VflagTmp = table2array(Data(irow+1:end,icol));
            [~,icol] = max(contains(HeaderLine, 'Temp.') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
            TTmp = table2array(Data(irow+1:end,icol));
            [~,icol] = max(contains(HeaderLine, 'Pres.') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
            PTmp = table2array(Data(irow+1:end,icol));
            [~,icol] = max(contains(HeaderLine, 'Actual Volume') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
            VTmp = table2array(Data(irow+1:end,icol));
            if j ==1
               SiteData = SiteTmp;
               DateData = DateTmp;
               TypeData = TypeTmp;
               PMData = DataTmp;
               MDLData = MDLTmp;
               VData = VflagTmp;
               PData = PTmp;
               TData = TTmp;
               VolData = VTmp;
            else
               SiteData = cat(1,SiteData,SiteTmp);  %concatenate data
               DateData = cat(1,DateData,DateTmp);
               TypeData = cat(1,TypeData,TypeTmp);
               PMData = cat(1,PMData,DataTmp);
               MDLData = cat(1,MDLData,MDLTmp);
               VData = cat(1,VData,VflagTmp);
               PData = cat(1,PData,PTmp);
               TData = cat(1,TData,TTmp);
               VolData = cat(1,VolData,VTmp);
            end
       end
       %filter for routine sampling values that are vaild or historical,
       %and above MDL
       PMData = str2double(PMData);
       MDLData = str2double(MDLData);
       VolData = str2double(VolData);
       PData = str2double(PData);
       TData = str2double(TData);
       PMData(~(contains(TypeData,'R') & (contains(VData,'H1') | contains(VData,'V0') | contains(VData,'' ) ) & PMData >= MDLData)) = NaN;
    
       %save to final waves that are collecting data from each site
       if i==1
           NAPS_PMData.Site = SiteData;
           NAPS_PMData.Date = DateData;
           NAPS_PMData.PM = PMData;
       else
           NAPS_PMData.Site = cat(1,NAPS_PMData(:).Site,SiteData);
           NAPS_PMData.Date = cat(1,NAPS_PMData(:).Date,DateData);
           NAPS_PMData.PM = cat(1,NAPS_PMData(:).PM,PMData);
       end
       
       %Load XRF Data
       if LoadXRF==1
            clear Data
            Data = readtable(filename,'Sheet','Elements_EDXRF','ReadVariableNames',false); %load data from XRF sheet 
            [~,irow] = max(contains(table2cell(Data(:,1)),'Sampler'));           %find the sampler header line
            SamplerLine = table2cell(Data(irow,:));                            %record all the sampler headers
            SamplerLine(cellfun(@isempty,SamplerLine)) = {'NaN'};          % convert all non-character elements into characters
            Samplers = SamplerLine(~contains(SamplerLine,'NaN') & ~contains(SamplerLine,'Sampler')); %find all sample types
            Samplers = unique(Samplers,'stable');
            numSamplers = size(Samplers,2);                         %get the number of samplers
            [~,irow] = max(contains(table2cell(Data(:,1)),'NAPS Site ID'));   %find the data header line
            HeaderLine = table2cell(Data(irow,:));                                %record all the data headers
            HeaderLine(cellfun(@isempty,HeaderLine)) = {'NaN'};
            
            Data2 = readtable(filename,'Sheet','PM2.5','ReadVariableNames',false); %load data from XRF sheet 
            [~,irow] = max(contains(table2cell(Data2(:,1)),'Sampler'));           %find the sampler header line
            SamplerLine2 = table2cell(Data2(irow,:));                            %record all the sampler headers
            SamplerLine2(cellfun(@isempty,SamplerLine2)) = {'NaN'};          % convert all non-character elements into characters
            Samplers2 = SamplerLine2(~contains(SamplerLine2,'NaN') & ~contains(SamplerLine2,'Sampler')); %find all sample types
            Samplers2 = unique(Samplers2,'stable');
            [~,irow] = max(contains(table2cell(Data2(:,1)),'NAPS Site ID'));   %find the data header line
            HeaderLine2 = table2cell(Data2(irow,:));                                %record all the data headers
            HeaderLine2(cellfun(@isempty,HeaderLine2)) = {'NaN'};
            
            for j=1:numSamplers
                [~,icol] = max(contains(HeaderLine,'NAPS Site ID'));
                SiteTmp = table2array(Data(irow+1:end,icol));      %extract all the site code data
                [~,icol] = max(contains(HeaderLine,'Sampling Date'));
                DateTmp = table2cell(Data(irow+1:end,icol));       %extract all sample times
                [~,icol] = max(contains(HeaderLine,'Sampling Type'));
                TypeTmp = table2cell(Data(irow+1:end,icol));        % extract all sample types
                [~,icol] = max(strcmp(HeaderLine, 'Aluminum (Al)') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
                AlTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'Al-MDL') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
                AlMDLTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'Al-VFlag') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
                AlVflagTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(strcmp(HeaderLine, 'Silicon (Si)') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
                SiTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'Si-MDL') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
                SiMDLTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'Si-VFlag') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
                SiVflagTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(strcmp(HeaderLine, 'Sulphur (S)') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
                STmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'S-MDL') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
                SMDLTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'S-VFlag') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
                SVflagTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(strcmp(HeaderLine, 'Calcium (Ca)') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
                CaTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'Ca-MDL') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
                CaMDLTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'Ca-VFlag') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
                CaVflagTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(strcmp(HeaderLine, 'Titanium (Ti)') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
                TiTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'Ti-MDL') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
                TiMDLTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'Ti-VFlag') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
                TiVflagTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(strcmp(HeaderLine, 'Iron (Fe)') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
                FeTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'Fe-MDL') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
                FeMDLTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'Fe-VFlag') & contains(SamplerLine, Samplers(j)));  %find the column of PM2.5 data from each sampler
                FeVflagTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(strcmp(HeaderLine2, 'Pres.') & contains(SamplerLine2, Samplers2(j)));  %find the column of PM2.5 data from each sampler
                PTmp = table2array(Data2(irow+1:end,icol));
                [~,icol] = max(strcmp(HeaderLine2, 'Temp.') & contains(SamplerLine2, Samplers2(j)));  %find the column of PM2.5 data from each sampler
                TTmp = table2array(Data2(irow+1:end,icol));
                [~,icol] = max(strcmp(HeaderLine2, 'Actual Volume') & contains(SamplerLine2, Samplers2(j)));  %find the column of PM2.5 data from each sampler
                VTmp = table2array(Data2(irow+1:end,icol));
                if size(TTmp,1)> size(FeTmp,1)
                    diff = size(TTmp,1) - size(FeTmp,1);
                    TTmp(end-diff+1:end) = [];
                    PTmp(end-diff+1:end) = [];
                    VTmp(end-diff+1:end) = [];
                end
                if j ==1
                    SiteData = SiteTmp;
                    DateData = DateTmp;
                    TypeData = TypeTmp;
                    AlData = AlTmp;
                    AlMDLData = AlMDLTmp;
                    AlVData = AlVflagTmp;
                    SiData = SiTmp;
                    SiMDLData = SiMDLTmp;
                    SiVData = SiVflagTmp;
                    SData = STmp;
                    SMDLData = SMDLTmp;
                    SVData = SVflagTmp;
                    CaData = CaTmp;
                    CaMDLData = CaMDLTmp;
                    CaVData = CaVflagTmp;
                    TiData = TiTmp;
                    TiMDLData = TiMDLTmp;
                    TiVData = TiVflagTmp;
                    FeData = FeTmp;
                    FeMDLData = FeMDLTmp;
                    FeVData = FeVflagTmp;
                    TData = TTmp;
                    PData = PTmp;
                    VData = VTmp;
                else
                    SiteData = cat(1,SiteData,SiteTmp);  %concatenate data
                    DateData = cat(1,DateData,DateTmp);
                    TypeData = cat(1,TypeData,TypeTmp);
                    AlData = cat(1,AlData,AlTmp);
                    AlMDLData = cat(1,AlMDLData,AlMDLTmp);
                    AlVData = cat(1,AlVData,AlVflagTmp);
                    SiData = cat(1,SiData,SiTmp);
                    SiMDLData = cat(1,SiMDLData,SiMDLTmp);
                    SiVData = cat(1,SiVData,SiVflagTmp);
                    SData = cat(1,SData,STmp);
                    SMDLData = cat(1,SMDLData,SMDLTmp);
                    SVData = cat(1,SVData,SVflagTmp);
                    CaData = cat(1,CaData,CaTmp);
                    CaMDLData = cat(1,CaMDLData,CaMDLTmp);
                    CaVData = cat(1,CaVData,CaVflagTmp);
                    TiData = cat(1,TiData,TiTmp);
                    TiMDLData = cat(1,TiMDLData,TiMDLTmp);
                    TiVData = cat(1,TiVData,TiVflagTmp);
                    FeData = cat(1,FeData,FeTmp);
                    FeMDLData = cat(1,FeMDLData,FeMDLTmp);
                    FeVData = cat(1,FeVData,FeVflagTmp);
                    TData = cat(1,TData,TTmp);
                    PData = cat(1,PData,PTmp);
                    VData = cat(1,VData,VTmp);
                end
            end
            %filter for routine sampling values that are vaild or historical,
            %and above MDL
            AlData = str2double(AlData);
            AlMDLData = str2double(AlMDLData);
            SiData = str2double(SiData);
            SiMDLData = str2double(SiMDLData);
            SData = str2double(SData);
            SMDLData = str2double(SMDLData);
            CaData = str2double(CaData);
            CaMDLData = str2double(CaMDLData);
            TiData = str2double(TiData);
            TiMDLData = str2double(TiMDLData);
            FeData = str2double(FeData);
            FeMDLData = str2double(FeMDLData);
            TData = str2double(TData);
            PData = str2double(PData);
            VData = str2double(VData);

            AlData(~((contains(AlVData,'H1') | contains(AlVData,'V0')|contains(AlVData,'')) & AlData >= AlMDLData & contains(TypeData,'R'))) = NaN;
            SiData(~((contains(SiVData,'H1') | contains(SiVData,'V0') | contains(SiVData,'')) & SiData >= SiMDLData & contains(TypeData,'R'))) = NaN;
            SData(~((contains(SVData,'H1') | contains(SVData,'V0') | contains(SVData,'')) & SData >= SMDLData & contains(TypeData,'R'))) = NaN;
            CaData(~((contains(CaVData,'H1') | contains(CaVData,'V0') | contains(CaVData,'')) & CaData >= CaMDLData & contains(TypeData,'R'))) = NaN;
            TiData(~((contains(TiVData,'H1') | contains(TiVData,'V0') | contains(TiVData,'')) & TiData >= TiMDLData & contains(TypeData,'R'))) = NaN;
            FeData(~((contains(FeVData,'H1') | contains(FeVData,'V0') | contains(FeVData,'')) & FeData >= FeMDLData & contains(TypeData,'R'))) = NaN;

            %save to final waves that are collecting data from each site
            if xrfcounter==1
                NAPS_XRFData.Site = SiteData;
                NAPS_XRFData.Date = DateData;
                NAPS_XRFData.Al = AlData;
                NAPS_XRFData.Si = SiData;
                NAPS_XRFData.S = SData;
                NAPS_XRFData.Ca = CaData;
                NAPS_XRFData.Ti = TiData;
                NAPS_XRFData.Fe = FeData;
            else
                NAPS_XRFData.Site = cat(1,NAPS_XRFData(:).Site,SiteData);
                NAPS_XRFData.Date = cat(1,NAPS_XRFData(:).Date,DateData);
                NAPS_XRFData.Al = cat(1,NAPS_XRFData(:).Al,AlData);
                NAPS_XRFData.Si = cat(1,NAPS_XRFData(:).Si,SiData);
                NAPS_XRFData.S = cat(1,NAPS_XRFData(:).S,SData);
                NAPS_XRFData.Ca = cat(1,NAPS_XRFData(:).Ca,CaData);
                NAPS_XRFData.Ti = cat(1,NAPS_XRFData(:).Ti,TiData);
                NAPS_XRFData.Fe = cat(1,NAPS_XRFData(:).Fe,FeData);
            end
       xrfcounter=xrfcounter+1;
       end
       
       %IC Data
       if LoadIC ==1
            clear Data
            Data = readtable(filename,'Sheet','Ions-Spec_IC','ReadVariableNames',false);%load data from XRF sheet 
            [~,irow] = max(contains(table2cell(Data(:,1)),'Sampler'));                  %find the sampler header line
            SamplerLine = table2cell(Data(irow,:));                                     %record all the sampler headers
            SamplerLine(cellfun(@isempty,SamplerLine)) = {'NaN'};                       % convert all non-character elements into characters
            Samplers = SamplerLine(~contains(SamplerLine,'NaN') & ~contains(SamplerLine,'Sampler')); %find all sample types
            Samplers = unique(Samplers,'stable');
            numSamplers = size(Samplers,2);                                             %get the number of samplers
            [~,irow] = max(contains(table2cell(Data(:,1)),'NAPS Site ID'));             %find the data header line
            HeaderLine = table2cell(Data(irow,:));                                      %record all the data headers
            HeaderLine(cellfun(@isempty,HeaderLine)) = {'NaN'};
            
            Data3 = readtable(filename,'Sheet','PM2.5','ReadVariableNames',false);      %load data from XRF sheet 
            [~,irow] = max(contains(table2cell(Data3(:,1)),'Sampler'));                 %find the sampler header line
            SamplerLine3 = table2cell(Data3(irow,:));                                   %record all the sampler headers
            SamplerLine3(cellfun(@isempty,SamplerLine3)) = {'NaN'};                     % convert all non-character elements into characters
            Samplers3 = SamplerLine3(~contains(SamplerLine3,'NaN') & ~contains(SamplerLine3,'Sampler')); %find all sample types
            Samplers3 = unique(Samplers3,'stable');
            [~,irow] = max(contains(table2cell(Data3(:,1)),'NAPS Site ID'));            %find the data header line
            HeaderLine3 = table2cell(Data3(irow,:));                                    %record all the data headers
            HeaderLine3(cellfun(@isempty,HeaderLine3)) = {'NaN'};
            
            LoadvNO3 = max(strcmp(table2cell(MetaData(:,2)),'Volatile Nitrate_IC'));
            if LoadvNO3 ==1
                Data2 = readtable(filename,'Sheet','Volatile Nitrate_IC','ReadVariableNames',false); %load data from volatile nitrate sheet 
                [~,irow] = max(contains(table2cell(Data2(:,1)),'Sampler'));             %find the sampler header line
                SamplerLine2 = table2cell(Data2(irow,:));                               %record all the sampler headers
                SamplerLine2(cellfun(@isempty,SamplerLine2)) = {'NaN'};                 % convert all non-character elements into characters
                Samplers2 = SamplerLine2(~contains(SamplerLine2,'NaN') & ~contains(SamplerLine2,'Sampler')); %find all sample types
                Samplers2 = unique(Samplers2,'stable');
                [~,irow] = max(contains(table2cell(Data2(:,1)),'NAPS Site ID'));        %find the data header line
                HeaderLine2 = table2cell(Data2(irow,:));                                %record all the data headers
                HeaderLine2(cellfun(@isempty,HeaderLine2)) = {'NaN'};
            end
            for j=1:numSamplers
                [~,icol] = max(contains(HeaderLine,'NAPS Site ID'));
                SiteTmp = table2array(Data(irow+1:end,icol));                           %extract all the site code data
                [~,icol] = max(contains(HeaderLine,'Sampling Date'));
                DateTmp = table2cell(Data(irow+1:end,icol));                            %extract all sample times
                [~,icol] = max(contains(HeaderLine,'Sampling Type'));
                TypeTmp = table2cell(Data(irow+1:end,icol));                            % extract all sample types
                [~,icol] = max(strcmp(HeaderLine, 'Chloride') & contains(SamplerLine, Samplers(j)));  %find the column of Chloride data from each sampler
                ChlTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'Chloride-MDL') & contains(SamplerLine, Samplers(j)));
                ChlMDLTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'Chloride-VFlag') & contains(SamplerLine, Samplers(j))); 
                ChlVflagTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(strcmp(HeaderLine, 'Sulphate') & contains(SamplerLine, Samplers(j)));  %find the column of sulfate data from each sampler
                SO4Tmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'Sulphate-MDL') & contains(SamplerLine, Samplers(j)));  
                SO4MDLTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'Sulphate-VFlag') & contains(SamplerLine, Samplers(j)));  
                SO4VflagTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(strcmp(HeaderLine, 'Nitrate') & contains(SamplerLine, Samplers(j)));     %find the column of nitrate data from each sampler
                NO3Tmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'Nitrate-MDL') & contains(SamplerLine, Samplers(j)));  
                NO3MDLTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'Nitrate-VFlag') & contains(SamplerLine, Samplers(j)));  
                NO3VflagTmp = table2array(Data(irow+1:end,icol));
                if LoadvNO3 ==1 %combine volatile and non-volatile nitrates
                    [~,icol] = max(strcmp(HeaderLine2, 'Nitrate') & contains(SamplerLine2, Samplers2(j))); 
                    vNO3Tmp = table2array(Data2(irow+1:end,icol));
                    [~,icol] = max(strcmp(HeaderLine2, 'Nitrite') & contains(SamplerLine2, Samplers2(j)));  
                    vNO2Tmp = table2array(Data2(irow+1:end,icol));
                    vNO3Tmp = str2double(vNO3Tmp)+str2double(vNO2Tmp);
                    NO3Tmp = str2double(NO3Tmp);
                else
                    vNO3Tmp = NaN(size(NO3Tmp,1),1);
                    NO3Tmp = NaN(size(NO3Tmp,1),1);
                end
                [~,icol] = max(strcmp(HeaderLine, 'Ammonium') & contains(SamplerLine, Samplers(j)));  %find the column of NH4 data from each sampler
                NH4Tmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'Ammonium-MDL') & contains(SamplerLine, Samplers(j)));  
                NH4MDLTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'Ammonium-VFlag') & contains(SamplerLine, Samplers(j)));  
                NH4VflagTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(strcmp(HeaderLine3, 'Pres.') & contains(SamplerLine3, Samplers3(j)));  %find the column of pressure data from each sampler
                PTmp = table2array(Data3(irow+1:end,icol));
                [~,icol] = max(strcmp(HeaderLine3, 'Temp.') & contains(SamplerLine3, Samplers3(j)));  %find the column of temperature data from each sampler
                TTmp = table2array(Data3(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine3, 'Actual Volume') & contains(SamplerLine3, Samplers3(j)));  %find the column of volume data from each sampler
                VTmp = table2array(Data3(irow+1:end,icol));
                if size(TTmp,1)> size(NH4Tmp,1)
                    diff = size(TTmp,1) - size(NH4Tmp,1);   %resize the met data to match the number of data points
                    TTmp(end-diff+1:end) = [];
                    PTmp(end-diff+1:end) = [];
                    VTmp(end-diff+1:end) = [];
                end
                if j ==1
                    SiteData = SiteTmp;
                    DateData = DateTmp;
                    TypeData = TypeTmp;
                    ChlData = ChlTmp;
                    ChlMDLData = ChlMDLTmp;
                    ChlVData = ChlVflagTmp;
                    SO4Data = SO4Tmp;
                    SO4MDLData = SO4MDLTmp;
                    SO4VData = SO4VflagTmp;
                    NO3Data = NO3Tmp;
                    NO3MDLData = NO3MDLTmp;
                    NO3VData = NO3VflagTmp;
                    vNO3Data = vNO3Tmp;
                    NH4Data = NH4Tmp;
                    NH4MDLData = NH4MDLTmp;
                    NH4VData = NH4VflagTmp;
                    PData = PTmp;
                    TData = TTmp;
                    VData = VTmp;
                else
                    SiteData = cat(1,SiteData,SiteTmp);  %concatenate data
                    DateData = cat(1,DateData,DateTmp);
                    TypeData = cat(1,TypeData,TypeTmp);
                    ChlData = cat(1,ChlData,ChlTmp);
                    ChlMDLData = cat(1,ChlMDLData,ChlMDLTmp);
                    ChlVData = cat(1,ChlVData,ChlVflagTmp);
                    SO4Data = cat(1,SO4Data,SO4Tmp);
                    SO4MDLData = cat(1,SO4MDLData,SO4MDLTmp);
                    SO4VData = cat(1,SO4VData,SO4VflagTmp);
                    NO3Data = cat(1,NO3Data,NO3Tmp);
                    NO3MDLData = cat(1,NO3MDLData,NO3MDLTmp);
                    NO3VData = cat(1,NO3VData,NO3VflagTmp);
                    vNO3Data = cat(1,vNO3Data,vNO3Tmp);
                    NH4Data = cat(1,NH4Data,NH4Tmp);
                    NH4MDLData = cat(1,NH4MDLData,NH4MDLTmp);
                    NH4VData = cat(1,NH4VData,NH4VflagTmp);
                    PData = cat(1,PData,PTmp);
                    TData = cat(1,TData,TTmp);
                    VData = cat(1,VData,VTmp);
                end
            end
            %filter for routine sampling values that are vaild or historical,
            %and above MDL
            ChlData = str2double(ChlData);
            ChlMDLData = str2double(ChlMDLData);
            SO4Data = str2double(SO4Data);
            SO4MDLData = str2double(SO4MDLData);
            NO3Data = (NO3Data)+(vNO3Data); %add volatile and particulate nitrate
            NO3MDLData = str2double(NO3MDLData);
            NH4Data = str2double(NH4Data);
            NH4MDLData = str2double(NH4MDLData);
            PData = str2double(PData);
            TData = str2double(TData);
            VData = str2double(VData);
 
            ChlData(~((contains(ChlVData,'H1') |contains(ChlVData,'V0') | contains(ChlVData,'')) & ChlData >= ChlMDLData & contains(TypeData,'R'))) = NaN;
            SO4Data(~((contains(SO4VData,'H1') |contains(SO4VData,'V0') | contains(SO4VData,'')) & SO4Data >= SO4MDLData & contains(TypeData,'R'))) = NaN;
            NO3Data(~((contains(NO3VData,'H1') |contains(NO3VData,'V0') | contains(NO3VData,'')) & NO3Data >= NO3MDLData & contains(TypeData,'R'))) = NaN;
            NH4Data(~((contains(NH4VData,'H1') |contains(NH4VData,'V0') | contains(NH4VData,'')) & NH4Data >= NH4MDLData & contains(TypeData,'R'))) = NaN;

            %save to final waves that are collecting data from each site
            if iccounter==1
                NAPS_ICData.Site = SiteData;
                NAPS_ICData.Date = DateData;
                NAPS_ICData.Chl = ChlData;
                NAPS_ICData.SO4 = SO4Data;
                NAPS_ICData.NO3 = NO3Data;
                NAPS_ICData.NH4 = NH4Data;
            else
                NAPS_ICData.Site = cat(1,NAPS_ICData(:).Site,SiteData);
                NAPS_ICData.Date = cat(1,NAPS_ICData(:).Date,DateData);
                NAPS_ICData.Chl = cat(1,NAPS_ICData(:).Chl,ChlData);
                NAPS_ICData.SO4 = cat(1,NAPS_ICData(:).SO4,SO4Data);
                NAPS_ICData.NO3 = cat(1,NAPS_ICData(:).NO3,NO3Data);
                NAPS_ICData.NH4 = cat(1,NAPS_ICData(:).NH4,NH4Data);
            end
            iccounter = iccounter+1;
       end
       
       %Organic Data
       if LoadOrganics ==1
            clear Data
            Data = readtable(filename,'Sheet','OCEC & OC Artifacts','ReadVariableNames',false); %load data from XRF sheet 
            [~,irow] = max(contains(table2cell(Data(:,1)),'Cartridge'));                        %find the sampler header line
            CartridgeLine = table2cell(Data(irow,:));                                           %record all the sampler headers
            CartridgeLine(cellfun(@isempty,CartridgeLine)) = {'NaN'};                           % convert all non-character elements into characters
            Cartridges = CartridgeLine(~contains(CartridgeLine,'NaN') & ~contains(CartridgeLine,'Cartridge')); %find all cartridge types
            Data2 = readtable(filename,'Sheet','PM2.5','ReadVariableNames',false);              %load data from XRF sheet 
            [~,irow] = max(contains(table2cell(Data2(:,1)),'NAPS Site ID'));                    %find the data header line
            HeaderLine2 = table2cell(Data2(irow,:));                                            %record all the data headers
            HeaderLine2(cellfun(@isempty,HeaderLine2)) = {'NaN'};
            if max(contains(Cartridges,'A,B')) ==1
                [~,irow] = max(contains(table2cell(Data(:,1)),'NAPS Site ID'));                 %find the data header line
                HeaderLine = table2cell(Data(irow,:));                                          %record all the data headers
                HeaderLine(cellfun(@isempty,HeaderLine)) = {'NaN'};  
                [~,icol] = max(contains(HeaderLine,'NAPS Site ID'));
                SiteTmp = table2array(Data(irow+1:end,icol));                                   %extract all the site code data
                [~,icol] = max(contains(HeaderLine,'Sampling Date'));
                DateTmp = table2cell(Data(irow+1:end,icol));                                    %extract all sample times
                [~,icol] = max(contains(HeaderLine,'Sampling Type'));
                TypeTmp = table2cell(Data(irow+1:end,icol));                                    % extract all sample types
                [~,icol] = max(strcmp(HeaderLine, 'OC(corr)'));                                 %find the column of OC data from each sampler
                OCTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'OC(corr)-MDL'));  
                OCMDLTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'OC(corr)-VFlag'));  
                OCVflagTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(strcmp(HeaderLine, 'TC(corr)'));                                 %find the column of total carbon data from each sampler
                TCTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'TC(corr)-MDL'));  
                TCMDLTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine, 'TC(corr)-VFlag'));  
                TCVflagTmp = table2array(Data(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine2, 'Pres.'));                                 %find the column of pressure data from each sampler
                PTmp = table2array(Data2(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine2, 'Temp.'));                                 %find the column of temperature data from each sampler
                TTmp = table2array(Data2(irow+1:end,icol));
                [~,icol] = max(contains(HeaderLine2, 'Actual Volume'));                         %find the column of volume data from each sampler
                VTmp = table2array(Data2(irow+1:end,icol));
                 if size(TTmp,1)> size(OCTmp,1)
                    diff = size(TTmp,1) - size(OCTmp,1);                                        %resize met data to correct number of points
                    TTmp(end-diff+1:end) = [];
                    PTmp(end-diff+1:end) = [];
                    VTmp(end-diff+1:end) = [];
                end
                SiteData = SiteTmp;
                DateData = DateTmp;
                TypeData = TypeTmp;
                OCData = str2double(OCTmp);
                OCMDLData = str2double(OCMDLTmp);
                TCData = str2double(TCTmp);
                TCMDLData = str2double(TCMDLTmp);
                PData = str2double(PTmp);
                VData = str2double(VTmp);
                TData = str2double(TTmp);

                OCData(~((contains(OCVflagTmp,'H1') |contains(OCVflagTmp,'V0') | contains(OCVflagTmp,'')) & OCData > OCMDLData & contains(TypeData,'R'))) = NaN;
                TCData(~((contains(TCVflagTmp,'H1') |contains(TCVflagTmp,'V0') | contains(TCVflagTmp,'')) & TCData > TCMDLData & contains(TypeData,'R'))) = NaN;
                ECData = TCData-OCData; %calculate the total corrected EC from OC and TC
                
                 %save to final waves that are collecting data from each site
                if orgcounter ==1
                    NAPS_CData.Site = SiteData;
                    NAPS_CData.Date = DateData;
                    NAPS_CData.OC = OCData;
                    NAPS_CData.EC = ECData;
                else
                    NAPS_CData.Site = cat(1,NAPS_CData(:).Site,SiteData);
                    NAPS_CData.Date = cat(1,NAPS_CData(:).Date,DateData);
                    NAPS_CData.OC = cat(1,NAPS_CData(:).OC,OCData);
                    NAPS_CData.EC = cat(1,NAPS_CData(:).EC,ECData);
                end
            else
                disp("**Organic Error: No corrected Data**")
            end 
       orgcounter = orgcounter+1;
       end
   end

    [~,~,allStations] = xlsread(sprintf('%sStations2017_v3.xlsx',DataFileLoc),'Stations_13DEC2017');
    allSiteID = uint32(cell2mat(allStations(2:end,1)));
    %allSiteLat = single(cell2mat(allStations(2:end,9)));
    %allSiteLon = single(cell2mat(allStations(2:end,10)));
    %allSiteElev = single(cell2mat(allStations(2:end,11)));
    allSiteType = allStations(2:end,52);
    allSiteLU = allStations(2:end,54);
    
    %%%Load Hourly files (ug / m3) - volumetric. Transformations are uncertain %%%%%%
    if Hourly ==1
        NAPS_HLY.Site = [];
        NAPS_HLY.Date = [];
        NAPS_HLY.Inst = [];
        NAPS_HLY.Lat = [];
        NAPS_HLY.Lon = [];
        NAPS_HLY.Elev = [];
        NAPS_HLY.PM = [];

        HrlyFiles = dir(sprintf('%s/HLY/PM25_%d.csv',DataFileLoc,year));
        for i = 1:numel(HrlyFiles)
            Data = readtable(HrlyFiles(i).name,'ReadVariableNames',false); %load data from csv file 
            [~,irow] = max(contains(table2cell(Data(:,1)),'PM2.5'));           %find the first line of data
            InstIDTmp = table2cell(Data(irow:end,2));  
            SiteIDTmp = table2cell(Data(irow:end,3)); 
            SiteLatTmp = table2cell(Data(irow:end,6));
            SiteLonTmp = table2cell(Data(irow:end,7));
            DateTmp = table2cell(Data(irow:end,8));  
            Tmp = str2double(table2array(Data(irow:end,9:end)));  
            Tmp(Tmp<0)=NaN;
            nummeas = sum(~isnan(Tmp(:,:)),2);      %number of measurements each day
            PMTmp= nanmean(Tmp(:,:),2);              %daily averages
            PMTmp(nummeas < 24)=NaN;                %report daily averages with at least 20 points
            NAPS_HLY.Site = cat(1,NAPS_HLY(:).Site,SiteIDTmp);
            NAPS_HLY.Date = cat(1,NAPS_HLY(:).Date,DateTmp);
            NAPS_HLY.Inst = cat(1,NAPS_HLY(:).Inst,InstIDTmp);
            NAPS_HLY.Lat = cat(1,NAPS_HLY(:).Lat,SiteLatTmp);
            NAPS_HLY.Lon = cat(1,NAPS_HLY(:).Lon,SiteLonTmp);
            NAPS_HLY.PM = cat(1,NAPS_HLY(:).PM,PMTmp); 
            notmember= ~ismember(str2double(NAPS_HLY.Site),NAPS_MonAvg_SiteInfo(:,4));   %set to 1 if site is not already in site list
            newSiteLat = str2double(unique(SiteLatTmp(notmember==1),'stable'));
            newSiteLon = str2double(unique(SiteLonTmp(notmember==1),'stable'));
            newSiteIDTmp = str2double(unique(SiteIDTmp(notmember==1),'stable'));
            for isite=1:size(newSiteLat,1)
                j = size(NAPS_MonAvg_SiteInfo,1)+1;                                    %add new sites to site list
                NAPS_MonAvg_SiteInfo(j,1) = double(newSiteLat(isite)); 
                NAPS_MonAvg_SiteInfo(j,2) = double(newSiteLon(isite));
                NAPS_MonAvg_SiteInfo(j,3) = NaN;%double(ElevTmp);
                NAPS_MonAvg_SiteInfo(j,4) = double(newSiteIDTmp(isite));
            end
            fprintf('%s: added.\n',HrlyFiles(i).name);               
        end
        % See correspondence from Randall from Aaron from EC on Sept 29, 2014: "Fwd: Re: TOEM bias equations"
        if (ApplyTransformations == 1)
            coldtmp = string(num2cell(NAPS_HLY.Date));
            coldtmp = str2double(extractBetween(coldtmp,5,6)); %extract date
            coldspot = (coldtmp <=3 | coldtmp >=10);
            warmspot = (coldtmp <=9 & coldtmp >=4);
                TEOMspot=(strcmp(NAPS_HLY.Inst,'706') | strcmp(NAPS_HLY.Inst,'760')); %TEOM or TEOM+FDMS
                SHARPspot=(strcmp(NAPS_HLY.Inst,'184'));
                    NAPS_HLY.PM(TEOMspot&coldspot) = NAPS_HLY.PM(TEOMspot&coldspot) .*1.3 + 0.94;
                    NAPS_HLY.PM(TEOMspot&warmspot) = NAPS_HLY.PM(TEOMspot&warmspot) .*0.94 + 1.72;
                    NAPS_HLY.PM(SHARPspot&coldspot) = NAPS_HLY.PM(SHARPspot&coldspot) .*0.95 - 0.33;
                    NAPS_HLY.PM(SHARPspot&warmspot) = NAPS_HLY.PM(SHARPspot&warmspot) .*1.01 -0.22; 
        end
   
        %add to PM Data structure
        NAPS_PMData.Site = cat(1,NAPS_PMData(:).Site,string(num2cell(NAPS_HLY.Site)));
        NAPS_PMData.Date = cat(1,NAPS_PMData(:).Date, string(datetime(string(num2cell(NAPS_HLY.Date)),'InputFormat','yyyyMMdd')));
        NAPS_PMData.PM = cat(1,NAPS_PMData(:).PM,NAPS_HLY.PM);
    end
   
    %filter to remove industrial sites and point source-influenced sites
    numsites = size(allStations,1);
    for i=1:numsites-1
       SiteIDTmp = string(num2cell(allSiteID(i)));
       SiteLUTmp = allSiteLU(i);
       SiteTypeTmp = allSiteType(i);
       NAPS_PMData.PM(strcmp(NAPS_PMData.Site,SiteIDTmp) & (strcmp(SiteLUTmp,'I') | strcmp(SiteTypeTmp, 'PS'))) = NaN;
       NAPS_ICData.NO3(strcmp(NAPS_PMData.Site,SiteIDTmp) & (strcmp(SiteLUTmp,'I') | strcmp(SiteTypeTmp, 'PS'))) = NaN;
       NAPS_ICData.NH4(strcmp(NAPS_ICData.Site,SiteIDTmp) & (strcmp(SiteLUTmp,'I') | strcmp(SiteTypeTmp, 'PS'))) = NaN;
       NAPS_ICData.SO4(strcmp(NAPS_ICData.Site,SiteIDTmp) & (strcmp(SiteLUTmp,'I') | strcmp(SiteTypeTmp, 'PS'))) = NaN;
       NAPS_ICData.Chl(strcmp(NAPS_ICData.Site,SiteIDTmp) & (strcmp(SiteLUTmp,'I') | strcmp(SiteTypeTmp, 'PS'))) = NaN;
       NAPS_XRFData.Al(strcmp(NAPS_XRFData.Site,SiteIDTmp) & (strcmp(SiteLUTmp,'I') | strcmp(SiteTypeTmp, 'PS'))) = NaN;
       NAPS_XRFData.Si(strcmp(NAPS_XRFData.Site,SiteIDTmp) & (strcmp(SiteLUTmp,'I') | strcmp(SiteTypeTmp, 'PS'))) = NaN;
       NAPS_XRFData.S(strcmp(NAPS_XRFData.Site,SiteIDTmp) & (strcmp(SiteLUTmp,'I') | strcmp(SiteTypeTmp, 'PS'))) = NaN;
       NAPS_XRFData.Ca(strcmp(NAPS_XRFData.Site,SiteIDTmp) & (strcmp(SiteLUTmp,'I') | strcmp(SiteTypeTmp, 'PS'))) = NaN;
       NAPS_XRFData.Ti(strcmp(NAPS_XRFData.Site,SiteIDTmp) & (strcmp(SiteLUTmp,'I') | strcmp(SiteTypeTmp, 'PS'))) = NaN;
       NAPS_XRFData.Fe(strcmp(NAPS_XRFData.Site,SiteIDTmp) & (strcmp(SiteLUTmp,'I') | strcmp(SiteTypeTmp, 'PS'))) = NaN;
       NAPS_CData.OC(strcmp(NAPS_CData.Site,SiteIDTmp) & (strcmp(SiteLUTmp,'I') | strcmp(SiteTypeTmp, 'PS'))) = NaN;
       NAPS_CData.EC(strcmp(NAPS_CData.Site,SiteIDTmp) & (strcmp(SiteLUTmp,'I') | strcmp(SiteTypeTmp, 'PS'))) = NaN;
   end

   warning('on','all')
   
   %Make monthly averages for each site (variable avg: site x month x variable )
   MonStr = {strcat('Jan-',YearStr), strcat('Feb-',YearStr), strcat('Mar-',YearStr), strcat('Apr-',YearStr), ...
       strcat('May-',YearStr), strcat('Jun-',YearStr), strcat('Jul-',YearStr), strcat('Aug-',YearStr), ...
       strcat('Sep-',YearStr), strcat('Oct-',YearStr), strcat('Nov-',YearStr), strcat('Dec-',YearStr)};
   MonStr = cellstr(MonStr)';
   SiteStr = NAPS_MonAvg_SiteInfo(:,4);
   SiteStr = string(SiteStr);
   WaveStr = {NAPS_PMData.PM,NAPS_ICData.NO3,NAPS_ICData.NH4,NAPS_ICData.SO4,NAPS_ICData.Chl,...
       NAPS_XRFData.Al,NAPS_XRFData.Si,NAPS_XRFData.S,NAPS_XRFData.Ca,NAPS_XRFData.Ti,NAPS_XRFData.Fe,...
       NAPS_CData.OC,NAPS_CData.EC};
   VarStr = {'PM','NO3','NH4','SO4','Chl','Al','Si','S','Ca','Ti','Fe','OC','EC'};
   numVar = size(VarStr,2);

   numsites = size(SiteStr,1);
 
  %calc PM averages
  NAPS_MonAvg_ugsm3 = NaN(numsites,12,numVar);
  for i=1:numsites
    for j=1:12 %for each month, filter through the data from each instrument
        for k=1:numVar
        %filter for values for this site, at this month, for each
        %desired variable and set the MonAvg array accordingly
            switch k
                case 1
                    clear filtertmp DataTmp SiteTmp DateTmp
                    SiteTmp = NAPS_PMData.Site;
                    DateTmp = NAPS_PMData.Date;
                case num2cell(2:5)
                    clear filtertmp DataTmp SiteTmp DateTmp
                    SiteTmp = NAPS_ICData.Site;
                    DateTmp = NAPS_ICData.Date;
                case num2cell(6:11)
                    clear filtertmp DataTmp SiteTmp DateTmp
                    SiteTmp = NAPS_XRFData.Site;
                    DateTmp = NAPS_XRFData.Date;
                case num2cell(12:13)
                    clear filtertmp DataTmp SiteTmp DateTmp
                    SiteTmp = NAPS_CData.Site;
                    DateTmp = NAPS_CData.Date;
            end
            filtertmp = (contains(SiteTmp, SiteStr(i)) & contains(DateTmp, MonStr(j)));
            DataTmp = WaveStr{k};
            tmp = DataTmp(filtertmp);
            tmp(tmp<0) = NaN;
            if contains(MonStr(j),'Feb') 
                dayfilter2=4;
                else;dayfilter2=5;
            end
            if sum(~isnan(tmp),1) >= dayfilter2                              %only record data if >= 4 measurements per month
                NAPS_MonAvg_ugsm3(i,j,k) = nanmean(tmp);   %take average of all month data from each site for each variable
            else
                NAPS_MonAvg_ugsm3(i,j,k) = NaN;
            end
        end
    end
  end
    savefile = sprintf('/misc/data6/emcduffie/ExtData/NAPS_Data/NAPS_MonthlyAvgs_%s.mat',YearStr);
    save(savefile,'NAPS_MonAvg_ugsm3','NAPS_MonAvg_SiteInfo');
 end
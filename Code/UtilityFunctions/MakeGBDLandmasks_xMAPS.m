%%%% MakeGBDLandmasks.m
%%%% E.E. McDuffie, March 2019
%%%
%%% Take Country Border Shapefile and Make GBD Region LandMasks at choice resolution
%%% Input: Country boarder shapefile specified by S
%%% Output: sGBDCountries - structure with 204 country land masks
%%%         sGBDRegions - structure with 21 region land masks
%%% Takes a long time (~> week) to calculate 1x1 km masks!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% USER SETTINGS
resolution = [0.01 0.01];   %desired resolution (lat, lon degree) of the country masks
savedata =1;                %don't save data = 0, save new data = 1 (if res  = 0.01, data will always save)
plotdata=0;                 %plot data (0=no, 1=yes)? (only applies to res > 0.01)
calc_country=0;             %calc the country data (1=yes, 0 = no) or just the region data. WARNING: 0.01 res takes ~1-2 weeks for countries
outfile = '/misc/data6/emcduffie/Matlab_Scripts/GBD-MAPS_Scripts_PublicScriptsPackage/Masks/'; %base location where to save output data
%%%%
% COUNTRY BORDER SHAPEFILES (user needs to add own files)
S = shaperead('/path/to/file/Countries_WGS84.shp','UseGeoCoords', true);% from https://hub.arcgis.com/datasets/a21fdb46d23e4ef896f31475217cbb08_1)
S2 = shaperead('/path/to/file/GBD2020_mapping_final.shp','UseGeoCoords',true);
FileName = '/path/to/file/IHME_GBD_2019_CountriesList_Y2020M07D10.xlsx';% Load list of all GBD countries and Regions
%%%%%

%%% Step 0. Make list of GBD country and region names and set mask resolution
for doloop=1:1
Data = readtable(FileName);
SetID = cell2mat(table2cell(Data(:,1)));    %first column is data source (GBD regions = 423)
LocID = cell2mat(table2cell(Data(:,2)));    %location ID
LocName = table2cell(Data(:,3));            %location name
ParentID = cell2mat(table2cell(Data(:,4))); %Parent location ID
Level = cell2mat(table2cell(Data(:,5)));    %Location Level (1 = global, 2 = GBD Region (21), 3 = Country/territory (195))
LocID =LocID(SetID ==423);                  %filter list for GBD reported region
LocName = LocName(SetID==423);
ParentID = ParentID(SetID==423);
Level = Level(SetID ==423);
SetID = SetID(SetID==423);
GBD_Regions = LocName(Level==2);            % make list of GBD region names
GBD_Countries = LocName(Level==3);          % make list of GBD country names
ParentID=ParentID(Level==3);
numGBDc = size(GBD_Countries,1);
numGBDr = size(GBD_Regions,1);
GBD_Countries_Rlist = {'NaN'};
for icountry=1:numGBDc
    GBD_Countries_Rlist(icountry,1) = LocName(ParentID(icountry)==LocID);
end
%%%%

%calculate Lat, Lon grids
res = [resolution(1),resolution(2)];
numlat = 180/res(1);
numlon = 360/res(2);
Lat = (-90+(res(1)/2):res(1):90-(res(1)/2))';
Lon = (-180+(res(2)/2):res(2):180-(res(2)/2))';
[fLAT,fLON]=meshgrid(Lat,Lon);

% extract list of countries from the shapefiles
Countries = extractfield(S,'CNTRY_NAME');
Countries2 = extractfield(S2,'loc_name');
%change country names for this particular shapefile to match GBD country names
Countries2(strcmp(Countries2,'Russia'))={'Russian Federation'};
Countries2(strcmp(Countries2,'Gambia'))={'The Gambia'};
val = cell(numGBDc,1);

%%% Make Structure sGBDCountries that contain Name, GBD Region, Lat, Lon
%%% (border), and mask for all GBD countries
sGBDCountries = struct('Name',val,'Region',val,'Lat',val,'Lon',val,'Mask',val);
end

%%% Step 1. Calculate the country masks
if calc_country==1
 for icountry=1:numGBDc
     sGBDCountries(icountry).Name = GBD_Countries{icountry};
     sGBDCountries(icountry).Region = GBD_Countries_Rlist{icountry};
     if strcmp(GBD_Countries(icountry),'Monaco') || strcmp(GBD_Countries(icountry),'Nauru') ||...
        strcmp(GBD_Countries(icountry),'Marshall Islands') || strcmp(GBD_Countries(icountry),'Tokelau') ||...
        strcmp(GBD_Countries(icountry),'Tuvalu') || strcmp(GBD_Countries(icountry),'Maldives')  % not in GBD shapefiles, so take from ArcGIS shapefile
        sGBDCountries(icountry).Lat = S(strcmp(Countries,GBD_Countries(icountry))).Lat';
        sGBDCountries(icountry).Lon = S(strcmp(Countries,GBD_Countries(icountry))).Lon';
    else
        match = find(strcmp(Countries2,GBD_Countries(icountry)));                               %Georgia and Niger have multiple entries
        if strcmp(Countries2(match(1)),'Niger')
            match_num=2;
        else
            match_num=1;
        end
        sGBDCountries(icountry).Lat = S2(match(match_num)).Lat';                                %set lat and lon based on which occurance number
        sGBDCountries(icountry).Lon = S2(match(match_num)).Lon';
    end
    mask = int8(inpolygon(fLON,fLAT,sGBDCountries(icountry).Lon,sGBDCountries(icountry).Lat));  %use the Lat and Lon boundaries from the shapefiles to calculated gridded masks
    sGBDCountries(icountry).Mask = mask; 
    %for high-res data, not enough memory to store data in single structure, so save country masks as individual files in the specified location
    if res(1)==0.01
        FileName = sprintf('%sGBD_Country_Masks_%s_0.01.mat',outfile,GBD_Countries{icountry});  
        sCountryMask=sGBDCountries(icountry);
        save(FileName,'sCountryMask','-v7.3');
        sGBDCountries = struct('Name',val,'Region',val,'Lat',val,'Lon',val,'Mask',val);         %in this case, empty make placeholder
    end
     fprintf('Finished Mask: %s (%d)\n',GBD_Countries{icountry}, icountry)
 end
end

%Step 2. combine country masks to get masks for 21 regions
val = cell(numGBDr,1);
sGBDRegions = struct('name',val,'Lat',val,'Lon',val,'Mask',val);
if res(1)==0.01
    FileName = sprintf('%sGBD_Country_Masks_0.10.mat', outfile);   %if the high resolution, need to reload the low res sGBDCountries matrix (just for country names)
    load(FileName);
end
RegionTmp = extractfield(sGBDCountries,'Region');
CountryNames = extractfield(sGBDCountries,'Name');
for iregion=1:numGBDr
    sGBDRegions(iregion).name = GBD_Regions{iregion};
    CountryList = (strcmp(RegionTmp,GBD_Regions{iregion})); 
    numCountry = sum(CountryList);
    CountryIndex = find(CountryList==1); %find index of all the countries of the given region
    Lattmp = [];
    Lontmp=[];
    Masktmp=repmat(0,numlon,numlat);
    for icountry = 1:numCountry
        if res(1)==0.01
            nametmp = CountryNames{CountryIndex(icountry)};
            nametmp = strrep(nametmp," ","_");          %replace spaces with underscores
            nametmp = strrep(nametmp,"-","_");          %replace dashes with underscores
            nametmp = strrep(nametmp,"'","");           %replace apostrophes with nothing (Cote d'Iviore)
            nametmp = strrep(nametmp,",","");           %replace commas with nothing (Virgin Islands, U.S.)
            nametmp = strrep(nametmp,".","");           %replace periods with nothing (Virgin Islands, U.S.)
            FileName = sprintf('%sGBD_Country_Masks_%s_0.01.mat',outfile,nametmp);
            load(FileName);
            Masktmp(sCountryMask.Mask>0)=1;
            Lattmp = cat(1,Lattmp,sCountryMask.Lat);
            Lontmp=cat(1,Lontmp,sCountryMask.Lon);
        else
            Lattmp = cat(1,Lattmp,sGBDCountries(CountryIndex(icountry)).Lat);
            Lontmp = cat(1,Lontmp,sGBDCountries(CountryIndex(icountry)).Lon);
            Masktmp(sGBDCountries(CountryIndex(icountry)).Mask>0)=1;    %set the mask to 1 for each country
        end
    end
    sGBDRegions(iregion).Lat =Lattmp;
    sGBDRegions(iregion).Lon =Lontmp;
    sGBDRegions(iregion).Mask =Masktmp;
    % if high resolution, not enough memory to store all masks in one
    % structure, so saved as separate files. **Need to rename to replace
    % spaces with underscores for the GBD_MAPS_Global_Analysis_Main Script**
    if res(1)==0.01
        FileName = sprintf('%sGBD_Region_Masks_%s_0.01.mat',outfile,GBD_Regions{iregion});
        sRegionMask=sGBDRegions(iregion);
        save(FileName,'sRegionMask','-v7.3');
        sGBDRegions = struct('name',val,'Lat',val,'Lon',val,'Mask',val);
    end
    fprintf('Finished Mask: %s\n',GBD_Regions{iregion})
end

%%% Step 3. save data
gLAT=Lat;
gLON=Lon;
FileName = sprintf('%sGBD_Country_Masks_%.2f.mat',outfile,res(1));
FileName2 = sprintf('%sGBD_Region_Masks_%.2f.mat',outfile,res(1));
if savedata==1 && res(1) > 0.01
    save(FileName,'sGBDCountries','gLAT','gLON','-v7.3');
    fprintf('Saving data to: %s\n',FileName);
    save(FileName2,'sGBDRegions','gLAT','gLON','-v7.3');
    fprintf('Saving data to: %s\n',FileName2);
end

%%% Step 4. Plot region masks and borders
for doloop=1:1
if plotdata==1 && res(1) > 0.01
    figure
    load coast;
    worldmap('world');
    set(gcf,'color','white');
    setm(gca,'ffacecolor',[0.9 0.9 0.9])
    setm(gca,'Grid','off','MapProjection','miller','parallellabel','off','meridianlabel','off','MapLatLimit',[-65 85],'MapLonLimit',[-179 179]) 

    for iregion=1:numGBDr
        tmp = sGBDRegions(iregion).Mask;    
        tmp(tmp<1)=NaN;
        tmp(~isnan(tmp)) = iregion;
        c(iregion) = surfm(fLAT,fLON,tmp);
    end
    RegionColors = [230 76 0;255 212 128;230 153 0;255 0 212;
        255 204 246;179 25 255;128 0 106;179 0 149;166 255 77;
        17 102 0;30 179 0;38 230 0;50 50 50;221 255 51;0 0 102;
        128 234 255;77 136 255;255 128 102;255 191 179;128 0 0;
        255 0 0];   
    RegionColors =RegionColors./255;
    colormap(RegionColors)

    for iregion=1:numGBDr
        plotm(sGBDRegions(iregion).Lat,sGBDRegions(iregion).Lon,'Color',[0 0 0]); %plot region borders
    end
    set(gcf,'Position',[342 378 850 396]);
    set(gca,'Position',[0.05 .05 0.9 0.9]);
    t = ntitle('GBD Countries and Regions','FontWeight','bold','FontSize',16);
    lg = legend(c,GBD_Regions,'Location','eastoutside','FontSize',12);
    title(lg,'Regions','FontWeight','bold','FontSize',12);
    print('-dpng','-r300','-loose','/data6/emcduffie/Figures/GBDRegions_Map.png');
end
end

fprintf('COMPLETE\n')
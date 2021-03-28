function [C,D,Met] = LoadData_xMAPS(GridRes,Nest,Data_Directory)
%clear
%---------------------------------------------------------------------
%                     LoadData_xMAPS.m, Oct. 16 2020
%          Script to Load Monthly PM2.5-relevant GEOS-Chem Data, 
%         Monthly OM/OC Rations, and calculate PM2.5 and aerosol components in ug m3
%                          Adpated for GBD-MAPS
%-------------------------------------------------------------------
% E. McDuffie
% General Overview
% Loads data sets for the specified months and days (month,dat,lon,lat,alt)
% Loads OMOC data (month,day,lon,lat,alt)
% Calculates PM2.5 from either the standard mechanism (PM2.5_ugsm3_calc)

%PROGRAM OPTIONS
inputFileType = 2;      %net CDF File type (1 = daily files, 2 = monthly files)
Month = 1:1:12;         %Month matrix
Day = 1;%1:1:31;        %Day matrix (set to 1 for monthly files)
Year = '2017';
h5file = 0;             %Are the data h5 files? (0 = no, 1 = yes) - only matters for inputfiletype 1
cSOA = 0;               %Did the sim use the Compelx SOA mechansim? (1 = yes, 0 = no)
PMLoad = 1;             % Load PM data (0 = no, 1 = yes)
PMCalc =1;              %Calculate PM2.5
GammaLoad = 0;          %Load N2O5 Uptake variables (0 = no, 1 = yes)
MetOnly = 0;            %only load state met file (1=yes, 0 = no)
METVAR = {'all'};       %specify the parameters you want to load, or set to all for pre-defined list (if PMLoad ==1)
Scripts_Directory = '/misc/data6/emcduffie/Matlab_Scripts';  %location of Matlab Scripts

%add appropriate paths. 1) make cell of paths, 2) add paths if necesssary
pathCell = regexp(path,pathsep,'split');
OnPath = any(strcmpi(Scripts_Directory,pathCell));
if OnPath==0;addpath(Scripts_Directory);end
DataPath = Data_Directory;

fprintf('Loading Data from %s\n for Month(s):%d-%d and Day(s):%d-%d, FileType =%d\n',...
    Data_Directory,Month(1),Month(end),Day(1),Day(end),inputFileType);

%0) Determine grid resolution
NumCoords = NaN(1,2); %(lat,lon)
switch GridRes
    case 1
        NumCoords(1)=91;
        NumCoords(2) =144;
    case 2
        switch Nest
            case 'na'
                NumCoords(1) = 121;
                NumCoords(2) = 161;
            case 'as'
                NumCoords(1) = 133;
                NumCoords(2) = 145;
            case 'eu'
                NumCoords(1) = 81;
                NumCoords(2) = 129;
        end
end

%1) Load all Chemical and Met Data
if PMLoad == 1
    disp('Loading Chemical Species and Met. Data...')
    [Lat,Lon,Met,D,~,~,~,~,~] = ExtractDataAll(inputFileType,DataPath,h5file,NumCoords,Month,Day,Year,cSOA,GammaLoad,MetOnly,METVAR);
    if GammaLoad ==0, clear yield_ClNO2 gamma1_N2O5 gamma2_N2O5 gamma3_N2O5; end   
else
    disp('Skipping PM Data Load...')
end
if PMCalc ==1
    %2) Load and Grid OM/OC Ratio Data
    disp('Loading OM/OC Data...')
    [mOMOC] = OMOC_Grid(GridRes,Nest,inputFileType);
    %3) Calculate PM2.5 (ugsm3), as well as components
    disp('Calculating PM2.5 and aerosol component concentrations...')
    [C] = PM25_Calc(D,Met,mOMOC,cSOA);
else
    disp('Skipping PM Calc...')
end
    
Met.Lat = Lat;
Met.Lon = Lon;
%clear extra variables
clear cSOA Day formatSpec inputFileType Month str GridRes AvgTmp c c2 Datatmp...
    Scripts_Directory Data_Directory cISOA EMLoad OnPath pathCell PMLoad time...
    PMCalc DataPath h5file GammaLoad SVPOA NumCoords Nest METVAR MetOnly ESPEC
disp('Load Complete')

%LOCAL FUNCTIONS
function [Lat,Lon,Met,D,gamma1_N2O5,gamma2_N2O5,gamma3_N2O5,yield_ClNO2,Aerosol_pH] = ExtractDataAll(inputFileType,DataPath,h5file,Res,Month,Day,Year,cSOA,GammaLoad,MetOnly,METVAR)
%extracts data for the specified file at the desired month and day
%Write 5D waves that for each parameter by: (Month, Day, Lat, Lon, Lev)
%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS
%Month and Day can be single values or vectors. The months and days must be
%listed in sequential order and be row vectors
%input file type: 1 or 2
% 1 = ts_24h (daily), 2 = net CDF diagnostic files (monthly)
% h5file (1 = yes , 0 = no) (only matters for input file type ==1)
%Res ([lat, lon] number)
%cSOA: 0 or 1
%1=complex SOA, 0 = simple SOA
%GammaLoad = 1, yes, 0 - no, Load N2O5 uptake variables?
%OUTPUTS
%Lat, lon, and specified Met parameters, and chemical concentrations
%Temp [K] - Temp interpolated to current time (from TMPU1 and TMPU2)
%RH [%]
%BoxHeight [m] - height of grid box
%PEdge [hPa] - surface pressure at level edges (based on moist air)
%PEdgeDry [hPa] - Surface pressure at level edges (based on dry air)
%pMid [hPa] - Pressure at mid point of each level (based on moist air)
%PMidDry [hPa] = pressure at mid-height of each level (based on dry air)
%AirVol [m3] - volume of dry air in a grid box 
%AirMass [kg] - mass of dry air in a grid box
%-PM25_ugsm3_nc [ug sm-3] - concentration of PM2.5 from netCDF diagnostic file
%Remaining are chemical compound mixing ratios in mol mol-1 air
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Initialize Matrices
lenM = 12;
if inputFileType ==1
    lenD = 31;
else
    lenD =1;
end

GCFileMat = MakeFileNameMatrixFn_xMAPS(inputFileType,DataPath,h5file, Month, Day,Year);  %get matrix of file names

%initialize Matrices in structures
%%% Met Parameters
if strcmp(METVAR,'all')
    MetList = {'AREAM2','T','AIRVOL','AD','PMIDDRY','PMID','PEDGEDRY','PEDGE','TropP','TropHt','PBLH','RH','TS','PSC2WET'};%,'RH','PEDGE','PEDGEDRY'};
else
    MetList = METVAR;
end
numMet = size(MetList,2);
for imet = 1:numMet
    if strcmp(MetList{imet},'PEDGE') || strcmp(MetList{imet},'PEDGEDRY')
        Met.(MetList{imet}) = repmat(NaN,[lenM,lenD,Res(2),Res(1),48]);
    elseif strcmp(MetList{imet},'TropP') || strcmp(MetList{imet},'TropHt')|| strcmp(MetList{imet},'TS')|| strcmp(MetList{imet},'PSC2WET')
        Met.(MetList{imet}) = repmat(NaN,[lenM,lenD,Res(2),Res(1),1]);
    else
        Met.(MetList{imet}) = repmat(NaN,[lenM,lenD,Res(2),Res(1),47]);
    end
end
%%% SPECIES
SpecList = {'NH4','NIT','SO4','BC','DST1','DST2','SALA','N2O5',...
                'ClNO2','HNO3','O3','NO','NO2','NH3','SO2','CO'};
numspec = size(SpecList,2);
for ispec=1:numspec
    D.(SpecList{ispec}) = repmat(NaN,[lenM,lenD,Res(2),Res(1),47]);
end
if cSOA ==1
	OrgList = {'TSOA','ISOA','ASOA','SOAGX','SOAMG','SOAIE','SOAME','LVOCOA','ISN1OA','INDIOL'};
    if SVPOA ==1
        POAList = {'OPOA1','OPOA2','POA1','POA2'};
    else
        POAList = {'OCPI','OCPO'};
    end
elseif cSOA ==0
	OrgList = {'SOAS'};
	POAList = {'OCPI','OCPO'};
end            
numOrg = size(OrgList,2);
numPOA = size(POAList,2);
for iorg=1:numOrg
    D.(OrgList{iorg}) = repmat(NaN,[lenM,lenD,Res(2),Res(1),47]);
end
for ipoa=1:numPOA
    D.(POAList{ipoa}) = repmat(NaN,[lenM,lenD,Res(2),Res(1),47]);
end

D.PM25_ugsm3_nc = repmat(NaN,[lenM,lenD,Res(2),Res(1),47]);

gamma1_N2O5 = repmat(NaN,[lenM,lenD,Res(2),Res(1),47]);   %N2O5+H2O --> 2HNO3
gamma2_N2O5 = repmat(NaN,[lenM,lenD,Res(2),Res(1),47]);   %N2O5+ HCl--> CLNO2+HNO3 (strat only)
gamma3_N2O5 = repmat(NaN,[lenM,lenD,Res(2),Res(1),47]);   %N2O5+sea salt --> ClNO2 + HNO3
yield_ClNO2 = repmat(NaN,[lenM,lenD,Res(2),Res(1),47]);
Aerosol_pH = repmat(NaN,[lenM,lenD,Res(2),Res(1),47]);

%%% Load Data
%only open data you actually want
lenM=size(Month,2);
lenD=size(Day,2);


if inputFileType == 2
    %If the Input files are the netCDF Diagnostic files...
    %The lat and lon will always be the same, so just pull from the first file
    if MetOnly ==1
        MetGCFileMat = strrep(GCFileMat(:,:),'SpeciesConc','StateMet');
        PARAM = sprintf('lat');
        Lat = double(ncread(MetGCFileMat{1,1}, PARAM));
        PARAM = sprintf('lon');
        Lon = double(ncread(MetGCFileMat{1,1}, PARAM));
    else
        PARAM = sprintf('lat');
        Lat = double(ncread(GCFileMat{1,1}, PARAM));
        PARAM = sprintf('lon');
        Lon = double(ncread(GCFileMat{1,1}, PARAM));
    end
    for iM=1:lenM
        if lenD >1
        end
        for iD=1:lenD
            %take the entire 3D matrix (all lat,lon,alt) of Temp from the file specified
            %by iM and iD in GCFileMat and write it to the Temp Matrix in the
            %correct month-row and day-column - result will be
            %5D matrix: month, day, SPEC Conc at: Spec_lat, Spec_lon, Spec_alt
            %Species
            if ~MetOnly==1
            for ispec=1:numspec
                if strcmp(SpecList{ispec},'BC')
                    SPEC = sprintf('SpeciesConc_%s','BCPI');
                    tmp1 = double(ncread(GCFileMat{iM,iD},SPEC));
                    SPEC = sprintf('SpeciesConc_%s','BCPO');
                    tmp2 = double(ncread(GCFileMat{iM,iD},SPEC));
                    D.BC(Month(iM),Day(iD),:,:,:) = nanmean(tmp1(:,:,:,1:end),4)+nanmean(tmp2(:,:,:,1:end),4);
                else
                    SPEC = sprintf('SpeciesConc_%s',SpecList{ispec});
                    tmp = double(ncread(GCFileMat{iM,iD},SPEC));
                    D.(SpecList{ispec})(Month(iM),Day(iD),:,:,:) = nanmean(tmp(:,:,:,1:end),4); %average time to 1-mo average
                end
            end

            for iorg=1:numOrg
                if strcmp(OrgList{iorg},'TSOA')
                    SPEC = sprintf('SpeciesConc_%s','TSOA0');
                    tmp1 = double(ncread(GCFileMat{iM,iD},SPEC));
                    SPEC = sprintf('SpeciesConc_%s','TSOA1');
                    tmp2 = double(ncread(GCFileMat{iM,iD},SPEC));
                    SPEC = sprintf('SpeciesConc_%s','TSOA2');
                    tmp3 = double(ncread(GCFileMat{iM,iD},SPEC));
                    SPEC = sprintf('SpeciesConc_%s','TSOA3');
                    tmp4 = double(ncread(GCFileMat{iM,iD},SPEC));
                    D.(OrgList{iorg})(Month(iM),Day(iD),:,:,:) = tmp1+tmp2+tmp3+tmp4;
                elseif strcmp(OrgList{iorg},'ISOA')
                     SPEC = sprintf('SpeciesConc_%s','ISOA1');
                    tmp1 = double(ncread(GCFileMat{iM,iD},SPEC));
                    SPEC = sprintf('SpeciesConc_%s','ISOA2');
                    tmp2 = double(ncread(GCFileMat{iM,iD},SPEC));
                    SPEC = sprintf('SpeciesConc_%s','ISOA3');
                    tmp3 = double(ncread(GCFileMat{iM,iD},SPEC));
                    D.(OrgList{iorg})(Month(iM),Day(iD),:,:,:) = tmp1+tmp2+tmp3;
                elseif strcmp(OrgList{iorg},'ASOA')
                    SPEC = sprintf('SpeciesConc_%s','ASOAN');
                    tmp1 = double(ncread(GCFileMat{iM,iD},SPEC));
                    SPEC = sprintf('SpeciesConc_%s','ASOA1');
                    tmp2 = double(ncread(GCFileMat{iM,iD},SPEC));
                    SPEC = sprintf('SpeciesConc_%s','ASOA2');
                    tmp3 = double(ncread(GCFileMat{iM,iD},SPEC));
                    SPEC = sprintf('SpeciesConc_%s','ASOA3');
                    tmp4 = double(ncread(GCFileMat{iM,iD},SPEC));
                    D.(OrgList{iorg})(Month(iM),Day(iD),:,:,:) = tmp1+tmp2+tmp3+tmp4; 
                else
                    SPEC = sprintf('SpeciesConc_%s',OrgList{iorg});
                    D.(OrgList{iorg})(Month(iM),Day(iD),:,:,:) = double(ncread(GCFileMat{iM,iD}, SPEC));
                end
            end
            for ipoa=1:numPOA
                SPEC = sprintf('SpeciesConc_%s',POAList{ipoa});
                D.(POAList{ipoa})(Month(iM),Day(iD),:,:,:) = double(ncread(GCFileMat{iM,iD}, SPEC));
            end
            
            %Aerosol Mass
            AerosolGCFileMat = strrep(GCFileMat(:,:), 'SpeciesConc', 'Aerosols');
            if GammaLoad ==1
                SPEC = sprintf('Chem_GammaN2O5H2O');
                gamma1_N2O5(Month(iM),Day(iD),:,:,:) = double(ncread(AerosolGCFileMat{iM,iD}, SPEC));
                SPEC = sprintf('Chem_GammaN2O5HCl');
                gamma2_N2O5(Month(iM),Day(iD),:,:,:) = double(ncread(AerosolGCFileMat{iM,iD}, SPEC));
                SPEC = sprintf('Chem_GammaN2O5SS');
                gamma3_N2O5(Month(iM),Day(iD),:,:,:) = double(ncread(AerosolGCFileMat{iM,iD}, SPEC));
                SPEC = sprintf('Chem_YieldClNO2');
                yield_ClNO2(Month(iM),Day(iD),:,:,:) = double(ncread(AerosolGCFileMat{iM,iD}, SPEC));
            end
            end
            
            %Met Data
            MetGCFileMat = strrep(GCFileMat(:,:),'SpeciesConc','StateMet');
            PresGCFileMat = strrep(GCFileMat(:,:), 'SpeciesConc', 'LevelEdgeDiags');
            for imet =1:numMet 
                PARAM = sprintf('Met_%s',MetList{imet});
                if strcmp(MetList{imet},'PEDGEDRY') || strcmp(MetList{imet},'PEDGE')
                    tmp = double(ncread(PresGCFileMat{iM,iD}, PARAM));
                    Met.(MetList{imet})(Month(iM),Day(iD),:,:,:) = squeeze(nanmean(tmp(:,:,:,1:end),4));
                elseif strcmp(MetList{imet},'PBLH') || strcmp(MetList{imet},'TS') || strcmp(MetList{imet},'PSC2WET')
                    tmp = double(ncread(MetGCFileMat{iM,iD}, PARAM));
                    Met.(MetList{imet})(Month(iM),Day(iD),:,:,1) = squeeze(nanmean(tmp(:,:,1:end),3));%has only 1 alt dimension
                else
                    tmp = double(ncread(MetGCFileMat{iM,iD}, PARAM));
                    Met.(MetList{imet})(Month(iM),Day(iD),:,:,:) = squeeze(nanmean(tmp(:,:,:,1:end),4));
                end
            end
        end
    end
end
end
function [C] = PM25_Calc(D,Met,mOMOC,cSOA)
%------------------------------------------------
%PM25_CalcFn reads individual specie data from the inout file and calcualtes PM2.5 at STP or ambient (volumetric) conditions - E. McDuffie, 10/17/2018
%General Equation: 
%   PM2.5[μg m^-3]=1.10 × (SO_4^(2-)+NO_3^-+NH_4^+ )+BC+DST1+0.38×DST2+1.86×SALA+1.02×SOAS+OM:OC × (OCPI+1.02×OCPO)
%Simulated mixing ratios of organic matter are converted into mass of organic carbon
%   using spatially and temporally explicit OM:OC ratios from Philip, et al. Atmos. Environ (2014).
%To account for aerosol hygroscopic growth, we use hygroscopicity factors at 35% RH
%   from Latimer and Martin, ACP (2019) for secondary inorganic (1.1) and organic 
%   (1.02) aerosol, and use the GEOS-Chem recommended value for sea salt (1.86). 
% 38% of the DST2 compoennt is added following GEOS-Chem Wiki recommendations
% Individual components are provided as output from the GEOS-Chem model as 
%   average monthly mixing ratios and are converted into units of micrograms
%   per m3 using ambient surface temperate (Tamb) and pressure (Pamb) from the MERRA-2 dataset

%INPUTS
%   NH4, NIT, SO4, BC, OCPO,OCPI,SOAS,DST1,DST2,SALA - 5D time matrices
%   cSOA; 0 = simple SOA (trop chem mech), 1=complexSOA mech
%   Model data in the input file is in mixing ratio at ambinet conditions
%   PM2.5 = 
%OUTPUTS
%   C, structure with PM2.5 component concentrations in ugsm3 - 5D matrix
%---------------------------------------------
clear PM25_ugsm3_calc C

%Define STP as 298K, 1 atm
VOL=1;              %VOL = 0, Standard conditions, VOL=1, Local Conditions
STD_P = 1013.25;    %hPa = Pa/100
STD_T = 298.0;      %K, ask why this isnt 273K
R = 8.314;          %Gas constant, [m3 Pa K-1 mol-1]
Amb_P = Met.PMID;   %hPa; pressure at grid box center (wet air) ~58m
Amb_T = Met.T;      %K; temperature at grid box center ~58m

%account for datatype 1 products are in ppbv, while datatype 2 are already in mixing ratio
tmp = squeeze(D.NH4(:,1,:,:,1));
tmp = nanmean(squeeze(nanmean(squeeze(nanmean(tmp,1)),1)),2);
if tmp > 1e-5
    if VOL ==1  
        factor = (Amb_P./(R.*Amb_T)).*100*.1e6.*1e-9;
    else
        factor = (STD_P/(R*STD_T))*100*1e6*1e-9;  %conversion factor in mol/m3, accounting for the 1e6 ug/g conversion factor needed to convert to ug
    end
else
    if VOL ==1
        factor = (Amb_P./(R.*Amb_T)).*100.*1e6;
    else
        factor = (STD_P/(R*STD_T))*100*1e6;
    end
end

if VOL==1
    vend = 'ugm3';
else
    vend='ugsm3';
end

%factor to convert mol mol-1 of each compound into ugm3
%factor = (AirMass./AirVol).*(1/MW_air).*1e6;  
%  mol spc                  (g/mol spc)        1
%  -------- * AirMass [kg] * ---------  * ------------ * 1e9 [kg to ug] or 1e6 for g to ug
%  mol air                  (g/mol air)    AirVol [m3]

%volumetric (ambient) to standard conversion factor
%Vol2Std = 1;%(STD_P./PMid).*(Temp./STD_T);
%disp(Vol2Std(12,1,73,46,1))

%MW in g/mol of NH4, Carbon, Nitrogen, Sulfate, Dust, SeaSalt, and SOA (SOAGX, SOAMG, SOAIE, SOAME, INDIOL, LVOCOA, ISN1OA)
MW_Aerosol =[18,12,62,96,29,36,150,58,72,118,102,102,154,226]; %changed SS MW from 31.4 to 36 g/mol (12/7/18)

%Hygroscopic Growth Factors
IGF = 1.1;%1.2%1.33; %inorganic hygroscopicity growth factor at 35% RH (1.51 at 50%)
OGF = 1.02;%1.16; %organic hygroscopicity growth factor at 35% RH (1.24 at 50%)
SSGF = 1.86; %sea salt hygroscopicity growth factor at 35% RH (2.42 at 50%)
OMOC = mOMOC;% 1.8; %global mean OM/OC ratio, recommended by Aerosol working group in 2016
fprintf('Hygrowscopic factors: SIA: %.2f, OGF: %.2f, SSA: %.2f\n',IGF,OGF,SSGF)

%PM2.5 Calculation - convert into ambient ugm3
C.(sprintf('NH4_%s',vend)) = D.NH4.*factor.*MW_Aerosol(1);          %Ammonium
C.(sprintf('BC_%s',vend)) = D.BC.*factor.*MW_Aerosol(2);            %sum of Hydrophillic and Hydrophobic Black Carbon
C.(sprintf('NIT_%s',vend)) = D.NIT.*factor.*MW_Aerosol(3);          %nitrate
C.(sprintf('SO4_%s',vend)) = D.SO4.*factor.*MW_Aerosol(4);          %sulfate
C.(sprintf('DST1_%s',vend)) = D.DST1.*factor.*MW_Aerosol(5);        %Dust  (accumulation mode)
C.(sprintf('DST2_%s',vend)) = D.DST2.*factor.*MW_Aerosol(5).*0.38;  %Dust  (coarse mode), only inlcude 38%
C.(sprintf('SALA_%s',vend)) = D.SALA.*factor.*MW_Aerosol(6);        %Sea Salt Aerosol (accumulation mode)

%inorganic (+black carbon) contribution to PM2.5
C.(sprintf('PM25_%s',vend)) = IGF.*(C.(sprintf('NH4_%s',vend))+C.(sprintf('NIT_%s',vend))+C.(sprintf('SO4_%s',vend)))...
    +C.(sprintf('BC_%s',vend))+C.(sprintf('DST1_%s',vend))+C.(sprintf('DST2_%s',vend))+SSGF.*C.(sprintf('SALA_%s',vend));

%OA
if cSOA ==0
    %Simple SOA
    SOAS = D.SOAS.*factor.*MW_Aerosol(7);   
    OCPO = D.OCPO.*factor.*MW_Aerosol(2);        %Hydrophobic Organic Carbon
    OCPI = D.OCPI.*factor.*MW_Aerosol(2);        %Hydrophillic Organic Carbon
    %Simple contribution to PM2.5
    C.(sprintf('PM25_%s',vend)) = C.(sprintf('PM25_%s',vend))+(OGF.*SOAS)+(OMOC.*(OCPO+OGF.*OCPI));
    C.(sprintf('SOA_%s',vend)) = OGF.*SOAS;
    C.(sprintf('POA_%s',vend)) = OMOC.*(OCPO+OGF.*OCPI);
end

%convert species to wet concentrations (35%)
C.(sprintf('NH4_%s',vend)) = C.(sprintf('NH4_%s',vend)).*IGF;%.* Vol2Std;
C.(sprintf('NIT_%s',vend)) = C.(sprintf('NIT_%s',vend)).*IGF;% .* Vol2Std;
C.(sprintf('SO4_%s',vend)) = C.(sprintf('SO4_%s',vend)).*IGF;% .* Vol2Std;
C.(sprintf('SALA_%s',vend)) = C.(sprintf('SALA_%s',vend)).*SSGF;% .* Vol2Std;
C.(sprintf('DST_%s',vend)) = C.(sprintf('DST1_%s',vend))+C.(sprintf('DST2_%s',vend));
%POA_ugsm3 = POA_ugsm3;% .* Vol2Std; %already has water component
%SOA_ugsm3 = SOA_ugsm3;% .* Vol2Std; %already has water component

end
function [mOMOC] = OMOC_Grid(GridRes,Nest,inputFileType)
%----------------------------------------------------------
%Makes matrix of OM/OC ratios, interpolated to current Grid Res 
%11/20/18 - EEM
%INPUTS
%   GridRes- Grid Resolution, 1 = (2x2.5), 2 = 0.5x0.625, 3 = (0.01x0.01)
%   Nest - nested region
%   inputFileType: 1 = 24 hour data, 2 = monthly data
%   Also need location of GEOS-Chem OMOC input data and a separate .mat
%   file with Lat and Lon vectors corresponding to different model
%   resolutions
%OUTPUTS
%   mOMOC - (12:31:Lon:Lat) matrix of OM/OC Ratios
%---------------------------------------------------------

fprintf('OMOC Grid = %d\n',GridRes);

% construct Grid (centre points) as per GEOS-Chem manual, 'Horizontal Grids'
%First,choose desired region to evaluate surface observations
%LatRegion = [-90 90];
%LonRegion = [-180 178]; 
nLev = 47;  %number of elevation altitudes in GEOS-Chem simulation
grid_infile = '/data6/emcduffie/ExtData/MatFiles/Coordinates.mat';   %file that defines Lat and Lon vectors at different resolutions

switch GridRes
    case 1
        S=load(grid_infile, 'Lat2','Lon25');  %load the appropriate lat and lon from the coordinates file
        gcLat = S.Lat2;gcLon = S.Lon25;
    case 2
        switch Nest
            case 'na'
                S=load(grid_infile, 'Lat05_na','Lon0625_na');
                gcLat = S.Lat05_na;gcLon = S.Lon0625_na;
            case 'as'
                S=load(grid_infile, 'Lat05_as','Lon0625_as');
                gcLat = S.Lat05_as;gcLon = S.Lon0625_as;
            case 'eu'
                S=load(grid_infile, 'Lat05_eu','Lon0625_eu');
                gcLat = S.Lat05_eu;gcLon = S.Lon0625_eu;
        end
end
         
         nLat = numel(gcLat);
         nLon = numel(gcLon);

%Load OM/OC data 
FileLoc = '/misc/data10/ctm/HEMCO/OMOC/v2018-01/OMOC.';     %local location of GEOS-Chem OMOC input files
FileEnd = '.01x01.nc';

daysinmonth = [31 28 31 30 31 30 31 31 30 31 30 31];
%Make [M:D:Lat:Lon] matrix of OM/OC       
if inputFileType ==1
    mOMOC = NaN(12,31,nLon,nLat,nLev);
else
    mOMOC = NaN(12,1,nLon,nLat,nLev);
end
lenD = size(mOMOC,2);
MONTHS = 1:12;    
for Mn=1:numel(MONTHS)
	% Load monthly omoc conversion factor
    if find(any(Mn==[12 1 2])); Season = 'DJF';
    elseif find(any(Mn==3:5)); Season = 'MAM';
    elseif find(any(Mn==6:8)); Season = 'JJA';
    elseif find(any(Mn==9:11)); Season = 'SON';
    end

    fname = strcat(FileLoc,Season,FileEnd); % sprintf('/home/junmeng/project/junmeng/OMOC/OMOC.%s.01x01.nc',Season);
    latom = ncread(fname,'lat'); 
    lonom = ncread(fname,'lon');
    OMOC = ncread(fname,'OMOC')';
    
    if ~(GridRes ==3)           %strcmp(GRID,'01x01') - if the GC Data are not at the same resolution as the OMOC data, interp the OMOC data to the new grid
    	OMOC = interpfn(latom,lonom,OMOC,gcLat,gcLon);   
        OMOC(:,1)=OMOC(:,2);    %bug - first column is written as NaNs, set values to second column instead
        OMOC = OMOC';%
    end
    
    %fill in the day values, set extra days to zero, also set altitude data
    if lenD >1
        for iday=1:31
            if iday <= daysinmonth(Mn)
                for ialt=1:nLev
                    mOMOC(MONTHS(Mn),iday,:,:,ialt) = OMOC;
                end
            end
        end
    else
        for ialt=1:nLev
            mOMOC(MONTHS(Mn),1,:,:,ialt) = OMOC;
        end
    end
end
end
end
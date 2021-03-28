function [Frac_Ess_aavg]=B_Calc_FractionalSourceContributions_GridFn(Cases,C_merge,totalPM_file)
% E. McDuffie, last updated: Nov. 8, 2020
% Calculate the gridded (0.5x0.5) fractional contributions (for source
% sectors or total fuel types - looks to case names to make that decision)
%
% INPUTS:
%   Cases (sensitivity simulation names)
%   C_merge - merged PM2.5 data from GEOS-Chem
%   totalPM_file, (string), file that will store the sum
%   of PM2.5 mass from all sensitivity simulations. 
% OUTPUT: 
%   Frac_Ess_aavg: gridded fractional contributions for each sensitivty simulation case
%%%%%

%%%Step 0. Define all cases, calculate those not simulated, define 'Total_Contribution'
for doloop=1:1
use_fuel =1; %use detailed sector & fuel simulations
Csims = C_merge;
%calculate the cases that were not simulated (non-coal ENE, non-coal IND,non-coal and non-biofuel RCOR)
if contains('noENEcoal',Cases)==1 && use_fuel==1
    Csims.noENEother.PM25_ugm3 = Csims.base.PM25_ugm3-(Csims.noENEcoal.PM25_ugm3- Csims.noENE.PM25_ugm3);    %calculate the non-coal contribution (not quite right, but quicker than another sim)
    Csims.noINDother.PM25_ugm3 = Csims.base.PM25_ugm3-(Csims.noINDcoal.PM25_ugm3- Csims.noIND.PM25_ugm3);     %NEED to check
    Csims.noRCORother.PM25_ugm3 = 2.*Csims.base.PM25_ugm3+Csims.noRCOR.PM25_ugm3- Csims.noRCORcoal.PM25_ugm3-Csims.noRCORbiofuel.PM25_ugm3;
    allCases = {'Base','noAGR','noENEcoal','noENEother','noINDcoal',...
     'noINDother','noNRTR','noROAD','noRCORcoal','noRCORbiofuel','noRCORother','noRCOC',...
     'noRCOO','noSLV','noWST','noSHP','noGFEDagburn','noGFEDoburn','noAFCID','noWDUST','noNAT'};
elseif contains('noBIOFUEL',Cases)==1
    Csims.noOTHER.PM25_ugm3 = zeros(size(Csims.base.PM25_ugm3,3),size(Csims.base.PM25_ugm3,4));                     %calculate this remained below (calc'd product, not simulated value)
     allCases = {'Base','noBIOFUEL','noCOAL','noOILGAS','noWDUST','noAFCID','noGFEDagburn','noGFEDoburn','noOTHER'};
end

%clear any previous results
clear C_Abs Total_Contribution Frac_Ess_aavg
numcases = size(allCases,2);

%if calculating the total fuel specific results, need to load the
%'Total_Contribution' data that is calcualted from the detailed sectoral
%simualtions. This ensures that the contributions of the sub-sectors and
%total fuel types are each calculated relative to the same total.
if contains('noBIOFUEL',Cases)==1
    load(totalPM_file,'Total_Contribution') %(need to use the same values at the sectors or else the percentages don't compare - due to non-linear chem)
    Sum = zeros(size(Csims.base.PM25_ugm3,3),size(Csims.base.PM25_ugm3,4));
else
    Total_Contribution = zeros(size(Csims.base.PM25_ugm3,3),size(Csims.base.PM25_ugm3,4));
end
end

%Step 1) Calculate the gridded difference between the base and sensitivity simulation
for doloop=1:1
    %Grid =5; %data are 0.5 x0.5
    base = squeeze(nanmean(Csims.base.PM25_ugm3(:,1,:,:),1));                   %get annual average
    base(isnan(base))=0;
    C_Abs.(allCases{1}) = base;
    for icase=2:numcases
        %for all cases other than the base...
        fprintf('%s\n',allCases{icase})
        simtmp = squeeze(nanmean(Csims.(allCases{icase}).PM25_ugm3(:,1,:,:),1));%get annual average
        simtmp(isnan(simtmp))=0;                                                %set NaN values to 0
        AbsContribution = base - simtmp;                                        % calculate the simulation - base difference
        C_Abs.(allCases{icase}) = AbsContribution;                              % stroe absolute contribution
        if contains('noENEcoal',Cases)==1 && use_fuel==1
            Total_Contribution = Total_Contribution + AbsContribution;          %only recalculate the total frac if the total sectors w/fuel categories are being calc'd
        else 
            if icase < numcases
                Sum  = Sum + AbsContribution;                                   % keep a running sum of the absolute contributions, but don't include 'OTHER' (fuel cases only)
            end
        end
    end
end

% Step 2) Calculate Gridded Fractional Contributions (at 0.5 x 0.5 global gridded resolution)
% fractional contribution = absolute contribution / total contribution
for doloop=1:1
    for icase=1:numcases
        if icase==1
            nametmp = allCases{icase};
            Frac_Ess_aavg.(char(nametmp)) = single(Total_Contribution./Total_Contribution); %**base is set to 1
        else
            nametmp = allCases{icase}(3:end);
            Frac_Ess_aavg.(char(nametmp)) = single(C_Abs.(allCases{icase})./Total_Contribution); 
        end
    end
end

%Step 3) Save the total_contribution matrix for future use (when calculating sectors)
% OR calculate the remaining fraction for the fuel simulations
for doloop=1:1
if contains('noENEcoal',Cases)==1 && use_fuel==1
    save(totalPM_file,'Total_Contribution')
else
    Frac_Ess_aavg.OTHER = single((Total_Contribution - Sum)./Total_Contribution);       %calculate the remaining fraction for the fuel simulations
    Frac_Ess_aavg.OTHER(Frac_Ess_aavg.OTHER > 1) =1;
    Frac_Ess_aavg.OTHER(Frac_Ess_aavg.OTHER < 0) =0;                                    %make sure all fuel types + OTHER still totals to the Total_Contribution
end
end
end
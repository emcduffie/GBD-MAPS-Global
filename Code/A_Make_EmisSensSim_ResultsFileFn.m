function [C_merge,C_merge_aavg]=A_Make_EmisSensSim_ResultsFileFn(Directory,outFile,cases,grid_infile)

% Intermediate function where the interpolation setting and species list can be set. 
% E. McDuffie, last updated Nov. 8, 2020
%
% Dependencies:
%   Merge_NestedGlobalFn.m
%%%%%%%

interp=0;   %interpolate data when re-gridding? (0 = no, 1 = yes)
speclist = {'PM25'}; %specify the compounds to load and merge. 
% Complete list = %,{'PM2.5','NIT','NH4','BC','SO4','SALA','POA','SOA','DST'};

%Load and merge all sensitivity simulations to uniform 0.5x0.5 grid
[C_merge,C_merge_aavg]=A1_Merge_GlobalNestedFn(Directory,cases,speclist,interp,grid_infile);  %produces C_merge.{case}.{spec}

%save gridded PM.25 results for each source simulation
save(outFile,'C_merge','C_merge_aavg','-v7.3')
fprintf('COMPLETE: Saved File %s\n',outFile)
end


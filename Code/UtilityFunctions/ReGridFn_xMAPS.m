function [outData,outLat,outLon] = ReGridFn_xMAPS(Grid,Nest,interp,inData,inLat,inLon)
%Puts the data from the inLat/ iLon grid, onto a new specified grid

%E. McDuffie - Oct 16, 2020
%INPUTS
 % Grid = 0 = 0.01x0.01,1 = 0.1x0.1, 2 = 0.5x0.5, 3 = 2x2.5, 4 = 0.5x0.625 (nested)
 % interp = 0 or 1 (0 = retain original grid, 1 = interpolate )
 % in Data = data to regrid (***lon x lat***)
 % inLat = lat of input data
 % inLon = lon of input data
 %OUTPUTS
 % outData = data (lon x lat) to new grid
 % outLat = Lat of re-gridded data
 % outLon = Lon of re-gridded data
 % *** NOTE - If Interp ==0:
 %      When gridding to coarser resolution  - averages of sub box placed
 %      into larger box
 %      When gridding to finer resolution - coarse box places in all
 %      sub-boxes. This means that the AVERAGE of the inData and outData
 %      will be the same, but the SUM of outData will be a factor of X
 %      larger than the SUM of inData, where X = the number of sub-boxes
 %      within each coarse grid box.

switch Grid
    case 0
        numlat = 18000;
        numlon = 36000;
        res = [0.01,0.01];
        outLat = (-90+(res(1)/2):res(1):90-(res(1)/2))';
        outLon = (-180+(res(2)/2):res(2):180-(res(2)/2))';
    case 1
        numlat = 1800;
        numlon = 3600;
        res = [0.1,0.1];
        outLat = (-90+(res(1)/2):res(1):90-(res(1)/2))';
        outLon = (-180+(res(2)/2):res(2):180-(res(2)/2))';
    case 2
        numlat = 360;
        numlon = 720;
        res = [0.5,0.5];
        outLat = (-90+(res(1)/2):res(1):90-(res(1)/2))';
        outLon = (-180+(res(2)/2):res(2):180-(res(2)/2))';
    case 3
        numlat = 91;
        numlon = 144;
        res = [2,2.5];
        outLat = (-90:res(1):90)';
        outLat(1) = outLat(1)+0.5;
        outLat(end) = outLat(end)-0.5;
        outLon = (-180:res(2):180-(res(2)/2))';
    case 4
        res = [0.5,0.625];
        switch Nest
            case 'na'
                numlat = 121;
                numlon = 161;
                LatBounds = [10.0,70.0];
                LonBounds = [-140.0,-40.0];

            case 'as'
                numlat = 133;
                numlon = 145;
                LatBounds = [-11.0,55.0];
                LonBounds = [60.0,150.0];
                
            case 'eu'
                numlat = 81;
                numlon = 129;
                LatBounds = [30.0,70.0];
                LonBounds = [-30.0,50.0];
        end
        outLat = (LatBounds(1):res(1):LatBounds(2))';
        outLon = (LonBounds(1):res(2):LonBounds(2))';
end

outData = repmat(0,size(outLon,1),size(outLat,1)); %initiailize new data matrix

%only do calcualtion if the inData is at a different resolution than the requested out resolution
if ~(size(inLat,1) ==size(outLat,1)) || ~(size(inLon,1)==size(outLon,1)) 
    if size(inLat,1)==18000; inres = [0.01,0.01]; istart=0;
    elseif size(inLat,1)==1800;inres=[0.1,0.1];istart=0;
    elseif size(inLat,1)==360; inres=[0.5,0.5];istart=0;
    elseif size(inLat,1)==133; inres=[0.5,0.625];istart=10;
    elseif size(inLat,1)==81; inres=[0.5,0.625];istart=10;
    elseif size(inLat,1)==121; inres=[0.5,0.625];istart=10;
    elseif size(inLat,1)==91; inres=[2,2.5];istart=0;
    end
   
if interp ==0
    if inres(1)*inres(2) < res(1)*res(2)    %if scaling up...
        difflat = res(1)/2;
        difflon= res(2)/2;
     %for all data points in new grid, find original grid pnt that is closest,
     %then fill that value into the new 'outData' Grid
     % i.e. retain grid boxes at new resolution
    for i=istart+1:numlat-istart
        for j=istart+1:numlon-istart
            outData(j,i) = nanmean(nanmean(inData(abs(outLon(j)-inLon)<=difflon,abs(outLat(i)-inLat)<=difflat),1),2);     %average all the subgrid data and place into new grid cell
        end
    end
    elseif inres(1)*inres(2) > res(1)*res(2)    %if scaling down to finer resolution...
        difflat = inres(1)/2;
        difflon= inres(2)/2;
        for i=istart+1:size(inData,2)-istart
            latlim = find(abs(outLat-inLat(i))<=difflat);
            for j=istart+1:size(inData,1)-istart
                lonlim = find(abs(outLon-inLon(j))<=difflon);
                outData(lonlim(1):lonlim(end),latlim(1):latlim(end)) = inData(j,i);
            end
        end
    end
elseif interp ==1
    %interpolate data onto new grid. 
    outData = interpfn_xMAPS(inLat,inLon,inData',outLat,outLon)'; %in (lat x lon), out (lat x lon)
    %grid data is slower
    %[outLAT,outLON]=meshgrid(outLat,outLon);
    %outData = griddata(inLat,inLon,inData,outLAT,outLON,'linear');
end
else
    fprintf("Data already at requested resolution (%sx%s)\n",num2str(res(1),'%.2f'),num2str(res(2),'%.2f'));
    outData = inData;
    outLat=inLat;
    outLon=inLon;
end
end
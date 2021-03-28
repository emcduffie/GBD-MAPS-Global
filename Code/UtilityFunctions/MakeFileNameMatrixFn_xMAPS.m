function GCFileMat = MakeFileNameMatrixFn (inputFileType,DataPath,h5file,Month,Day,Year)
%makes matrix of file names that you want to load data for with the
%ExtractDataFn
%INPUTS
%Month and Day can be single values or scalars. The months and days must be
%listed in sequential order and be row vectors
%input file type: 1 or 2
% 1 = ts_24h, 2 = trac net CDF
%h5file = 0 or 1 (no or yes)
%OUTPUTS
%Matrix with File Names
if inputFileType ==1
    FileStart = strcat(DataPath,sprintf('ts_24h_avg.%s',Year));
    if h5file ==0
        FileEnd = '.nc';
    else
        FileEnd = '.h5';
    end
elseif inputFileType ==2
    FileStart = strcat(DataPath,sprintf('GEOSChem.SpeciesConc.%s',Year));
    FileEnd = '_0000z.nc4';
else
    error('ExtractData: File Type not supported');
end
M = Month;  
D = Day;
lenM = size(M,2);%get the row diemnsion of M and D
lenD = size(D,2);
GCFileMat= cell(lenM, lenD);   %initialize GCFile name matrix with the correct dimensions

if lenM > 1 %if month is a matrix, calculate the list of file names for all months
    for iM=1:lenM
        if M(iM) <=9
            strm=strcat('0',num2str(M(iM)));%convert month number to month string
        else
            strm=num2str(M(iM));
        end
    if lenD > 1 %if day is a matrix, 
        if M(iM) == 1||M(iM) ==3||M(iM) ==5||M(iM) ==7||M(iM) == 8||M(iM)==10||M(iM)==12
            lenD =31-D(1)+1;
        elseif M(iM) ==4||M(iM)==6||M(iM)==9||M(iM)==11
            lenD =30-D(1)+1;
        elseif M(iM)==2
            lenD=28-D(1)+1;
        end
        for iD=1:lenD
            if D(iD) <=9   %format for 09 and 10 and above
                GCFileMat{iM,iD} = strcat(FileStart,strm,sprintf('0%d',D(iD)),FileEnd);
            else
                GCFileMat{iM,iD} = strcat(FileStart,strm,sprintf('%d',D(iD)),FileEnd);
            end
        end
    elseif lenD==1
        iD = 1;
        if D(iD) <=9   %format for 09 and 10 and above
            GCFileMat{iM,iD} = strcat(FileStart,strm,sprintf('0%d',D(iD)),FileEnd);
        else
           GCFileMat{iM,iD} = strcat(FileStart,strm,sprintf('%d.nc',D(iD)),FileEnd);
        end
    else
        error('ExtractData: Day not vector or scalar');
    end
    end
elseif lenM ==1
    %iM = iM+1;
    iM = 1;
    if M(iM) <=9
        strm=strcat('0',num2str(M(iM))); %convert month number to month string
    else
        strm=num2str(M(iM));
    end
    if lenD > 1 %if day is a matrix, 
        for iD=D(1):lenD
            if D(iD) <=9   %format for 09 and 10 and above
                GCFileMat{iM,iD} = strcat(FileStart,strm,sprintf('0%d',D(iD)),FileEnd);
            else
                GCFileMat{iM,iD} = strcat(FileStart,strm,sprintf('%d',D(iD)),FileEnd);
            end
        end
    elseif lenD==1
        iD = 1;
        if D(iD) <=9   %format for 09 and 10 and above
            GCFileMat{iM,iD} = strcat(FileStart,strm,sprintf('0%d',D(iD)),FileEnd);
        else
           GCFileMat{iM,iD} = strcat(FileStart,strm,sprintf('%d',D(iD)),FileEnd);
        end
    else
        error('ExtractData: Day not vector or scalar');
    end
else
    error('ExtractData: Month is not a vector or scalar');
end
end

    
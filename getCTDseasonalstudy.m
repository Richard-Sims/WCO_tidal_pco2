function [CTDtemp, CTDdepthm,CTDsalinity ,CTDDT,CTDdepth_interp,CTDTemperature_interp,CTDsalinity_interp,CTDchl,CTDchl_interp] = getCTDseasonalstudy(fullpathCTD,Licor_Datetime,j,fList);
%get ctd for seasonal study

% % % %Add exception for first deployment using seabird
if fList(j,1:10)==('2016-04-27');
SeabirdData =dlmread(fullpathCTD,'\s',163,0); %open ctd text file
CTDdepthm= SeabirdData(1:end,1); %define CTD salinity PSU
CTDsalinity = SeabirdData(1:end,2); %define ctd temperature
CTDtemp = SeabirdData(1:end,3); %define ctd conductivity
strdate=('2016-04-27 11:02:12');dt=datenum(strdate,'yyyy-mm-dd HH:MM:SS');
timearray=(1:length(SeabirdData(:,1)))'; timearrayDT=timearray/(24*3600*4); %sampling frequency 4Hz
CTDDT=dt+timearrayDT;
CTDchl=CTDsalinity;CTDchl_interp=CTDsalinity; %note as no chl make it equal salinity here to stop code falling over later
else
%read rbr files
rbrData =dlmread(fullpathCTD,'',54,2);%read data
CTDsalinity= rbrData(1:end,6); %define CTD salinity
CTDtemp = rbrData(1:end,2); %define ctd ctemp
CTDdepthm= rbrData(1:end,5); %define CTD salinity PSU
CTDchl= rbrData(1:end,4); %define CTD schlorphyll
clear rbrData
%import starttime
delimiter = ' ';startRow = 51;endRow = 51;formatSpec = '%*s%s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
fileID = fopen(fullpathCTD,'r');textscan(fileID, '%[^\n\r]', startRow-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false);fclose(fileID);
CTDhms = cell2mat(dataArray{:, 1}); CTDymd = fList(j,[1:4 6:7 9:10]); CTDymdhms=[CTDymd,CTDhms];
ghj=str2num(CTDymdhms(18:20));CTDDT_start= datenum(CTDymdhms,'yyyymmddHH:MM:SS')+(ghj/(3600*24*100)); 
%populate times for data using start time
samp_frequency=6; 
deploytimesecs=length(CTDdepthm)/samp_frequency; CTDDT=[];CTDtimevec=[];
CTDtimevec=(1:1:length(CTDdepthm))';%make a matrix of time
CTDDT=(CTDtimevec*((1/samp_frequency)/(3600*24)))+CTDDT_start; %convert all times to datenum format
CTDchl_interp = interp1(CTDDT,CTDchl,Licor_Datetime);m = find(isnan(CTDchl_interp ))'; CTDchl_interp (m) = 0;
end


%exception needed for when sensor comes out of water and salinity drops
%NOTE- these corrections were written before the salinity offset was available.
CTDDTcopy=CTDDT;
if  fList(j,1:10)==('2016-06-15');
    [delind,~]=find(CTDsalinity<35.295);
    CTDsalinity(delind) = [];
    CTDDTcopy(delind) = [];
elseif fList(j,1:10)==('2016-05-11');
    [delind,~]=find(CTDsalinity<35.34);
    CTDsalinity(delind) = [];
    CTDDTcopy(delind) = [];
elseif fList(j,1:10)==('2016-05-26');
    [delind,~]=find(CTDsalinity<35.3);
    CTDsalinity(delind) = [];
    CTDDTcopy(delind) = [];       
elseif fList(j,1:10)==('2016-06-10');
    [delind,~]=find(CTDsalinity<35.445);
    CTDsalinity(delind) = [];
    CTDDTcopy(delind) = [];        
elseif fList(j,1:10)==('2016-06-30');
    [delind,~]=find(CTDsalinity<35.25);
    CTDsalinity(delind) = [];
    CTDDTcopy(delind) = [];
elseif fList(j,1:10)==('2016-07-13');
    [delind,~]=find(CTDsalinity<35.35);
    CTDsalinity(delind) = [];
    CTDDTcopy(delind) = [];
elseif fList(j,1:10)==('2016-08-24');
    [delind,~]=find(CTDsalinity<35.16);
    CTDsalinity(delind) = [];
    CTDDTcopy(delind) = [];
elseif fList(j,1:10)==('2016-06-22');
    [delind,~]=find(CTDsalinity<35.2);
    CTDsalinity(delind) = [];
    CTDDTcopy(delind) = []; 
elseif fList(j,1:10)==('2016-08-17');
    [delind,~]=find(CTDsalinity<35.2);
    CTDsalinity(delind) = [];
    CTDDTcopy(delind) = []; 
end

%use calibration coefficents to calibrate CTD- NO calibration for CTD yet
CTDtemp=(CTDtemp*1)-0.005;
CTDsalinity=(CTDsalinity*1)-0.1;   

%use interp1 to linearly interpolate the salinity, temperature and depth such that they are at the 1HZ frequency as the co2 data. Where interpolating out of bounds replace nan with 0.
CTDdepth_interp = interp1(CTDDT,CTDdepthm,Licor_Datetime);k = find(isnan(CTDdepth_interp ))'; CTDdepth_interp (k) = 0;
CTDTemperature_interp = interp1(CTDDT,CTDtemp,Licor_Datetime);l = find(isnan(CTDTemperature_interp ))'; CTDTemperature_interp (l) = 0;
CTDsalinity_interp = interp1(CTDDTcopy,CTDsalinity,Licor_Datetime);m = find(isnan(CTDsalinity_interp ))'; CTDsalinity_interp (m) = 0;

%clearvars k l m filename delimiter endRow formatSpec fileID dataArray ans CTDdepth CTDtimevec CTDymdhms CTDhms CTDDT_start CTDDT_doy_start ValeportData J1 CTDDT_doy_start CTDhhms CTDymd fname path samp_frequency deploytimesecs CTDsalinity CTDdepthm CTDtemp CTDDT;
end

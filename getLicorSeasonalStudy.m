function [LicorH20B,LicorPressure,LicorTemperature,Licor_Datetime,Licor_Datetime_doy,LicorxCO2B ] = getLicorSeasonalStudy(fullpathLicor,j,fList)
%%%%%%%%%%%%%%%%%%%   Licor data inport %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%importpt100 text docs into workspace be selecting folder location
dataimp=dlmread(fullpathLicor,'\t',3, 1); %open licor text file
% % % %Add exception for Licor files formatted in the non default format
if  fList(j,1:10)==('2009-07-01');  
    LicorxCO2B=dataimp(:,2);
    LicorPressure=dataimp(:,18);
    LicorH20B=dataimp(:,5);
    LicorTemperature=dataimp(:,19);
else 
    LicorxCO2B=dataimp(:,17);%define licor xco2
    LicorPressure=dataimp(:,48);%define licor pressure
    LicorH20B=dataimp(:,41);%define licor water vapour content
    LicorTemperature=dataimp(:,53); %define licor temperature   
end

if fList(j,1:10)==('2016-04-27'); 
    strdate=('2016-04-27 10:31:07')
    dt=datenum(strdate,'yyyy-mm-dd HH:MM:SS');
    timearray=(1:length(dataimp(:,1)))';
    timearrayDT=timearray/(24*3600);
    Licor_Datetime=dt+timearrayDT;
elseif fList(j,1:10)==('2016-05-11')
    strdate=('2016-05-11 09:42:50')
    dt=datenum(strdate,'yyyy-mm-dd HH:MM:SS');
    timearray=(1:length(dataimp(:,1)))';
    timearrayDT=timearray/(24*3600);
    Licor_Datetime=dt+timearrayDT;
elseif fList(j,1:10)==('2016-08-24')
    strdate=('2016-08-24 08:56:21')
    dt=datenum(strdate,'yyyy-mm-dd HH:MM:SS');
    timearray=(1:length(dataimp(:,1)))';
    timearrayDT=timearray/(24*3600);
    Licor_Datetime=dt+timearrayDT;
elseif fList(j,1:10)==('2016-09-21')
    strdate=('2016-09-21 08:28:29')
    dt=datenum(strdate,'yyyy-mm-dd HH:MM:SS');
    timearray=(1:length(dataimp(:,1)))';
    timearrayDT=timearray/(24*3600);
    Licor_Datetime=dt+timearrayDT;
elseif fList(j,1:10)==('2016-09-15')
    strdate=('2016-09-15 09:11:04')
    dt=datenum(strdate,'yyyy-mm-dd HH:MM:SS');
    timearray=(1:length(dataimp(:,1)))';
    timearrayDT=timearray/(24*3600);
    Licor_Datetime=dt+timearrayDT;
elseif fList(j,1:10)==('2016-08-17')
    strdate=('2016-08-17 09:00:24')
    dt=datenum(strdate,'yyyy-mm-dd HH:MM:SS');
    timearray=(1:length(dataimp(:,1)))';
    timearrayDT=timearray/(24*3600);
    Licor_Datetime=dt+timearrayDT;
elseif fList(j,1:10)==('2016-09-02')
    strdate=('2016-09-02 09:37:29')
    dt=datenum(strdate,'yyyy-mm-dd HH:MM:SS');
    timearray=(1:length(dataimp(:,1)))';
    timearrayDT=timearray/(24*3600);
    Licor_Datetime=dt+timearrayDT;    
elseif fList(j,1:10)==('2016-08-10')
    strdate=('2016-08-10 09:34:17')
    dt=datenum(strdate,'yyyy-mm-dd HH:MM:SS');
    timearray=(1:length(dataimp(:,1)))';
    timearrayDT=timearray/(24*3600);
    Licor_Datetime=dt+timearrayDT;
elseif fList(j,1:10)==('2016-07-07')
    strdate=('2016-07-07 09:24:33')
    dt=datenum(strdate,'yyyy-mm-dd HH:MM:SS');
    timearray=(1:length(dataimp(:,1)))';
    timearrayDT=timearray/(24*3600);
    Licor_Datetime=dt+timearrayDT;
elseif fList(j,1:10)==('2016-07-13')
    strdate=('2016-07-13 09:43:36')
    dt=datenum(strdate,'yyyy-mm-dd HH:MM:SS');
    timearray=(1:length(dataimp(:,1)))';
    timearrayDT=timearray/(24*3600);
    Licor_Datetime=dt+timearrayDT;
elseif fList(j,1:10)==('2016-06-30')
    strdate=('2016-06-30 09:54:43')
    dt=datenum(strdate,'yyyy-mm-dd HH:MM:SS');
    timearray=(1:length(dataimp(:,1)))';
    timearrayDT=timearray/(24*3600);
    Licor_Datetime=dt+timearrayDT;
elseif fList(j,1:10)==('2016-06-22')
    strdate=('2016-06-22 09:48:28')
    dt=datenum(strdate,'yyyy-mm-dd HH:MM:SS');
    timearray=(1:length(dataimp(:,1)))';
    timearrayDT=timearray/(24*3600);
    Licor_Datetime=dt+timearrayDT;
elseif fList(j,1:10)==('2016-06-15')
    strdate=('2016-06-15 08:12:42')
    dt=datenum(strdate,'yyyy-mm-dd HH:MM:SS');
    timearray=(1:length(dataimp(:,1)))';
    timearrayDT=timearray/(24*3600);
    Licor_Datetime=dt+timearrayDT;
elseif fList(j,1:10)==('2016-05-26')
    strdate=('2016-05-26 09:42:42')
    dt=datenum(strdate,'yyyy-mm-dd HH:MM:SS');
    timearray=(1:length(dataimp(:,1)))';
    timearrayDT=timearray/(24*3600);
    Licor_Datetime=dt+timearrayDT;    
elseif fList(j,1:10)==('2016-08-04')
    strdate=('2016-08-04 10:05:42')
    dt=datenum(strdate,'yyyy-mm-dd HH:MM:SS');
    timearray=(1:length(dataimp(:,1)))';
    timearrayDT=timearray/(24*3600);
    Licor_Datetime=dt+timearrayDT;    
elseif fList(j,1:10)==('2016-06-10')
    strdate=('2016-06-10 08:54:58')
    dt=datenum(strdate,'yyyy-mm-dd HH:MM:SS');
    timearray=(1:length(dataimp(:,1)))';
    timearrayDT=timearray/(24*3600);
    Licor_Datetime=dt+timearrayDT;   
elseif fList(j,1:10)==('2016-07-20')
    strdate=('2016-07-20 09:13:52')
    dt=datenum(strdate,'yyyy-mm-dd HH:MM:SS');
    timearray=(1:length(dataimp(:,1)))';
    timearrayDT=timearray/(24*3600);
    Licor_Datetime=dt+timearrayDT;
elseif fList(j,1:10)==('2016-07-27');
    strdate=('2016-07-27 09:27:41');
    dt=datenum(strdate,'yyyy-mm-dd HH:MM:SS');
    timearray=(1:length(dataimp(:,1)))';
    timearrayDT=timearray/(24*3600);
    Licor_Datetime=dt+timearrayDT;     
else   
    %import Licor date and time seperately
    delimiter = '\t'; formatSpec = '%s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]'; fileID = fopen(fullpathLicor,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false); fclose(fileID);
    DateTime = dataArray{:, 1}; %define licor date time
    c = cellstr(DateTime(3:end,:));                %convert cell array to cell vector (selecting only data of interest)
    strDate = char(c);                             %convert each row of the cell vector to a character array
    Licor_Datetime = datenum(strDate,'yyyy-mm-dd HH:MM:SS');    %convert to matlab date - make sure form of date is specified correctly.
end
 %note this changes if run in different years
Jan1_serial = datenum([2016, 1, 1, 0, 0, 0]); %define first day of year
Licor_Datetime_doy=Licor_Datetime - Jan1_serial + 1;%    %convert date to doy (day of year)

end
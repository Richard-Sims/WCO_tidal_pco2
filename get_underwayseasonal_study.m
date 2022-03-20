clear all; close all;
% Import Underway data from text file.
% Script for importing data from the following text file:%    C:\Users\rps207\Documents\Seasonal Study Data\Quest Underway
addpath('c:\Users\rps207\Documents\Matlab\Functions');
addpath('c:\Users\rps207\Documents\Matlab\Functions\Add_Axis');
addpath('c:\Users\rps207\Documents\Matlab\Functions\cbdate');
addpath('c:\Users\rps207\Documents\Matlab\Functions\mixing_library');
addpath('c:\Users\rps207\Documents\Matlab\Functions\despiking_tooblox');
addpath('c:\Users\rps207\Documents\Matlab\Functions\cm_and_cb_utilities');
addpath('c:\Users\rps207\Documents\Matlab\Functions\tsplot');

%loop through Licor files and 
path = 'C:\Users\rps207\Documents\Seasonal Study Data\Quest Underway\';%identify file path
cd(path); %change current directory to folder with files want to load
fList = ls('*.txt');  % whilst in files folder, generates list of LogData files in directory , ls('*LogData*.lvm');     

mfileDir = 'C:\Users\rps207\Documents\Matlab\CO2 NSOP output analysis\'; %path for main matlab analysis
cd(mfileDir);   %change currect directory original matlab folder

Latitude=[];Longitude=[];Windspeedabs=[];Winddirabs=[];AtmTemp=[];Parport=[];Humidity=[];Pressure=[];SeaTemp=[];Salinity=[];O2saturation=[];Atten=[];Chla=[];Cdom=[];Turbidity=[];Und_DT=[];
[r,c]=size(fList);
    for m = 1:r;
        filename= (['C:\Users\rps207\Documents\Seasonal Study Data\Quest Underway\' fList(m,:) ])
        dataimp=dlmread(filename,',',1, 2); %open licor text file
        %import date and time seperately
        delimiter = ','; formatSpec = '%s%s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]'; fileID = fopen(filename,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false); fclose(fileID);
        DMYdataarray = dataArray{:, 1}; %define licor date time
        DMY=DMYdataarray(2:end);
        HMSdataarray = dataArray{:, 2}; %define licor date time
        HMS=HMSdataarray(2:end);
        DMY_datenum=datenum(DMY,'dd/mm/yyyy')+datenum(2000,0,0,0,0,0);
        HMS_datenum=datenum(HMS,'HH:MM:SS')-datenum(2019,0,1,0,0,0);
        dt=HMS_datenum+DMY_datenum; %add 1 hour to get into UTC
        dt=dt+(1/24);
        % Allocate imported array to column variable names
        lat = dataimp(:, 1);
        long = dataimp(:, 2);
        windspdx = dataimp(:, 3);%knots
        windspd=0.514444*windspdx;%conversion to m/s
        
        winddir = dataimp(:, 4);%degrees
        atemp = dataimp(:, 5); %degrees
        hum = dataimp(:, 6); %percent
        pre = dataimp(:, 7); %millibars
        %8 inlet??
        stemp = dataimp(:, 9);
        sal = dataimp(:, 10); %PSU
        o2sat = dataimp(:, 11);
        att = dataimp(:, 12); %cm/s
        par = dataimp(:, 14);%W M^2%par raw13
        chlorophylla = dataimp(:, 16);%Mg m^3 %chla raw 15
        cdom = dataimp(:, 18);%Mg m^3 %cdom raw 17
        tub= dataimp(:, 20);%turbidity raw 19
        
        %concatenate individual files
        Latitude = [Latitude ; lat];%
        Longitude = [Longitude ; long];%
        Windspeedabs = [Windspeedabs ; windspd];%
        Winddirabs = [Winddirabs ; winddir];%
        AtmTemp = [AtmTemp ; atemp];%
        Humidity = [Humidity ; hum];%
        Pressure = [Pressure ; pre];%
        SeaTemp = [SeaTemp ; stemp];%
        Salinity = [Salinity ; sal];%
        O2saturation = [O2saturation ; o2sat];
        Atten = [Atten ; att];
        Parport = [Parport ; par];
        Chla = [Chla ; chlorophylla];
        Cdom = [Cdom ; cdom];
        Turbidity = [Turbidity ; tub];
        Und_DT = [Und_DT ; dt];
             
    end

    
    
 [C,ia,ic]=unique(Und_DT);
 Und_DT=Und_DT(ia);
 
 
        Latitude = Latitude(ia);
        Longitude = Longitude(ia);
        Windspeedabs = Windspeedabs(ia);
        Winddirabs = Winddirabs(ia);
        AtmTemp = AtmTemp(ia);
        Humidity = Humidity(ia);
        Pressure = Pressure(ia);
        SeaTemp = SeaTemp(ia);
        Salinity = Salinity(ia);
        O2saturation = O2saturation(ia);
        Atten = Atten(ia);
        Parport = Parport(ia);
        Chla = Chla(ia);
        Cdom = Cdom(ia);
        Turbidity = Turbidity(ia);
     
    
    
    
% Clear temporary variables
clearvars windspdx C ia ic dataimp filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;
clearvars atemp att c cdom chlorophylla dt fList hum r m lat long  o2sat  par pre sal stemp tub winddir windspd HMSdataarray HMS_datenum DMY_datenum DMYdataarray DMY  DMY_datenumYMD HMS AA BB AAA AAAA ABC ABCD filename delimiter startRow formatSpec fileID dataArray ans numericData rawNumericColumns rawCellColumns R;





% Clear temporary variables
clearvars underwayYear1 filename delimiter startRow endRow formatSpec fileID block dataArrayBlock col dataArray ans path mfileDir;






%%%%%%%%%%%%%%%%%%%%%%%%%%% Save ouput as loadable mat files
 
    if ~exist('Data', 'dir')
        mkdir('Data');
    end
   save([pwd '/Data/Seasonalstudyunderway.mat']);   
    



   



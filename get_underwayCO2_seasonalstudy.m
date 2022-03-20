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

[r,c]=size(fList);

%load Underway CO2
filename = 'C:\Users\rps207\Documents\Seasonal Study Data\Underway CO2\WCO_2016_data_for_MY_RS.csv';
delimiter = ',';
startRow = 2;
endRow = 455;
formatSpec = '%f%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'HeaderLines', startRow-1, 'ReturnOnError', false);
fclose(fileID);
underwayDOY = dataArray{:, 1};
underway_DT = dataArray{:, 2};
underway_DT=datenum(underway_DT(1:454),'dd/mm/yyyy HH:MM');
underwayYear1 = dataArray{:, 3};
underwayLat = dataArray{:, 4};
underwayLong = dataArray{:, 5};
underwayP_atm = dataArray{:, 6};
underwayT_sea = dataArray{:, 7};
underwayS_sea = dataArray{:, 8};
underwaypCO2_sw = dataArray{:, 9};
underwaypCO2_atm = dataArray{:, 10};
underwayfCO2_sw = dataArray{:, 11};
underwayfCO2_atm = dataArray{:, 12};
underwayDpCO2 = dataArray{:, 13};
underwayDfCO2 = dataArray{:, 14};
underwaychlorophyll = dataArray{:, 15};
underwayRel_Wind = dataArray{:, 16};
underwayWind_Dir = dataArray{:, 17};
underwayO2_optode = dataArray{:, 18};


 [C,ia,ic]=unique(underway_DT);

underway_DT=underway_DT(ia);
underwayDOY = underwayDOY(ia);
underwayYear1 = underwayYear1(ia);
underwayLat = underwayLat(ia);
underwayLong = underwayLong(ia);
underwayP_atm = underwayP_atm(ia);
underwayT_sea = underwayT_sea(ia);
underwayS_sea = underwayS_sea(ia);
underwaypCO2_sw = underwaypCO2_sw(ia);
underwaypCO2_atm = underwaypCO2_atm(ia);
underwayfCO2_sw = underwayfCO2_sw(ia);
underwayfCO2_atm = underwayfCO2_atm(ia);
underwayDpCO2 = underwayDpCO2(ia);
underwayDfCO2 = underwayDfCO2(ia);
underwaychlorophyll = underwaychlorophyll(ia);
underwayRel_Wind = underwayRel_Wind(ia);
underwayWind_Dir = underwayWind_Dir(ia);
underwayO2_optode = underwayO2_optode(ia);


underway_DT=underway_DT+(1/24); % convert from GMT to british summer time
underway_DT=underway_DT-(2220/(24*3600)); %account for apparent lag


% Clear temporary variables
clearvars dataimp filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;
clearvars min23 atemp att c cdom chlorophylla dt fList hum r m lat long  o2sat  par pre sal stemp tub winddir windspd HMSdataarray HMS_datenum DMY_datenum DMYdataarray DMY  DMY_datenumYMD HMS AA BB AAA AAAA ABC ABCD filename delimiter startRow formatSpec fileID dataArray ans numericData rawNumericColumns rawCellColumns R;





% Clear temporary variables
clearvars underwayYear1 filename delimiter startRow endRow formatSpec fileID block dataArrayBlock col dataArray ans path mfileDir;

clearvars C ia ic

%%%%%%%%%%%%%%%%%%%%%%%%%%% Save ouput as loadable mat files
 
    if ~exist('Data', 'dir')
        mkdir('Data');
    end
   save([pwd '/Data/underwaypCO2Seasonalstudy.mat']);   
    


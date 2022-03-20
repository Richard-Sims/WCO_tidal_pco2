%get underway seasonal study interp

%get underway seasonal study interp
clc;  clear all ; close all; %reset workspace
tstart = tic;
pFigs = true;
addpath('c:\Users\rps207\Documents\Matlab\Functions');
addpath('c:\Users\rps207\Documents\Matlab\Functions\Add_Axis');
mfileDir = 'C:\Users\rps207\Documents\Matlab\CO2 NSOP output analysis\'; %path for main matlab analysis

%loop through Licor files and 
path = 'C:\Users\rps207\Documents\Seasonal Study Data\LicorData\';%identify file path
cd(path); %change current directory to folder with files want to load
fList = ls('*.txt');  % whilst in files folder, generates list of LogData files in directory , ls('*LogData*.lvm');     
iRow=[];
for i=1:size(fList,1);    iRow = [iRow ; i];  end;
%   LINE TO EXCLUDE PROFILES for i=1:size(fList,1);   if strfind(fList(i,:),'2016-06-10.txt') > 0; elseif strfind(fList(i,:),'2016-05-26.txt') > 0;else iRow = [iRow ; i]; end;    end;
cd(mfileDir);   %change currect directory original matlab folder
fList = fList(iRow,:);

%work through profiles

for j = 1:length(fList(:,1));               % loops through each file in logdata list
    %%%%%%%%%%%%%%%%%%%   Licor data import %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %import Licor text docs into workspace be selecting folder location
    path = 'C:\Users\rps207\Documents\Seasonal Study Data\LicorData\';%identify file path
    fullpathLicor =[path fList(j,1:10) '.txt'];
    [LicorH20B,LicorPressure,LicorTemperature,Licor_Datetime,Licor_Datetime_doy,LicorxCO2B ] = getLicorSeasonalStudy(fullpathLicor,j,fList);
    
    C=load ('C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\Seasonalstudyunderway.mat')
    load ('C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\Seasonalstudyunderway.mat','Und_DT')
    D=load('C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\underwaypCO2Seasonalstudy.mat');
    load('C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\underwaypCO2Seasonalstudy.mat','underway_DT');
    
    
    
    names=fieldnames(C)
    names2=fieldnames(D);

    %create new names for variables with _interp suffix
 for t=1:numel(names)  
    rty=[names(t) '_interp'];x=horzcat(rty{:}) ;new_names(:,t)=cellstr(x);
 end
    % Interpolate variables in a loop using new variable names
   for i=1:numel(names)  
       interpvars.(new_names{i})=interp1(Und_DT,C.(names{i}),Licor_Datetime);
   end
    
       for t=1:numel(names2)  
    rty=[names2(t) '_interp'];x=horzcat(rty{:}) ;new_names2(:,t)=cellstr(x);
 end
    % Interpolate variables in a loop using new variable names
   for i=1:numel(names2)  
       interpvars.(new_names2{i})=interp1(underway_DT,D.(names2{i})',Licor_Datetime);
   end
      
   %%%%%%%%%%%%%%%%%%%%%%%%
   %create variables same longer than  eg 24 hours before and 4 horus after licor file
   
       %create new names for variables with _interp suffix
 for t=1:numel(names)  
    rty=[names(t) '_interplongtime'];x=horzcat(rty{:}) ;new_names(:,t)=cellstr(x);
 end
    % Interpolate variables in a loop using new variable names
   for i=1:numel(names)      
       Time_interp=Licor_Datetime(1)-(1):(1/(24*60*60)):Licor_Datetime(end)+(6/24);
       interpvars.(new_names{i})=interp1(Und_DT,C.(names{i}),Time_interp);
   end
       
    for t=1:numel(names2)  
    rty=[names2(t) '_interplongtime'];x=horzcat(rty{:}) ;new_names2(:,t)=cellstr(x);
 end
    % Interpolate variables in a loop using new variable names
   for i=1:numel(names2)      
       Time_interp=Licor_Datetime(1)-(1):(1/(24*60*60)):Licor_Datetime(end)+(6/24);
       interpvars.(new_names2{i})=interp1(underway_DT,D.(names2{i})',Time_interp);
   end
       

   v2struct(interpvars)

%    clearvars s interpvars LicorH20B LicorPressure s LicorTemperature Licor_Datetime Licor_Datetime_doy LicorxCO2B
%          clearvars  C N S3 Und-DT and car_new  fullpathLicor i irow  mfiledir names new_names Und_DT ans iRow mfileDir path newnames pFigs ptah result rty t tstart u x y
     if ~exist('Data/Underwayinterp', 'dir')
        mkdir('Data/Underwayinterp');
    end
    save([pwd '/Data/Underwayinterp/' fList(j,1:10) '.mat']);
end

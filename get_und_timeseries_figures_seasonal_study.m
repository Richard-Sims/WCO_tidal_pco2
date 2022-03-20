%get mooring data l4
clc;  clear all ; close all; %reset workspace
tstart = tic;
addpath('c:\Users\rps207\Documents\Matlab\Functions');
addpath('c:\Users\rps207\Documents\Matlab\Functions\Add_Axis');
addpath('c:\Users\rps207\Documents\Matlab\Functions\cbdate');
addpath('c:\Users\rps207\Documents\Matlab\Functions\mixing_library');
addpath('c:\Users\rps207\Documents\Matlab\Functions\despiking_tooblox');
addpath('c:\Users\rps207\Documents\Matlab\Functions\cm_and_cb_utilities');
addpath('c:\Users\rps207\Documents\Matlab\Functions\tsplot');
mfileDir = 'C:\Users\rps207\Documents\Matlab\CO2 NSOP output analysis\'; %path for main matlab analysis

degree_symbol= sprintf('%c', char(176));
micro_symbol= sprintf('%c', char(0181));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  LOAD IN L4 MOORING DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%



%find all subdirectories
d = dir('C:\Users\rps207\Documents\Seasonal Study Data\Buoy data\L4\2016');
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
nameFolds=nameFolds(3:end)
path= 'C:\Users\rps207\Documents\Seasonal Study Data\Buoy data\L4\2016\';
cd (path)

L4SST=[];L4Salinity=[];L4Oxygen=[];L4Flourescence=[];L4Turbidity=[];L4PHL=[];L4_dt=[];L4Nit=[]; L4Windspeed=[];L4Winddir=[]; L4Humidity=[];L4Temperature=[];L4Pressure=[];L4Latitude=[]; L4Longitude=[];
L4_dt_latlon=[  ];L4_dt_Nit=[];L4_dt_PHL=[]; L4PAR=[];L4_dt_PAR=[];

for j=1:length(nameFolds)
    folderpath=nameFolds(j);
    string = folderpath{1};
    cd(path)
    cd(string)
    fList = ls('*.txt');  % whilst in files folder, generates list of LogData files in directory , ls('*LogData*.lvm');     
    if isempty(fList) 
        continue
    else
    txtfile=fList(1,:);

    %% Initialize variables.
    path = 'C:\Users\rps207\Documents\Seasonal Study Data\Buoy data\L4\2016\';%identify file path
    filename =[path string '\' txtfile];delimiter = ' ';
    startRow = 10;
    formatSpec = '%*s%f%f%s%s%s%s%*s%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    fclose(fileID);
    col1 = dataArray{:, 1};
    col2 = dataArray{:, 2};
    col3 = dataArray{:, 3};
    col4 = dataArray{:, 4};
    col5 = dataArray{:, 5};
    col6 = dataArray{:, 6};
    
    delimiter = ' ';
    startRow = 10;
    formatSpec = '%s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    fclose(fileID);
    Reading = dataArray{:, 1};

    %%%%%%%%%%%%%%%%%%%%
    [ind,~]=find(ismember(Reading,'GPS'));
    Latitude=col2(ind);
    Longitude=cellfun(@str2num,col3(ind));
    %get time for data
    p=str2num(string(5:end));
    q=str2num(string(1:4));
    dmy=doy2date(p,q);
    dt_latlon=((0:1:(length(ind)-1))/(length(ind)-1))+dmy;

    
    %
    [ind,~]=find(ismember(Reading,'MPB'));
    Pressure=col2(ind);
    Temperature=cellfun(@str2num,col3(ind));
    Humidity=cellfun(@str2num,col4(ind));
    Winddir=cellfun(@str2num,col5(ind));
    Windspeed=cellfun(@str2num,col6(ind));
    
    %
    [ind,~]=find(ismember(Reading,'NLC'));
    Nit=col2(ind);
    dt_Nit=((0:1:(length(ind)-1))/(length(ind)-1))+dmy;

    %
    [ind,~]=find(ismember(Reading,'PHL'));
    PHL=col2(ind);
    %get time for data
    dt_PHL=((0:1:(length(ind)-1))/(length(ind)-1))+dmy;
    
        %
    [ind,~]=find(ismember(Reading,'PAR'));
    par=col2(ind);
    %get time for data
    dt_par=((0:1:(length(ind)-1))/(length(ind)-1))+dmy;
    
    
    %
    [ind,~]=find(ismember(Reading,'WQM'));
    SST=col2(ind);
    Salinity=cellfun(@str2num,col3(ind));
    Oxygen=cellfun(@str2num,col4(ind));
    Flourescence=cellfun(@str2num,col5(ind));
    Turbidity=cellfun(@str2num,col6(ind));
    
    %get time for data
    p=str2num(string(5:end));
    q=str2num(string(1:4));
    dmy=doy2date(p,q);
    dt=((0:1:(length(ind)-1))/(length(ind)-1))+dmy;

    L4SST=[L4SST;SST];
    L4Salinity=[L4Salinity;Salinity];
    L4Oxygen=[L4Oxygen;Oxygen];
    L4Flourescence=[L4Flourescence;Flourescence];
    L4Turbidity=[L4Turbidity;Turbidity];
    L4PHL=[L4PHL;PHL];
    L4_dt=[L4_dt;dt'];
    L4PAR=[L4PAR,par'];
    L4_dt_PAR=[L4_dt_PAR,dt_par];
    L4_dt_latlon=[L4_dt_latlon;dt_latlon'];
    L4_dt_Nit=[L4_dt_Nit;dt_Nit'];
    L4_dt_PHL=[L4_dt_PHL;dt_PHL'];
    L4Nit=[L4Nit;Nit];
    L4Windspeed=[L4Windspeed;Windspeed];
    L4Winddir=[L4Winddir; Winddir]; 
    L4Humidity=[L4Humidity;Humidity];
    L4Temperature=[L4Temperature;Temperature];
    L4Pressure=[L4Pressure;Pressure];
    L4Latitude=[L4Latitude;Latitude];
    L4Longitude=[L4Longitude;Longitude];
    
    clearvars filename delimiter startRow formatSpec fileID dataArray ans;
    clearvars filename delimiter startRow formatSpec fileID dataArray ans;
    end

end
 
clearvars par Salinity Turbidity ans Flourescence Humidity dt dmy fList folderpath ind isub j nameFolds p q string tstart txtfile Latitude Longitude Nit Oxygen PHL Pressure Reading SST Salinity Temperature Turbulence Winddir Windspeed col1 col2 col3 col4 col5 col6 d
%%%%%%%%%%%%%%%%%%%%%%%%%%% Save ouput as loadable mat files
cd(mfileDir)
if ~exist('Data', 'dir')
    mkdir('Data');
end
save([pwd '/Data/Seasonal_study_L4_mooring.mat']);
   

%load in Penlee data
filename = 'C:\Users\rps207\Documents\Seasonal Study Data\Penlee Wind\Penlee_Gas_Met_Data_hr_all.txt';delimiter = '\t';startRow = 2;
formatSpec = '%{MM/dd/yy HH:mm:ss}D%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

Penlee_MonthDayYear_TimeUTC = dataArray{:, 1};
Penlee_dt=datenum(Penlee_MonthDayYear_TimeUTC);
Penlee_Wspd_m_s = dataArray{:, 2};
Penlee_Wdir_deg = dataArray{:, 3};
Penlee_Pres_mb = dataArray{:, 4};
Penlee_RH_percent = dataArray{:, 5};
Penlee_Tair_C = dataArray{:, 6};
Penlee_Dewpt_C = dataArray{:, 7};
Penlee_RainRate_mm_hr = dataArray{:, 8};
Penlee_SO2_ppb = dataArray{:, 9};
Penlee_O3_ppb = dataArray{:, 10};
Penlee_CO2_ppm = dataArray{:, 11};
Penlee_CH4_ppm = dataArray{:, 12};

clearvars filename delimiter startRow formatSpec fileID dataArray ans Penlee_MonthDayYear_TimeUTC;

%load in rainfall rate data
filename = 'C:\Users\rps207\Documents\Seasonal Study Data\riverflow_gdf.csv';
delimiter = ',';
startRow = 20;
formatSpec = '%q%q%q%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
tamar_dt = dataArray{:, 1};
tamar_flowrate_ms3 = dataArray{:, 2};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;



degree_symbol= sprintf('%c', char(176));
micro_symbol= sprintf('%c', char(0181));
%%% Custom RGB colour vectors
colour_teal = [18 150 155] ./ 255;
colour_lightgreen = [94 250 81] ./ 255;
colour_green = [12 195 82] ./ 255;
colour_lightblue = [8 180 238] ./ 255;
colour_darkblue = [1 17 181] ./ 255;
colour_yellow = [251 250 48] ./ 255;
colour_peach = [251 111 66] ./ 255;
colour_crimsom = [192 11 52] ./ 255;
colour_greyshade= [212 229 208] ./ 255;

   
%Reinitalise the filelist as required for plotting
path = 'C:\Users\rps207\Documents\Matlab\CO2 NSOP output analysis\Data\';%identify file path
cd(path); %change current directory to folder with files want to load
fList = ls('*.mat');  % whilst in files folder, generates list of LogData files in directory , ls('*LogData*.lvm');     
iRow=[];
for i=1:size(fList,1);   if strfind(fList(i,:),'underway.txt') > 0; else iRow = [iRow ; i]; end;    end;
cd(mfileDir);   %change currect directory original matlab folder
fList = fList(iRow,:);
fList=fList(30:45,1:10);
[rows cols]=size(fList);
datenumdeploy=datenum(fList,'yyyy-mm-dd')


%load in wind from Penlee






load ('C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\Seasonalstudyunderway.mat');
load ('C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\underwaypCO2Seasonalstudy.mat');


deg = km2deg(1);
[lon,lat] = scircle1(-4.217,50.25,deg);
[in ~]= inpolygon(-1*Longitude,Latitude,lon,lat);
L4ind=find(in==1);

deg = km2deg(1);
[in2 ~]= inpolygon(underwayLong,underwayLat,lon,lat);
L4indco2=find(in2==1);


% plot(-1*Longitude,Latitude,'r*')
% hold on
% plot(lon,lat)


%Wind and Par
figure(1)  
subplot(3,1,1)
set(gcf, 'Color','w','Position', get(0,'Screensize'));
for i=1:rows;
    path = 'C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\';
    fListfilename = [path fList(i,1:10) '.mat'];
    load(fListfilename,'inds','indt','Licor_Datetime')
    Licor_Datetime(inds(1))
    r=indt(end);
    x=[r inds(1) inds(1) r];
    y= [ max([L4Windspeed;Windspeedabs]) max([L4Windspeed;Windspeedabs]) 0 0];
    patch(Licor_Datetime(x),y,colour_greyshade); 
    clearvars indslag indtlag Licor_Datetime

end
hold on;
despiked_Windspeedabs=func_despike_phasespace3d(Windspeedabs);
despiked_Windspeedabs=func_despike_phasespace3d(L4Windspeed);
plot(Penlee_dt,Penlee_Wspd_m_s,'-g','LineWidth',2);
plot(L4_dt,L4Windspeed,'b','LineWidth',2);
plot(Und_DT(L4ind),Windspeedabs(L4ind),'*r','LineWidth',2);
dynamicDateTicks([], [], 'mm/dd');
xlabel('Time (month/day)','FontSize',16);
ylabel({[' Wind speed {(ms^{-1})}']},'FontSize',16);
ylim([0 max(despiked_Windspeedabs)])
klp=datenum(2016,04,15,0,0,1)
jkl=datenum(2016,09,30,0,0,1)
xlim([klp jkl])
set(gca,'FontSize',16);
set(gca,'FontSize',16);
text(0.02,0.93,'(A)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized'); 


subplot(3,1,2)
for i=1:rows;
    path = 'C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\';
    fListfilename = [path fList(i,1:10) '.mat'];
    load(fListfilename,'inds','indt','Licor_Datetime')
    Licor_Datetime(inds(1))
    r=indt(end);
    x=[r inds(1) inds(1) r];
    y= [ max([Winddirabs;L4Winddir]) max([Winddirabs;L4Winddir]) 0 0];
    patch(Licor_Datetime(x),y,colour_greyshade); 
    clearvars indslag indtlag Licor_Datetime
hold on
end
hold on;
plot(Penlee_dt,Penlee_Wdir_deg,'-g','LineWidth',2);
plot(L4_dt,L4Winddir,'b','LineWidth',2);
plot(Und_DT(L4ind),Winddirabs(L4ind),'*r','LineWidth',2);
dynamicDateTicks([], [], 'mm/dd');
xlabel('Time (month/day)','FontSize',16);
ylabel({[' Wind direction(',num2str(degree_symbol) ')']},'FontSize',16);
ylim([0 max(L4Winddir)])
klp=datenum(2016,04,15,0,0,1)
jkl=datenum(2016,09,30,0,0,1)
xlim([klp jkl])
set(gca,'FontSize',16);
set(gca,'FontSize',16);
text(0.02,0.93,'(B)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized'); 

subplot(3,1,3)
for i=1:rows;
    path = 'C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\';
    fListfilename = [path fList(i,1:10) '.mat'];
    load(fListfilename,'inds','indt','Licor_Datetime')
    Licor_Datetime(inds(1))
    r=indt(end);
    x=[r inds(1) inds(1) r];
    y= [ max(L4PAR) max(L4PAR) 0 0];
    patch(Licor_Datetime(x),y,colour_greyshade); 
    clearvars indslag indtlag Licor_Datetime
    hold on
end
plot(L4_dt_PAR,L4PAR,'b','LineWidth',2);
hold on
plot(Und_DT(L4ind),Parport(L4ind),'*r','LineWidth',2);
dynamicDateTicks([], [], 'mm/dd');
xlabel('Time (month/day)','FontSize',16);
ylabel({['PAR (Wm^{-2})']},'FontSize',16);
ylim([0 max(L4PAR)])
klp=datenum(2016,04,15,0,0,1)
jkl=datenum(2016,09,30,0,0,1)
xlim([klp jkl])
set(gca,'FontSize',16);
set(gca,'FontSize',16);
text(0.02,0.93,'(C)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized'); 





%%%%%%%%%%%%%%%%%Temperature and salinity
figure(2)    
subplot(2,1,1)
for i=1:rows;
    path = 'C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\';
    fListfilename = [path fList(i,1:10) '.mat'];
    load(fListfilename,'inds','indt','Licor_Datetime')
    Licor_Datetime(inds(1))
    r=indt(end);
    x=[r inds(1) inds(1) r];
    y= [ max(L4SST) max(L4SST) 0 0];
    patch(Licor_Datetime(x),y,colour_greyshade); 
    clearvars indslag indtlag Licor_Datetime

end
hold on;
plot(L4_dt,L4SST,'LineWidth',2);
plot(Und_DT(L4ind),SeaTemp(L4ind),'*r','LineWidth',2);
dynamicDateTicks([], [], 'mm/dd');
ylabel({['Temperature (',num2str(degree_symbol),'C)']},'FontSize',22);
ylim([9.5 max(L4SST)])
klp=datenum(2016,04,15,0,0,1)
jkl=datenum(2016,09,30,0,0,1)
xlim([klp jkl])
set(gca,'FontSize',22);
set(gca,'FontSize',22);
text(0.02,0.93,'(A)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized');
xlabel('Time (month/day)','FontSize',22);

subplot(2,1,2)
for i=1:rows;
    path = 'C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\';
    fListfilename = [path fList(i,1:10) '.mat'];
    load(fListfilename,'inds','indt','Licor_Datetime')
    Licor_Datetime(inds(1))
    r=indt(end);
    x=[r inds(1) inds(1) r];
    y= [ max(L4Salinity) max(L4Salinity) 0 0];
    patch(Licor_Datetime(x),y,colour_greyshade); 
    clearvars indslag indtlag Licor_Datetime

end
hold on;
plot(L4_dt,L4Salinity,'LineWidth',2);
plot(Und_DT(L4ind),Salinity(L4ind),'*r','LineWidth',2);
dynamicDateTicks([], [], 'mm/dd');
xlabel('Time (month/day)','FontSize',22);
ylabel({['Salinity(PSU)']},'FontSize',22);
ylim([34 max(L4Salinity)])
klp=datenum(2016,04,15,0,0,1)
jkl=datenum(2016,09,30,0,0,1)
xlim([klp jkl])
set(gca,'FontSize',22);
set(gca,'FontSize',22);
text(0.02,0.93,'(B)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% oxygen and nitrate
figure(3)
subplot(2,1,1)
for i=1:rows;
    path = 'C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\';
    fListfilename = [path fList(i,1:10) '.mat'];
    load(fListfilename,'inds','indt','Licor_Datetime')
    Licor_Datetime(inds(1))
    r=indt(end);
    x=[r inds(1) inds(1) r];
    y= [ max(L4Oxygen) max(L4Oxygen) 0 0];
    patch(Licor_Datetime(x),y,colour_greyshade); 
    clearvars indslag indtlag Licor_Datetime

end
hold on;
plot(Und_DT(L4ind),O2saturation(L4ind),'*r','LineWidth',2);
plot(L4_dt,L4Oxygen,'LineWidth',2);
dynamicDateTicks([], [], 'mm/dd');
ylabel(['Oxygen' char(10) '(',num2str(micro_symbol),'m)'],'FontSize',22);ylim([210 310])
klp=datenum(2016,04,15,0,0,1)
jkl=datenum(2016,09,30,0,0,1)
xlim([klp jkl])
set(gca,'FontSize',22);
set(gca,'FontSize',22);
text(0.02,0.95,'(A)','color','k','Fontsize',22,'Fontweight','bold','units','normalized'); 
subplot(2,1,2)
plot(L4_dt_Nit,L4Nit,'LineWidth',2);
for i=1:rows;
    path = 'C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\';
    fListfilename = [path fList(i,1:10) '.mat'];
    load(fListfilename,'inds','indt','Licor_Datetime')
    Licor_Datetime(inds(1))
    r=indt(end);
    x=[r inds(1) inds(1) r];
    y= [ max(L4Nit) max(L4Nit) 0 0];
    patch(Licor_Datetime(x),y,colour_greyshade); 
    clearvars indslag indtlag Licor_Datetime

end
hold on;
text(0.02,0.95,'(B)','color','k','Fontsize',22,'Fontweight','bold','units','normalized'); 
dynamicDateTicks([], [], 'mm/dd');
xlabel('Time (month/day)','FontSize',22);
ylabel(['Nitrate' char(10) '(',num2str(micro_symbol),'m)'],'FontSize',22);ylim([6 18.5])
klp=datenum(2016,04,15,0,0,1)
jkl=datenum(2016,09,30,0,0,1)
xlim([klp jkl])
set(gca,'FontSize',22);
set(gca,'FontSize',22);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% flourescnence and turbidity
figure(4)
subplot(2,1,1)
for i=1:rows;
    path = 'C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\';
    fListfilename = [path fList(i,1:10) '.mat'];
    load(fListfilename,'inds','indt','Licor_Datetime')
    Licor_Datetime(inds(1))
    r=indt(end);
    x=[r inds(1) inds(1) r];
    y= [ max(L4Flourescence) max(L4Flourescence) 0 0];
    patch(Licor_Datetime(x),y,colour_greyshade); 
    clearvars indslag indtlag Licor_Datetime

end
hold on;
plot(L4_dt,L4Flourescence,'LineWidth',2);
plot(Und_DT(L4ind),Chla(L4ind),'*r','LineWidth',2);
dynamicDateTicks([], [], 'mm/dd');
ylabel({['Flourescence']},'FontSize',22);
ylim([0 max(L4Flourescence)])
klp=datenum(2016,04,15,0,0,1)
jkl=datenum(2016,09,30,0,0,1)
xlim([klp jkl])
set(gca,'FontSize',22);
set(gca,'FontSize',22);
text(0.02,0.93,'(A)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
subplot(2,1,2)
for i=1:rows;
    path = 'C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\';
    fListfilename = [path fList(i,1:10) '.mat'];
    load(fListfilename,'inds','indt','Licor_Datetime')
    Licor_Datetime(inds(1))
    r=indt(end);
    x=[r inds(1) inds(1) r];
    y= [ max(L4Turbidity) max(L4Turbidity) 0 0];
    patch(Licor_Datetime(x),y,colour_greyshade); 
    clearvars indslag indtlag Licor_Datetime

end
hold on;
plot(L4_dt,L4Turbidity,'LineWidth',2);
plot(Und_DT(L4ind),Turbidity(L4ind),'*r','LineWidth',2);
dynamicDateTicks([], [], 'mm/dd');
ylabel({['Turbidity']},'FontSize',22);
ylim([0 max(L4Turbidity)])
klp=datenum(2016,04,15,0,0,1)
jkl=datenum(2016,09,30,0,0,1)
xlim([klp jkl])
set(gca,'FontSize',22);
set(gca,'FontSize',22);
text(0.02,0.93,'(B)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized'); 

tamar_dtt=datenum(tamar_dt);
tamar_flowrate_mms3=str2double(tamar_flowrate_ms3);

[p o]=find(tamar_dtt>datenum('2016-04-01 00:00:00','yyyy-mm-dd HH:MM:SS') & tamar_dtt<datenum('2016-09-30 00:00:00','yyyy-mm-dd HH:MM:SS'));

figure(5)
plot(tamar_dtt(p),tamar_flowrate_mms3(p),'LineWidth',2);
for i=1:rows;
    path = 'C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\';
    fListfilename = [path fList(i,1:10) '.mat'];
    load(fListfilename,'inds','indt','Licor_Datetime')
    Licor_Datetime(inds(1))
    r=indt(end);
    x=[r inds(1) inds(1) r];
    y= [ max(L4Turbidity) max(L4Turbidity) 0 0];
    patch(Licor_Datetime(x),y,colour_greyshade); 
    clearvars indslag indtlag Licor_Datetime

end
set(gca,'FontSize',22);
set(gca,'FontSize',22);
klp=datenum(2016,04,15,0,0,1)
jkl=datenum(2016,09,30,0,0,1)
ylabel({['River Tamar flowrate (m^{3} s^{-1})']},'FontSize',22);
xlabel({['Time']},'FontSize',22);
dynamicDateTicks([], [], 'dd');


%timeline of deployments
figure(101)
subplot(6,1,1)
plot([datenumdeploy(1:end),datenumdeploy(1:end)],[0,1],'Color',colour_green,'LineWidth',2);
dynamicDateTicks([], [], 'yy:mm');
xlabel('Time','FontSize',22);
set(gca,'FontSize',22);
set(gca,'FontSize',22);
set(gca,'ytick',[])
set(gca,'yticklabel',[])


load('C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\Plymouth_tides.mat')

figure(102)
subplot(3,1,1)
for i=1:rows;
    path = 'C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\';
    fListfilename = [path fList(i,1:10) '.mat'];
    load(fListfilename,'inds','indt','Licor_Datetime')
    Licor_Datetime(inds(1))
    r=indt(end);
    x=[r inds(1) inds(1) r];
    y= [ max(underwayfCO2_sw) max(underwayfCO2_sw) 0 0];
    patch(Licor_Datetime(x),y,colour_greyshade); 
    clearvars indslag indtlag Licor_Datetime
end
hold on;
plot(underway_DT(L4indco2),underwayfCO2_sw(L4indco2),'*r','LineWidth',2);
dynamicDateTicks([], [], 'mm/dd');
ylabel(['fCO_{2(sw)}' char(10) '(',num2str(micro_symbol),'atm)'],'FontSize',22);
ylim([280 max(underwayfCO2_sw)])
klp=datenum(2016,04,15,0,0,1)
jkl=datenum(2016,09,30,0,0,1)
xlim([klp jkl])
set(gca,'FontSize',22);
set(gca,'FontSize',22);
text(0.02,0.93,'(A)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
yyaxis right
plot(Plymouth_tidal_height_dt,Plymouth_tidal_height,'-k','Linewidth',0.5)



subplot(3,1,2)
for i=1:rows;
    path = 'C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\';
    fListfilename = [path fList(i,1:10) '.mat'];
    load(fListfilename,'inds','indt','Licor_Datetime')
    Licor_Datetime(inds(1))
    r=indt(end);
    x=[r inds(1) inds(1) r];
    y= [ max(underwayfCO2_atm) max(underwayfCO2_atm) 0 0];
    patch(Licor_Datetime(x),y,colour_greyshade); 
    clearvars indslag indtlag Licor_Datetime
end
hold on;
plot(underway_DT(L4indco2),underwayfCO2_atm(L4indco2),'*r','LineWidth',2);
dynamicDateTicks([], [], 'mm/dd');
ylabel(['fCO_{2(atm)}' char(10) '(',num2str(micro_symbol),'atm)'],'FontSize',22);
ylim([380 440])
klp=datenum(2016,04,15,0,0,1)
jkl=datenum(2016,09,30,0,0,1)
xlim([klp jkl])
set(gca,'FontSize',22);
set(gca,'FontSize',22);
text(0.02,0.93,'(B)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized'); 

subplot(3,1,3)
for i=1:rows;
    path = 'C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\';
    fListfilename = [path fList(i,1:10) '.mat'];
    load(fListfilename,'inds','indt','Licor_Datetime')
    Licor_Datetime(inds(1))
    r=indt(end);
    x=[r inds(1) inds(1) r];
    y= [ max(underwayDfCO2) max(underwayDfCO2) -140 0];
    patch(Licor_Datetime(x),y,colour_greyshade); 
    clearvars indslag indtlag Licor_Datetime
end
hold on;
plot(underway_DT(L4indco2),underwayDfCO2(L4indco2),'*r','LineWidth',2);
dynamicDateTicks([], [], 'mm/dd');
xlabel('Time (month/day)','FontSize',22);
    ylabel(['{\Delta}fCO_{2(sw)}' char(10) '(',num2str(micro_symbol),'atm)'],'FontSize',22);
ylim([-140 40])
klp=datenum(2016,04,15,0,0,1)
jkl=datenum(2016,09,30,0,0,1)
xlim([klp jkl])
set(gca,'FontSize',22);
set(gca,'FontSize',22);
text(0.02,0.93,'(C)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized'); 

% subplot(2,1,2)
% set(gcf, 'Color','w','Position', get(0,'Screensize'));
% despiked_Downirrmax=func_despike_phasespace3d(Downirrmax);
% plot(Und_DT2,despiked_Downirrmax,'LineWidth',2);
% for i=1:rows;
%     path = 'C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\';
%     fListfilename = [path fList(i,1:10) '.mat'];
%     load(fListfilename,'inds','indt','Licor_Datetime')
%     Licor_Datetime(inds(1))
%     r=indt(end);
%     x=[r inds(1) inds(1) r];
%     y= [ max(despiked_Downirrmax) 1000 0 0];
%     patch(Licor_Datetime(x),y,colour_greyshade); 
% end
% hold on;
% dynamicDateTicks([], [], 'mm/dd');
% xlabel('Time (month/day)','FontSize',16);
% ylabel({['Irradiance (Wm^{-2})']},'FontSize',16);
% set(gca,'FontSize',16);
% set(gca,'FontSize',16);
% klp=datenum(2014,10,20,0,0,1)
% jkl=datenum(2014,11,05,0,0,1)
% xlim([klp jkl])
% ylim([min(despiked_Downirrmax) 1000])
% text(0.02,0.93,'(B)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized'); 




% %work out avergae variables
% 
% [ind1, ~]=find(Und_DT>(datenum('2014-10-20 01:00:00','yyyy-mm-dd HH:MM:SS')) & Und_DT<(datenum('2014-11-05 01:00:00','yyyy-mm-dd HH:MM:SS')));
% [ind2, ~]=find(Und_DT2>(datenum('2014-10-20 01:00:00','yyyy-mm-dd HH:MM:SS')) & Und_DT2<(datenum('2014-11-05 01:00:00','yyyy-mm-dd HH:MM:SS')));
% 
% mean(Windspeedabs(ind1))
% std(Windspeedabs(ind1))
% 
% mean(Downirrmax(ind2))
% std(Downirrmax(ind2))


%avergae par
mean(L4PAR(750:end-200))
std(L4PAR(750:end-200))


Penleewind_interp= interp1(Penlee_dt,Penlee_Wspd_m_s ,L4_dt); %interpolate both temperature probes so they are the same lemngth 

[c,v]=find(L4Windspeed>0);

figure(1000)
scatter(Penleewind_interp(c),L4Windspeed(c))

b = polyfit(Penleewind_interp(c), L4Windspeed(c), 1);
f = polyval(b, Penleewind_interp(c));

[r p] = corrcoef(Penleewind_interp(c), L4Windspeed(c));


%figure plots for L4 horizontal co2 paper 

%hoffmuller style plot for co2
%jumbo figure for paper combining all the transects as distance plots
opengl('save', 'software')
clc;  clear all ; close all; %reset workspace
tstart = tic;
addpath('c:\Users\rps207\Documents\Matlab\Functions');
addpath('c:\Users\rps207\Documents\Matlab\Functions\Add_Axis');
addpath('c:\Users\rps207\Documents\Matlab\Functions\cbdate');
addpath('c:\Users\rps207\Documents\Matlab\Functions\mixing_library');
addpath('c:\Users\rps207\Documents\Matlab\Functions\despiking_tooblox');
addpath('c:\Users\rps207\Documents\Matlab\Functions\cm_and_cb_utilities');
addpath('c:\Users\rps207\Documents\Matlab\Functions\m_map');
addpath('C:\Users\rps207\Documents\MATLAB\Functions\Colormaps\Colormaps (5)\Colormaps');
[cm_data]=viridis();
set(groot,'DefaultFigureColormap',cm_data)
mfileDir = 'C:\Users\rps207\Documents\Matlab\CO2 NSOP output analysis\'; %path for main matlab analysis

degree_symbol= sprintf('%c', char(176));
micro_symbol= sprintf('%c', char(0181));
load coast.dat;
markersty='p' ;
colormap(parula)
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
colour_crimson = [192 11 52] ./ 255;
colour_greyshade= [212 229 208] ./ 255;
colour_verylightblue = [204 255 229] ./ 255;
colour_yellowlight = [255 255 153] ./ 255;
colour_peachback = [255 255 204] ./ 255;

%loop through Licor files and 
path = 'C:\Users\rps207\Documents\Seasonal Study Data\LicorData\';%identify file path
cd(path); %change current directory to folder with files want to load
fList = ls('*.txt');  % whilst in files folder, generates list of LogData files in directory , ls('*LogData*.lvm');     
iRow=[];
for i=1:size(fList,1);    iRow = [iRow ; i];  end;
%   LINE TO EXCLUDE PROFILES 
for i=1:size(fList,1);   if strfind(fList(i,:),'2016-06-10.txt') > 0; elseif strfind(fList(i,:),'2016-05-26.txt') > 0;else iRow = [iRow ; i]; end;    end;
cd(mfileDir);   %change currect directory original matlab folder
% fList = fList(iRow,:);

% fList = fList(4,:);

%NOTES ON WHEN TO IGNORE UNDERWAY DATA

%04/27 no underway on , ignore in and out 
%05/11 no underway temp on  either in or out, ignore both journies
%05/26  only recorded way out  large spike suspicious, dont include either
%06/10 not loggging on way out
%06/15 mfc not working on way out
%06/22 bucket response time tests on way back ignore way back
%08/04 water in system on way out scrap
%08/17 water in system co2 values very low on way out scrap
%08/24 underway data cuts out on way back scrap
%09/02 no data on return journey- not clear from notes scrap way back
%09/21 did not log on way out scrap


% in
%  j==1 no underway temp data
% j==2 | no underway
% j==3 | no underway temp data
% j==6 | was running response tests
% j==15 | no underway temp data
% j==16 no underway temp data


%out
% j==1 no underway temp data
% j==2 | no underway
% j==3 | no underway temp data
% j==4; didnt log
% j==5 %mfc not working on way out
% j==12 %water bad agreement
% j==14 %extremely low data must be wrong
% j==16 %no data
% j==18 didnt log successfully









figure(1)
subsamplefig1=50;
ss_start= datenum('2016-06-09 00:00:00','yyyy-mm-dd HH:MM:SS');
ss_end= datenum('2016-09-22 00:00:00','yyyy-mm-dd HH:MM:SS');hold on

subplot(1, 3, 1)  
for j = 1:length(fList(:,1)); 
 %plot all data on journeys back from L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
else    hold on
load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'undind1','undind2','undind3','undind4','Salinity_interp','Licor_Datetime','dist','inds','indt');
%skip if underway indexes are empty
scatter(dist(undind3:subsamplefig1:undind4),Licor_Datetime(undind3:subsamplefig1:undind4),25, Salinity_interp(undind3:subsamplefig1:undind4),'^')
hold off
end
end
 
for j = 1:length(fList(:,1)); 
%plot all data on journeys out to L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')
else    hold on
load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'undind1','undind2','undind3','undind4','Salinity_interp','Licor_Datetime','dist','inds','indt');
%skip if underway indexes are empty
scatter(dist(undind1:subsamplefig1:undind2),Licor_Datetime(undind1:subsamplefig1:undind2),25, Salinity_interp(undind1:subsamplefig1:undind2),'v')
hold off
end
end
set(gca,'fontsize',12)
set(gca,'fontsize',12)
xlabel(['Distance from L4 (km)'],'FontSize',12);
ylabel(['Date'],'FontSize',12);
set(gca, 'YTick', 736482:7:736598)
set(gca, 'XTick', 0:2:12)
datetick('y','mmm/dd', 'keepticks')
u=colorbar;title(u,['Salinity (PSU)'],'FontSize',12);
hold on
plot(zeros(length(34:0.1:35.4)),34:0.1:35.4,'k','LineWidth',2);
plot(10*ones(length(34:0.1:35.4)),34:0.1:35.4,'k','LineWidth',2);
set(gca,'Xlim',[0 12]);
set(gca,'Ylim',[ss_start ss_end]);
text(0.02,1.05,'(a)','color','k','Fontsize',13,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
hold on

subplot(1, 3, 2)  
temptime=[];tempstats=[];
for j = 1:length(fList(:,1)); 
 %plot all data on journeys back from L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
else    hold on    
load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'undind1','undind2','undind3','undind4','SeaTemp_interp','Licor_Datetime','dist','inds','indt');
%crop out 04-27,05-11 and 05-26 udnerway was not working
scatter(dist(undind3:subsamplefig1:undind4),Licor_Datetime(undind3:subsamplefig1:undind4),25, SeaTemp_interp(undind3:subsamplefig1:undind4),'^')
temptime=[temptime;Licor_Datetime(undind3:subsamplefig1:undind4)];
tempstats=[tempstats;SeaTemp_interp(undind3:subsamplefig1:undind4)];
hold off
end
end

for j = 1:length(fList(:,1)); 
%plot all data on journeys out to L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')
else    hold on
load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'undind1','undind2','undind3','undind4','SeaTemp_interp','Licor_Datetime','dist','inds','indt');
%crop out 04-27,05-11 and 05-26 udnerway was not working
scatter(dist(undind1:subsamplefig1:undind2),Licor_Datetime(undind1:subsamplefig1:undind2),25, SeaTemp_interp(undind1:subsamplefig1:undind2),'v')
hold off
end
end

set(gca,'fontsize',12)
set(gca,'fontsize',12)
xlabel(['Distance from L4 (km)'],'FontSize',12);
set(gca, 'YTick', 736482:7:736598)
set(gca, 'XTick', 0:2:12)
datetick('y','mmm/dd', 'keepticks')
u=colorbar;title(u,['Temperature (',num2str(degree_symbol),'C)'],'FontSize',12);
hold on
plot(zeros(length(34:0.1:35.4)),34:0.1:35.4,'k','LineWidth',2);
plot(10*ones(length(34:0.1:35.4)),34:0.1:35.4,'k','LineWidth',2);
set(gca,'Xlim',[0 12]);
set(gca,'Ylim',[ss_start ss_end]);
text(0.02,1.05,'(b)','color','k','Fontsize',13,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
hold on

subplot(1, 3, 3) 
for j = 1:length(fList(:,1)); 
hold on
%plot all data on journeys back from L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
else       
load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'undind1','undind2','undind3','undind4','Licor_fco2_combined','Licor_Datetime','dist','underway_dist','underway_DT','underwayfCO2_sw','underwayLat','underwayLong');
%skip if underway indexes are empty
scatter(dist(undind3:subsamplefig1:undind4),Licor_Datetime((undind3:subsamplefig1:undind4)),25, Licor_fco2_combined(undind3:subsamplefig1:undind4),'^')
 hold on
 %overlay showerhead data
[a ~]=find(underway_DT>Licor_Datetime(1) & underway_DT<Licor_Datetime(end));
scatter(underway_dist(a), underway_DT(a),50,underwayfCO2_sw(a),'d','filled')
end
end
for j = 1:length(fList(:,1)); 
%plot all data on journeys out to L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')
else    
load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'undind1','undind2','undind3','undind4','Licor_fco2_combined','Licor_Datetime','dist','underway_dist','underway_DT','underwayfCO2_sw','underwayLat','underwayLong');
%skip if underway indexes are empty
scatter(dist(undind1:subsamplefig1:undind2),Licor_Datetime(undind1:subsamplefig1:undind2),25, Licor_fco2_combined(undind1:subsamplefig1:undind2),'v')
hold on
%overlay showerhead data
[b ~]=find(underway_DT>Licor_Datetime(1) & underway_DT<Licor_Datetime(end));
scatter(underway_dist(a), underway_DT(a),40,underwayfCO2_sw(a),'d','filled')
end
end
set(gca,'fontsize',12)
set(gca,'fontsize',12)
xlabel(['Distance from L4 (km)'],'FontSize',12);
set(gca, 'YTick', 736482:7:736598)
set(gca, 'XTick', 0:2:12)
datetick('y','mmm/dd', 'keepticks')
u=colorbar;title(u,['fCO_2_(_s_w_)' '(',num2str(micro_symbol),'atm)'],'FontSize',12);
plot(zeros(length(34:0.1:35.4)),34:0.1:35.4,'k','LineWidth',2);
plot(10*ones(length(34:0.1:35.4)),34:0.1:35.4,'k','LineWidth',2);
text(0.02,1.05,'(c)','color','k','Fontsize',13,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
set(gca,'Xlim',[0 12]);
set(gca,'Ylim',[ss_start ss_end]);

NSOP_overlap_=[]; Shower_overlap =[]; NSOP_overlap_und_ =[]; Shower_overlap_und =[];


%plot points where the two systems overlap
figure(800)
for j = 1:length(fList(:,1)); 
    j
    hold on
load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'shower_overlap','NSOP_overlap','und_dist_L4','shower_overlap_und','NSOP_overlap_und','und_dist_L4_und');
scatter(shower_overlap,NSOP_overlap,50,und_dist_L4,'filled')
 %add date label
%     text(dist(undind4),Salinity_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
    hold on
    klp=scatter(shower_overlap_und,NSOP_overlap_und,80,und_dist_L4_und,'s','filled') 
%combine for stats
   Shower_overlap =[Shower_overlap ; shower_overlap];
   NSOP_overlap_ =[NSOP_overlap_ ; NSOP_overlap];
   NSOP_overlap_und_ =[NSOP_overlap_und_ ; NSOP_overlap_und];
   Shower_overlap_und =[Shower_overlap_und ; shower_overlap_und];
end
xlabel(['Showerhead fCO_2(',num2str(micro_symbol),'atm)'],'fontsize',22)
ylabel(['Membrane fCO_2(',num2str(micro_symbol),'atm)'],'fontsize',22)
set(gca,'fontsize',22)
set(gca,'fontsize',22)
set(gca,'Xlim',[320 440]);
set(gca,'Ylim',[320 440]);
ld=320:1:440;
hold on
plot(ld,ld,'-');
b=colorbar
ylabel(b, 'Distance from L4 (km)','fontsize',22);

% 
% %probably for the appendix
% figure(11)
% hold on
% subplot(1, 3, 1)  
% for j = 1:length(fList(:,1)); 
%     hold on
%     load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'undind1','undind2','undind3','undind4','Salinity_interp','Licor_Datetime','dist');
%     scatter(dist(undind3:undind4), Salinity_interp(undind3:undind4),25,Licor_Datetime((undind3:undind4)))
%     %add date label
% %     text(dist(undind4),Salinity_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
%     hold off
% end
% xlabel('Distance from L4 (km)');
% set(gca,'fontsize',22)
% set(gca,'fontsize',22)
% ylabel(['Salinity'],'FontSize',24);
% set(gca,'Xlim',[-1 11]);
% u=colorbar; cbdate(u,'dd/mm');
% title(u,'Date(dd/mm)','FontSize',22);    hold on
% plot(zeros(length(34:0.1:35.4)),34:0.1:35.4,'k','LineWidth',4)
% plot(10*ones(length(34:0.1:35.4)),34:0.1:35.4,'k','LineWidth',4)
%          
%          
% hold on
% subplot(1, 3, 2)  
% for j = 1:length(fList(:,1)); 
%     hold on
% load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'undind1','undind2','undind3','undind4','SeaTemp_interp','Licor_Datetime','dist');
% %skip if underway indexes are empty
% scatter(dist(undind3:undind4), SeaTemp_interp(undind3:undind4),25,Licor_Datetime((undind3:undind4)))
% % text(dist(undind4),SeaTemp_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% hold off
% end
% xlabel('Distance from L4 (km)');
% set(gca,'fontsize',22)
% set(gca,'fontsize',22)
% ylabel(['Temperature (',num2str(degree_symbol),'C)'],'FontSize',24);
% set(gca,'Xlim',[-1 11]);
% u=colorbar; cbdate(u,'dd/mm');
% title(u,'Date(dd/mm)','FontSize',22);    hold on
% plot(zeros(length(13:1:18)),13:1:18,'k','LineWidth',4);
% plot(10*ones(length(13:1:18)),13:1:18,'k','LineWidth',4)
% 
% 
% hold on
% subplot(1, 3, 3)  
% for j = 1:length(fList(:,1)); 
%     hold on
% load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'undind1','undind2','undind3','undind4','Licor_fco2_combined','Licor_Datetime','dist');
% %skip if underway indexes are empty
% scatter(dist(undind3:undind4), Licor_fco2_combined(undind3:undind4),25,Licor_Datetime((undind3:undind4)))
% % text(dist(undind4),SeaTemp_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% hold off
% end
% xlabel('Distance from L4 (km)');
% set(gca,'fontsize',22)
% set(gca,'fontsize',22)
% ylabel(['fCO_2_(_s_w_)' '(',num2str(micro_symbol),'atm)'],'FontSize',24);
% set(gca,'Xlim',[-1 11]);
% u=colorbar; cbdate(u,'dd/mm');
% title(u,'Date(dd/mm)','FontSize',22);    hold on
% plot(zeros(length(250:1:550)),250:1:550,'k','LineWidth',4)
% plot(10*ones(length(250:1:550)),250:1:550,'k','LineWidth',4)

figure(12)
%plot river Tamar bottle data
load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\Tamar_bottles.mat'],'Tamar_Date','Tamar_DICumolkg','Tamar_fco2_calc_ta_dic','Tamar_LatN','Tamar_LongW','Tamar_Salinity','Tamar_TAumolkg','Tamar_Temp1','Salinity_interp','Licor_Datetime','dist');
%on october 1st salinity was 35.15, march 11 35.18 , this is 12 days
%after(best i can do though)unless i can get it from quest underway?. 
subplot(1,2,2)
        Tamarlabels={'T1' 'T2' 'T3' 'T4' 'T5' 'T6' 'T7' 'T9' 'T8'}
        Tamar_distL4=[];
     for t=1:length(Tamar_LongW)
        Tamar_distL4(t,:)= pos2dist(50.25, -4.217,Tamar_LatN(t),-Tamar_LongW(t),2);
     end
     scatter(Tamar_fco2_calc_ta_dic([1:6])-Tamar_fco2_calc_ta_dic(7),Tamar_Salinity([1:6])-Tamar_Salinity(7),30,Tamar_distL4([1:6]),'d','filled')
hold on
        
    labs=text(Tamar_fco2_calc_ta_dic([1:6])-Tamar_fco2_calc_ta_dic(7)+50,Tamar_Salinity([1:6])-Tamar_Salinity(7), Tamarlabels([1:6]),'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
     set(labs,'BackgroundColor', 'none');
    xlabel(['\Delta fCO_2(',num2str(micro_symbol),'atm)'],'fontsize',22);
    set(gca,'fontsize',22);
    set(gca,'fontsize',22);
    ylabel('\Delta Salinity','fontsize',22);
    hj=colorbar;
    ylabel(hj,'Distance from L4 (km)','FontSize',16);
    set(gca,'fontsize',16);
    set(gca,'fontsize',16);
    hold on
    y=0:-0.1:-25;
    x=(y*-39.83)+ 5.50;
    plot(x,y)
    xlim([ 0 500])
subplot(1,2,1)
subsample=30;
hold on
for j = 1:length(fList(:,1));
        %plot all data on journeys back from L4 , if listed below ignore, wrong or didnt work for reasons above
% if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02');
if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
    
elseif strcmpi(fList(j,1:10),'2016-07-07')
        hold on
    load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'Licor_fco2_combined','undind1','undind2','undind3','undind4','Salinity_interp','Licor_Datetime','dist');
    scatter(Licor_fco2_combined(undind3:undind4)-Licor_fco2_combined(undind3),Salinity_interp(undind3:undind4)-Salinity_interp(undind3),8,dist(undind3:undind4),'MarkerEdgeColor',colour_greyshade)
    xlabel(['\Delta fCO_2(',num2str(micro_symbol),'atm)'],'fontsize',22)
    set(gca,'fontsize',22)
    set(gca,'fontsize',22)
    ylabel('\Delta Salinity','fontsize',22)
    hj=colorbar;
    ylabel(hj,'Distance from L4 (km)','FontSize',16);
    set(gca,'fontsize',16)
    set(gca,'fontsize',16)
    hold off
    
else
    hold on
    load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'Licor_fco2_combined','undind1','undind2','undind3','undind4','Salinity_interp','Licor_Datetime','dist');
    scatter(Licor_fco2_combined(undind3:undind4)-Licor_fco2_combined(undind3),Salinity_interp(undind3:undind4)-Salinity_interp(undind3),8,dist(undind3:undind4))
    xlabel(['\Delta fCO_2(',num2str(micro_symbol),'atm)'],'fontsize',22)
    set(gca,'fontsize',22)
    set(gca,'fontsize',22)
    ylabel('\Delta Salinity','fontsize',22)
    hj=colorbar;
    ylabel(hj,'Distance from L4 (km)','FontSize',16);
    set(gca,'fontsize',16)
    set(gca,'fontsize',16)
    hold off
end
     
    %plot all data on journeys out to L4 , if listed below ignore, wrong or didnt work for reasons above
% if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-04')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')|strcmpi(fList(j,1:10),'2016-07-07');
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')

else
    hold on
    load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'Licor_fco2_combined','undind1','undind2','undind3','undind4','Salinity_interp','Licor_Datetime','dist');
    scatter(Licor_fco2_combined(undind1:undind2)-Licor_fco2_combined(undind2),Salinity_interp(undind1:undind2)-Salinity_interp(undind2),8,dist(undind1:undind2))
    xlabel(['\Delta fCO_2(',num2str(micro_symbol),'atm)'],'fontsize',22)
    set(gca,'fontsize',22)
    set(gca,'fontsize',22)
    ylabel('\Delta Salinity','fontsize',22)
    hj=colorbar;
    ylabel(hj,'Distance from L4 (km)','FontSize',16);
    set(gca,'fontsize',16)
    set(gca,'fontsize',16)
    hold off
        end  
  hold on
%        plot(Tamar_fco2_calc_ta_dic([6:7,9])-Tamar_fco2_calc_ta_dic(8),Tamar_Salinity([6:7,9])-Tamar_Salinity(8),'d','MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor','r')
% 
%         
%     labs=text(20+(Tamar_fco2_calc_ta_dic([6:7,9])-Tamar_fco2_calc_ta_dic(8)),(Tamar_Salinity([6:7,9])-Tamar_Salinity(8)), Tamarlabels([6:7,9]),'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
%     
end




% (8,3,[1 4 7 10 13 16])
% (8,3,[2 5 8 11 14 17])
% (8,3,[3 6 9 12 15 18])
% (8,3,[19:24])
%plot representative of low water +0:00hrs
figure(801)
set(gcf, 'Color','w','Position', get(0,'Screensize'));

%%
load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\2016-07-07.mat'],'Longitude_interp','Latitude_interp','undind1','undind2','undind3','undind4','Licor_fco2_combined','SeaTemp_interp','Salinity_interp','Licor_Datetime','dist');
    subsample=30;
    ax1=subplot(8,3,[1 4 7 10 13 16])%plot patch first underneath.
            text(0.02,1.07,'(a)','color','k','Fontsize',20,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-4.25:0.05:-4.1],'ytick',[50.25:0.025:50.35],'fontsize',13)

m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
% m_ruler([0.75 0.90],0.1,[0 1000 2000])
 m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,Salinity_interp(undind3:subsample:undind4),'filled')
    y=ylabel('Latitude','Fontsize',13);
    set(y,'Units','Normalized','Position',[-0.2,0.5,0]);
    x=xlabel('Longitude','Fontsize',13);
    set(x,'Units','Normalized','Position',[0.5,-0.05,0]); set(gca,'fontsize',13);
    set(gca,'fontsize',13);
    [C,D]=m_ll2xy(-4.217,50.25);    %add L4
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.217,50.25);
    text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.188801,50.318280);       %add penlee
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.188801,50.318280);
    text(CC,D,'PPAO','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.158932,50.334499);    %add breakwater
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.158932,50.334499);
    text(CC,D,'Breakwater','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
%     ring1=m_range_ring(-4.217,50.25,0.5); set(ring1,'color','r')ring2=m_range_ring(-4.217,50.25,1); set(ring2,'color','b');
    hj=colorbar;
    title(hj,['Salinity (PSU)'],'FontSize',13);   
    colormap(ax1,flipud(cm_data))

    
    ax2=subplot(8,3,[2 5 8 11 14 17])%plot patch first underneath.
m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-4.25:0.05:-4.1],'ytick',[50.25:0.025:50.35],'fontsize',13)
m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
% m_ruler([0.75 0.90],0.1,[0 1000 2000])
text(0.02,1.07,'(b)','color','k','Fontsize',20,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
markersty='p' ;
    m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,SeaTemp_interp(undind3:subsample:undind4),'filled'); hold on;
        y=ylabel('Latitude','Fontsize',13);
    set(y,'Units','Normalized','Position',[-0.2,0.5,0]);
    x=xlabel('Longitude','Fontsize',13);
    set(x,'Units','Normalized','Position',[0.5,-0.05,0]);
    set(gca,'fontsize',13)
    set(gca,'fontsize',13)
    [C,D]=m_ll2xy(-4.217,50.25);     %add L4
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.217,50.25);
    text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.188801,50.318280);    %add penlee
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.188801,50.318280);
    text(CC,D,'PPAO','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.158932,50.334499);    %add breakwater
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.158932,50.334499);
    text(CC,D,'Breakwater','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
%      ring1=m_range_ring(-4.217,50.25,0.5);
%     set(ring1,'color','r');ring2=m_range_ring(-4.217,50.25,1);set(ring2,'color','b')
     hj=colorbar;
    title(hj,['Temperature (',num2str(degree_symbol),'C)'],'FontSize',13);
    colormap(ax2,cm_data)

    
    ax3=subplot(8,3,[3 6 9 12 15 18]);%plot patch first underneath.
m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-4.25:0.05:-4.1],'ytick',[50.25:0.025:50.35],'fontsize',13)
m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
% m_ruler([0.75 0.90],0.1,[0 1000 2000])
m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,Licor_fco2_combined(undind3:subsample:undind4),'filled')
text(0.02,1.07,'(c)','color','k','Fontsize',20,'Fontweight','bold','BackgroundColor','w','units','normalized');     
 y=ylabel('Latitude','Fontsize',13);
    set(y,'Units','Normalized','Position',[-0.2,0.5,0]);
    x=xlabel('Longitude','Fontsize',13);
    set(x,'Units','Normalized','Position',[0.5,-0.05,0]);
ax3.YAxis.TickLabelFormat = '%.3f';
 set(gca,'fontsize',13)
    set(gca,'fontsize',13)
    [C,D]=m_ll2xy(-4.217,50.25);    %add L4
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.217,50.25);
    text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.188801,50.318280);    %add penlee
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.188801,50.318280);
    text(CC,D,'PPAO','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.158932,50.334499);    %add breakwater
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.158932,50.334499);
    text(CC,D,'Breakwater','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
%     ring1=m_range_ring(-4.217,50.25,0.5); set(ring1,'color','r');ring2=m_range_ring(-4.217,50.25,1);set(ring2,'color','b')
    hj=colorbar;
    title(hj,['fCO_2_(_s_w_)(',num2str(micro_symbol),'atm)'],'FontSize',13);
        colormap(ax3,cm_data)

        ax4=subplot(8,3,[19:24])%plot timeseies beneath
            colormap(ax4,cm_data)
plot(Licor_Datetime(undind3:subsample:undind4),Salinity_interp(undind3:subsample:undind4),'LineWidth',4,'color','black')
        hold on
        xlabel(['Time'],'FontSize',13);
        ylabel('Salinity (PSU)','Fontsize',13);
        set(gca, 'XTick', 736518+(15/24):1/(24*12):736518+(15/24)+10/(24*12),'FontSize',13)%15:00 to 15:50, every 5 minutes
        datetick('x','HH:MM', 'keepticks')
        addaxissubplot(Licor_Datetime(undind3:subsample:undind4),SeaTemp_interp(undind3:subsample:undind4),'LineWidth',4)
        pk=addaxislabel(2,['Temperature (',num2str(degree_symbol),'C)'],'FontSize',13);
        set(pk,'FontSize',13)
        hold on
        addaxissubplot(Licor_Datetime(undind3:subsample:undind4),Licor_fco2_combined(undind3:subsample:undind4),'LineWidth',4)
        pk2=addaxislabel(3,['fCO_2_(_s_w_)(',num2str(micro_symbol),'atm)'],'FontSize',13);
        set(pk2,'FontSize',13)
        AX=findall(0,'type','axes'); 
        set(AX(2),'fontsize',13)
        set(AX(3),'fontsize',13)
        set(ax4,'fontsize',13)
        set(ax4,'fontsize',13)
        text(0.02,1.16,'(d)','color','k','Fontsize',20,'Fontweight','bold','BackgroundColor','w','units','normalized');
        set(gca, 'ColorOrder', cm_data, 'NextPlot', 'replacechildren');

 
%%

%plot representative of low water +3:00hrs
figure(802)
set(gcf, 'Color','w','Position', get(0,'Screensize'));

%%
load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\2016-06-15.mat'],'Longitude_interp','Latitude_interp','undind1','undind2','undind3','undind4','Licor_fco2_combined','SeaTemp_interp','Salinity_interp','Licor_Datetime','dist');
    subsample=30;
    ax1=subplot(8,3,[1 4 7 10 13 16])%plot patch first underneath.
m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-4.25:0.05:-4.1],'ytick',[50.25:0.025:50.35],'fontsize',13)
m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
% m_ruler([0.75 0.90],0.1,[0 1000 2000])
 m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,Salinity_interp(undind3:subsample:undind4),'filled')
text(0.02,1.07,'(a)','color','k','Fontsize',20,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
 y=ylabel('Latitude');
    set(y,'Units','Normalized','Position',[-0.2,0.5,0]);
    x=xlabel('Longitude');
    set(x,'Units','Normalized','Position',[0.5,-0.05,0]);
    set(gca,'fontsize',13);
    set(gca,'fontsize',13);
    [C,D]=m_ll2xy(-4.217,50.25);    %add L4
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.217,50.25);
    text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.188801,50.318280);       %add penlee
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.188801,50.318280);
    text(CC,D,'PPAO','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.158932,50.334499);    %add breakwater
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.158932,50.334499);
    text(CC,D,'Breakwater','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
%     ring1=m_range_ring(-4.217,50.25,0.5); set(ring1,'color','r')ring2=m_range_ring(-4.217,50.25,1); set(ring2,'color','b');
    hj=colorbar;
    title(hj,['Salinity (PSU)'],'FontSize',13);
    colormap(ax1,flipud(cm_data))

    
    ax2=subplot(8,3,[2 5 8 11 14 17])%plot patch first underneath.
m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-4.25:0.05:-4.1],'ytick',[50.25:0.025:50.35],'fontsize',13)
m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
% m_ruler([0.75 0.90],0.1,[0 1000 2000])
markersty='p' ;            
text(0.02,1.07,'(b)','color','k','Fontsize',20,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
    m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,SeaTemp_interp(undind3:subsample:undind4),'filled'); hold on;
            y=ylabel('Latitude');
    set(y,'Units','Normalized','Position',[-0.2,0.5,0]);
    x=xlabel('Longitude');
    set(x,'Units','Normalized','Position',[0.5,-0.05,0]);
    set(gca,'fontsize',13)
    set(gca,'fontsize',13)
    [C,D]=m_ll2xy(-4.217,50.25);     %add L4
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.217,50.25);
    text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.188801,50.318280);    %add penlee
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.188801,50.318280);
    text(CC,D,'PPAO','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.158932,50.334499);    %add breakwater
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.158932,50.334499);
    text(CC,D,'Breakwater','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
%      ring1=m_range_ring(-4.217,50.25,0.5);
%     set(ring1,'color','r');ring2=m_range_ring(-4.217,50.25,1);set(ring2,'color','b')
     hj=colorbar;
    title(hj,['Temperature (',num2str(degree_symbol),'C)'],'FontSize',13);
    colormap(ax2,cm_data)

    
    ax3=subplot(8,3,[3 6 9 12 15 18]);%plot patch first underneath.
m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-4.25:0.05:-4.1],'ytick',[50.25:0.025:50.35],'fontsize',13)
m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
% m_ruler([0.75 0.90],0.1,[0 1000 2000])
 m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,Licor_fco2_combined(undind3:subsample:undind4),'filled')
  text(0.02,1.07,'(c)','color','k','Fontsize',20,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
       y=ylabel('Latitude','FontSize',13);
    set(y,'Units','Normalized','Position',[-0.2,0.5,0]);
    x=xlabel('Longitude','FontSize',13);
    set(x,'Units','Normalized','Position',[0.5,-0.05,0]);
 set(gca,'fontsize',13)
    set(gca,'fontsize',13)
    [C,D]=m_ll2xy(-4.217,50.25);    %add L4
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.217,50.25);
    text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.188801,50.318280);    %add penlee
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.188801,50.318280);
    text(CC,D,'PPAO','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.158932,50.334499);    %add breakwater
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.158932,50.334499);
    text(CC,D,'Breakwater','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
%     ring1=m_range_ring(-4.217,50.25,0.5); set(ring1,'color','r');ring2=m_range_ring(-4.217,50.25,1);set(ring2,'color','b')
    hj=colorbar;
    title(hj,['fCO_2_(_s_w_)(',num2str(micro_symbol),'atm)'],'FontSize',13);
    colormap(ax3,cm_data)

    
        ax4=subplot(8,3,[19:24])%plot timeseies beneath
            colormap(ax4,cm_data)
plot(Licor_Datetime(undind3:subsample:undind4),Salinity_interp(undind3:subsample:undind4),'LineWidth',4,'color','black')
        hold on
        xlabel(['Time'],'FontSize',13);
        ylabel('Salinity (PSU)','Fontsize',13);
        set(gca, 'XTick', 736496+(12/24):1/(24*12):736496+(12/24)+8/(24*12),'FontSize',13)%12:00 to 12:40, every 5 minutes
        datetick('x','HH:MM', 'keepticks')
        addaxissubplot(Licor_Datetime(undind3:subsample:undind4),SeaTemp_interp(undind3:subsample:undind4),'LineWidth',4)
        pk=addaxislabel(2,['Temperature (',num2str(degree_symbol),'C)'],'FontSize',13);
        set(pk,'FontSize',13)
        hold on
        addaxissubplot(Licor_Datetime(undind3:subsample:undind4),Licor_fco2_combined(undind3:subsample:undind4),'LineWidth',4)
        pk2=addaxislabel(3,['fCO_2_(_s_w_)(',num2str(micro_symbol),'atm)'],'FontSize',13);
        set(pk2,'FontSize',13)
        AX=findall(0,'type','axes'); 
        set(AX(2),'fontsize',13)
        set(AX(3),'fontsize',13)
        set(ax4,'fontsize',13)
        set(ax4,'fontsize',13)
        text(0.02,1.16,'(d)','color','k','Fontsize',20,'Fontweight','bold','BackgroundColor','w','units','normalized');

    
%%
%plot representative of low water +6:00hrs
figure(803)
set(gcf, 'Color','w','Position', get(0,'Screensize'));

%%
load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\2016-06-30.mat'],'Longitude_interp','Latitude_interp','undind1','undind2','undind3','undind4','Licor_fco2_combined','SeaTemp_interp','Salinity_interp','Licor_Datetime','dist');
      subsample=30;
    ax1=subplot(8,3,[1 4 7 10 13 16])%plot patch first underneath.
m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-4.25:0.05:-4.1],'ytick',[50.25:0.025:50.35],'fontsize',13)
m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
% m_ruler([0.75 0.90],0.1,[0 1000 2000])
 m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,Salinity_interp(undind3:subsample:undind4),'filled')
text(0.02,1.07,'(a)','color','k','Fontsize',20,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
 y=ylabel('Latitude');
    set(y,'Units','Normalized','Position',[-0.2,0.5,0]);
    x=xlabel('Longitude');
    set(x,'Units','Normalized','Position',[0.5,-0.05,0]);
    set(gca,'fontsize',13);
    set(gca,'fontsize',13);
    [C,D]=m_ll2xy(-4.217,50.25);    %add L4
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.217,50.25);
    text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.188801,50.318280);       %add penlee
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.188801,50.318280);
    text(CC,D,'PPAO','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.158932,50.334499);    %add breakwater
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.158932,50.334499);
    text(CC,D,'Breakwater','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
%     ring1=m_range_ring(-4.217,50.25,0.5); set(ring1,'color','r')ring2=m_range_ring(-4.217,50.25,1); set(ring2,'color','b');
    hj=colorbar;
    title(hj,['Salinity (PSU)'],'FontSize',13);
    colormap(ax1,flipud(cm_data))

    
    ax2=subplot(8,3,[2 5 8 11 14 17])%plot patch first underneath.
m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-4.25:0.05:-4.1],'ytick',[50.25:0.025:50.35],'fontsize',13)
m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
% m_ruler([0.75 0.90],0.1,[0 1000 2000])
markersty='p' ;            
text(0.02,1.07,'(b)','color','k','Fontsize',20,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
    m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,SeaTemp_interp(undind3:subsample:undind4),'filled'); hold on;
            y=ylabel('Latitude');
    set(y,'Units','Normalized','Position',[-0.2,0.5,0]);
    x=xlabel('Longitude');
    set(x,'Units','Normalized','Position',[0.5,-0.05,0]);
    set(gca,'fontsize',13)
    set(gca,'fontsize',13)
    [C,D]=m_ll2xy(-4.217,50.25);     %add L4
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.217,50.25);
    text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.188801,50.318280);    %add penlee
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.188801,50.318280);
    text(CC,D,'PPAO','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.158932,50.334499);    %add breakwater
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.158932,50.334499);
    text(CC,D,'Breakwater','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
%      ring1=m_range_ring(-4.217,50.25,0.5);
%     set(ring1,'color','r');ring2=m_range_ring(-4.217,50.25,1);set(ring2,'color','b')
     hj=colorbar;
    title(hj,['Temperature (',num2str(degree_symbol),'C)'],'FontSize',13);
    colormap(ax2,cm_data)

    
    ax3=subplot(8,3,[3 6 9 12 15 18]);%plot patch first underneath.
m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-4.25:0.05:-4.1],'ytick',[50.25:0.025:50.35],'fontsize',13)
m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
% m_ruler([0.75 0.90],0.1,[0 1000 2000])
 m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,Licor_fco2_combined(undind3:subsample:undind4),'filled')
  text(0.02,1.07,'(c)','color','k','Fontsize',20,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
       y=ylabel('Latitude','FontSize',13);
    set(y,'Units','Normalized','Position',[-0.2,0.5,0]);
    x=xlabel('Longitude','FontSize',13);
    set(x,'Units','Normalized','Position',[0.5,-0.05,0]);
 set(gca,'fontsize',13)
    set(gca,'fontsize',13)
    [C,D]=m_ll2xy(-4.217,50.25);    %add L4
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.217,50.25);
    text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.188801,50.318280);    %add penlee
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.188801,50.318280);
    text(CC,D,'PPAO','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.158932,50.334499);    %add breakwater
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.158932,50.334499);
    text(CC,D,'Breakwater','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
%     ring1=m_range_ring(-4.217,50.25,0.5); set(ring1,'color','r');ring2=m_range_ring(-4.217,50.25,1);set(ring2,'color','b')
    hj=colorbar;
    title(hj,['fCO_2_(_s_w_)(',num2str(micro_symbol),'atm)'],'FontSize',13);
    colormap(ax3,cm_data)

    
        ax4=subplot(8,3,[19:24])%plot timeseies beneath
            colormap(ax4,cm_data)
plot(Licor_Datetime(undind3:subsample:undind4),Salinity_interp(undind3:subsample:undind4),'LineWidth',4,'color','black')
        hold on
        xlabel(['Time'],'FontSize',13);
        ylabel('Salinity (PSU)','Fontsize',13);
        set(gca, 'XTick', 736511+(14/24)+4/(24*12):1/(24*12):736511+(14/24)+10/(24*12),'FontSize',13)%14:20 to 14:50, every 5 minutes
        datetick('x','HH:MM', 'keepticks')
        addaxissubplot(Licor_Datetime(undind3:subsample:undind4),SeaTemp_interp(undind3:subsample:undind4),'LineWidth',4)
        pk=addaxislabel(2,['Temperature (',num2str(degree_symbol),'C)'],'FontSize',13);
        set(pk,'FontSize',13)
        hold on
        addaxissubplot(Licor_Datetime(undind3:subsample:undind4),Licor_fco2_combined(undind3:subsample:undind4),'LineWidth',4)
        pk2=addaxislabel(3,['fCO_2_(_s_w_)(',num2str(micro_symbol),'atm)'],'FontSize',13);
        set(pk2,'FontSize',13)
        AX=findall(0,'type','axes'); 
        set(AX(2),'fontsize',13)
        set(AX(3),'fontsize',13)
        set(ax4,'fontsize',13)
        set(ax4,'fontsize',13)
        text(0.02,1.16,'(d)','color','k','Fontsize',20,'Fontweight','bold','BackgroundColor','w','units','normalized');

%%
    
%plot representative of low water +9:00hrs
figure(804)
set(gcf, 'Color','w','Position', get(0,'Screensize'));

%%
load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\2016-06-10.mat'],'Longitude_interp','Latitude_interp','undind1','undind2','undind3','undind4','Licor_fco2_combined','SeaTemp_interp','Salinity_interp','Licor_Datetime','dist');
        subsample=30;
    ax1=subplot(8,3,[1 4 7 10 13 16])%plot patch first underneath.
m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-4.25:0.05:-4.1],'ytick',[50.25:0.025:50.35],'fontsize',13)
m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
% m_ruler([0.75 0.90],0.1,[0 1000 2000])
 m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,Salinity_interp(undind3:subsample:undind4),'filled')
text(0.02,1.07,'(a)','color','k','Fontsize',20,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
 y=ylabel('Latitude');
    set(y,'Units','Normalized','Position',[-0.2,0.5,0]);
    x=xlabel('Longitude');
    set(x,'Units','Normalized','Position',[0.5,-0.05,0]);
    set(gca,'fontsize',13);
    set(gca,'fontsize',13);
    [C,D]=m_ll2xy(-4.217,50.25);    %add L4
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.217,50.25);
    text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.188801,50.318280);       %add penlee
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.188801,50.318280);
    text(CC,D,'PPAO','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.158932,50.334499);    %add breakwater
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.158932,50.334499);
    text(CC,D,'Breakwater','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
%     ring1=m_range_ring(-4.217,50.25,0.5); set(ring1,'color','r')ring2=m_range_ring(-4.217,50.25,1); set(ring2,'color','b');
    hj=colorbar;
    title(hj,['Salinity (PSU)'],'FontSize',13);
    colormap(ax1,flipud(cm_data))

    
    ax2=subplot(8,3,[2 5 8 11 14 17])%plot patch first underneath.
m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-4.25:0.05:-4.1],'ytick',[50.25:0.025:50.35],'fontsize',13)
m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
% m_ruler([0.75 0.90],0.1,[0 1000 2000])
markersty='p' ;            
text(0.02,1.07,'(b)','color','k','Fontsize',20,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
    m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,SeaTemp_interp(undind3:subsample:undind4),'filled'); hold on;
            y=ylabel('Latitude');
    set(y,'Units','Normalized','Position',[-0.2,0.5,0]);
    x=xlabel('Longitude');
    set(x,'Units','Normalized','Position',[0.5,-0.05,0]);
    set(gca,'fontsize',13)
    set(gca,'fontsize',13)
    [C,D]=m_ll2xy(-4.217,50.25);     %add L4
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.217,50.25);
    text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.188801,50.318280);    %add penlee
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.188801,50.318280);
    text(CC,D,'PPAO','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.158932,50.334499);    %add breakwater
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.158932,50.334499);
    text(CC,D,'Breakwater','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
%      ring1=m_range_ring(-4.217,50.25,0.5);
%     set(ring1,'color','r');ring2=m_range_ring(-4.217,50.25,1);set(ring2,'color','b')
     hj=colorbar;
    title(hj,['Temperature (',num2str(degree_symbol),'C)'],'FontSize',13);
    colormap(ax2,cm_data)

    
    ax3=subplot(8,3,[3 6 9 12 15 18]);%plot patch first underneath.
m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-4.25:0.05:-4.1],'ytick',[50.25:0.025:50.35],'fontsize',13)
m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
% m_ruler([0.75 0.90],0.1,[0 1000 2000])
 m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,Licor_fco2_combined(undind3:subsample:undind4),'filled')
  text(0.02,1.07,'(c)','color','k','Fontsize',20,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
       y=ylabel('Latitude','FontSize',13);
    set(y,'Units','Normalized','Position',[-0.2,0.5,0]);
    x=xlabel('Longitude','FontSize',13);
    set(x,'Units','Normalized','Position',[0.5,-0.05,0]);
 set(gca,'fontsize',13)
    set(gca,'fontsize',13)
    [C,D]=m_ll2xy(-4.217,50.25);    %add L4
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.217,50.25);
    text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.188801,50.318280);    %add penlee
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.188801,50.318280);
    text(CC,D,'PPAO','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
    [C,D]=m_ll2xy(-4.158932,50.334499);    %add breakwater
    line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
    [CC,D]=m_ll2xy(-4.158932,50.334499);
    text(CC,D,'Breakwater','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',13);
%     ring1=m_range_ring(-4.217,50.25,0.5); set(ring1,'color','r');ring2=m_range_ring(-4.217,50.25,1);set(ring2,'color','b')
    hj=colorbar;
    title(hj,['fCO_2_(_s_w_)(',num2str(micro_symbol),'atm)'],'FontSize',13);
    colormap(ax3,cm_data)

    
        ax4=subplot(8,3,[19:24])%plot timeseies beneath
            colormap(ax4,cm_data)
plot(Licor_Datetime(undind3:subsample:undind4),Salinity_interp(undind3:subsample:undind4),'LineWidth',4,'color','black')
        hold on
        xlabel(['Time'],'FontSize',13);
        ylabel('Salinity (PSU)','Fontsize',13);
        set(gca, 'XTick', 736491+(12/24)+9/(24*12):1/(24*12):736491+(13/24)+3/(24*12),'FontSize',13)%12:45 to 13:15, every 5 minutes
        datetick('x','HH:MM', 'keepticks')
        addaxissubplot(Licor_Datetime(undind3:subsample:undind4),SeaTemp_interp(undind3:subsample:undind4),'LineWidth',4)
        pk=addaxislabel(2,['Temperature (',num2str(degree_symbol),'C)'],'FontSize',13);
        set(pk,'FontSize',13)
        hold on
        addaxissubplot(Licor_Datetime(undind3:subsample:undind4),Licor_fco2_combined(undind3:subsample:undind4),'LineWidth',4)
        pk2=addaxislabel(3,['fCO_2_(_s_w_)(',num2str(micro_symbol),'atm)'],'FontSize',13);
        set(pk2,'FontSize',13)
        AX=findall(0,'type','axes'); 
        set(AX(2),'fontsize',13)
        set(AX(3),'fontsize',13)
        set(ax4,'fontsize',13)
        set(ax4,'fontsize',13)
        text(0.02,1.16,'(d)','color','k','Fontsize',20,'Fontweight','bold','BackgroundColor','w','units','normalized');

%%


%plot of all lat/long coloured by date
acaba=figure(116)
subsample=30;
m_proj('mercator','lon',[-4.3 -4.05],'lat',[50.235 50.50]);
m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-4.3:0.05:-4.05],'ytick',[50.225:0.025:50.55],'fontsize',22)
m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
% m_gshhs_f('patch',[.5 .5 .5]);
% m_ruler([0.75 0.90],0.1,[0 1000 2000])
hold on
for j = 1:length(fList(:,1));
    %plot all data on journeys back from L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
    else
        hold on
        load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_ret','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','Salinity_interp','undind1','undind2','undind3','undind4','ind''Salinity_interp','Licor_Datetime','dist');
        m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,Licor_Datetime(undind3:subsample:undind4))
        hold on
%         m_text(Longitude_interp(undind4)*-1,Latitude_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
        hold off
    end
    %plot all data on journeys out to L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')
    else
        hold on
        load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_ret','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','undind1','undind2','undind3','undind4','ind''Salinity_interp','Licor_Datetime','dist');
        m_scatter(Longitude_interp((undind1:subsample:undind2))*-1,Latitude_interp((undind1:subsample:undind2)),25,Licor_Datetime(undind1:subsample:undind2))
        hold on
%         m_text(Longitude_interp(undind2)*-1,Latitude_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
        hold off
    end
end
set(gca,'fontsize',22)
y=ylabel('Latitude');
set(y,'Units','Normalized','Position',[-0.28,0.5,0]);
set(gca,'fontsize',22);
x=xlabel('Longitude');
set(x,'Units','Normalized','Position',[0.5,-0.055,0]);
set(gca,'fontsize',22);
set(gca,'fontsize',22); hold on;
[C,D]=m_ll2xy(-4.217,50.25);    %add L4
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.219,50.25);
text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);

[C,D]=m_ll2xy(-4.188801,50.318280);    %add penlee
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.194801,50.318280);
text(CC,D,'PPAO','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% 
% [C,D]=m_ll2xy(-4.140932,50.334499);    %add breakwater
% line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% [CC,D]=m_ll2xy(-4.153932,50.334499);
% text(CC,D,'Breakwater','HorizontalAlignment', 'left','VerticalAlignment', 'bottom','fontsize',18);
% 
% [C,D]=m_ll2xy(-4.195,50.375);    %add breakwater
% line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% [CC,D]=m_ll2xy(-4.195,50.375);
% text(CC,D,{['River' char(10) 'Tamar']},'HorizontalAlignment', 'right','VerticalAlignment', 'top','fontsize',18);

annotation(acaba,'textarrow',[0.535937499999998 0.538541666666666],...
    [0.50434188034188 0.571581196581196],'String',{'Breakwater'},'FontSize',18);

annotation(acaba,'textarrow',[0.432812499999997 0.48125],...
    [0.792803418803418 0.725427350427351],'String',{'River','Tamar'},...
    'FontSize',18);

annotation(acaba,'textarrow',[0.432812499999997 0.48125],...
    [0.792803418803418 0.725427350427351],'String',{'Cawsand','Bay'},...
    'FontSize',18);

[C,D]=m_ll2xy(-4.2230,50.3124);    %add Rame
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.2180,50.3114);
text(CC,D,{['Rame' char(10) 'Head']},'HorizontalAlignment', 'right','VerticalAlignment', 'top','fontsize',18);

% [C,D]=m_ll2xy(-4.222,50.531);    %add Gunnislake gauge
% line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% [CC,D]=m_ll2xy(-4.222,50.531);
% text(CC,D,{['Gunnislake' char(10) 'flowrate gauge']},'HorizontalAlignment', 'left','VerticalAlignment', 'top','fontsize',18);

%this is off map
[C,D]=m_ll2xy(-4.185,50.368);    %add Gunnislake gauge
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.185,50.368);
text(CC,D,{['Devonport' char(10) 'tidal gauge']},'HorizontalAlignment', 'left','VerticalAlignment', 'bottom','fontsize',18);

[C,D]=m_ll2xy(-Tamar_LongW,Tamar_LatN);    %add bott;e locations
m_scatter(-Tamar_LongW,Tamar_LatN,25,'rd','filled');
[CC,D]=m_ll2xy(-Tamar_LongW,Tamar_LatN);
for V=1:length(CC)-3
text(CC(V),D(V),Tamarlabels(V),'HorizontalAlignment', 'left','VerticalAlignment', 'top','fontsize',12);
end
u=colorbar; 
set(u, 'YTick', 736482:7:736598)
cbdate(u,'mmm/dd');
ylabel(u,'Date(mmm/dd)','FontSize',22);    


%plot of all lat/long coloured by date
acab=figure(1161)
subsample=30;
m_proj('mercator','lon',[-4.25 -4.10],'lat',[50.235 50.4]);
m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-4.3:0.05:-4.05],'ytick',[50.225:0.025:50.55],'fontsize',22)
m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
% m_gshhs_f('patch',[.5 .5 .5]);
% m_ruler([0.75 0.90],0.1,[0 1000 2000])
hold on
for j = 1:length(fList(:,1));
    %plot all data on journeys back from L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
    else
        hold on
        load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_ret','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','Salinity_interp','undind1','undind2','undind3','undind4','ind''Salinity_interp','Licor_Datetime','dist');
        m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,Licor_Datetime(undind3:subsample:undind4))
        hold on
%         m_text(Longitude_interp(undind4)*-1,Latitude_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
        hold off
    end
    %plot all data on journeys out to L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')
    else
        hold on
        load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_ret','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','undind1','undind2','undind3','undind4','ind''Salinity_interp','Licor_Datetime','dist');
        m_scatter(Longitude_interp((undind1:subsample:undind2))*-1,Latitude_interp((undind1:subsample:undind2)),25,Licor_Datetime(undind1:subsample:undind2))
        hold on
%         m_text(Longitude_interp(undind2)*-1,Latitude_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
        hold off
    end
end
set(gca,'fontsize',22)
y=ylabel('Latitude');
set(y,'Units','Normalized','Position',[-0.28,0.5,0]);
set(gca,'fontsize',22);
x=xlabel('Longitude');
set(x,'Units','Normalized','Position',[0.5,-0.055,0]);
set(gca,'fontsize',22);
set(gca,'fontsize',22); hold on;
[C,D]=m_ll2xy(-4.217,50.25);    %add L4
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.219,50.25);
text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);

[C,D]=m_ll2xy(-4.188801,50.318280);    %add penlee
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.194801,50.318280);
text(CC,D,'PPAO','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);

% [C,D]=m_ll2xy(-4.140932,50.331499);    %add breakwater
% line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% [CC,D]=m_ll2xy(-4.153932,50.3311);
% text(CC,D,'Breakwater','HorizontalAlignment', 'left','VerticalAlignment', 'bottom','fontsize',18);
annotation(acab,'textarrow',[0.535937499999998 0.538541666666666],...
    [0.50434188034188 0.571581196581196],'String',{'Breakwater'},'FontSize',18);

% [C,D]=m_ll2xy(-4.195,50.375);    %add tamar
% line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% [CC,D]=m_ll2xy(-4.1536,50.360802);
% text(CC,D,{['River' char(10) 'Tamar']},'HorizontalAlignment', 'right','VerticalAlignment', 'top','fontsize',18);
annotation(acab,'textarrow',[0.432812499999997 0.48125],...
    [0.792803418803418 0.725427350427351],'String',{'River','Tamar'},...
    'FontSize',18);

annotation(acab,'textarrow',[0.432812499999997 0.48125],...
    [0.792803418803418 0.725427350427351],'String',{'Cawsand','Bay'},...
    'FontSize',18);

[C,D]=m_ll2xy(-4.2230,50.3124);    %add Rame
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.2180,50.3114);
text(CC,D,{['Rame' char(10) 'Head']},'HorizontalAlignment', 'right','VerticalAlignment', 'top','fontsize',18);

%this is off map
[C,D]=m_ll2xy(-4.185,50.368);    %add Gunnislake gauge
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.185,50.368);
text(CC,D,{['Devonport' char(10) 'tidal gauge']},'HorizontalAlignment', 'left','VerticalAlignment', 'bottom','fontsize',18);
u=colorbar; 
set(u, 'YTick', 736482:7:736598)
cbdate(u,'mmm/dd');
ylabel(u,'Date','FontSize',22);    hold on
ylabel(u,'Date(mmm/dd)','FontSize',22);    hold on
[C,D]=m_ll2xy(-Tamar_LongW(6:9),Tamar_LatN(6:9));    %add bott;e locations
m_scatter(-Tamar_LongW(6:9),Tamar_LatN(6:9),25,'rd','filled');
[CC,D]=m_ll2xy(-Tamar_LongW,Tamar_LatN);
for V=6:9
text(CC(V),D(V),Tamarlabels(V),'HorizontalAlignment', 'left','VerticalAlignment', 'top','fontsize',12);
end


%plot of all lat/long coloured by date
acaba=figure(1162)
subsample=30;
m_proj('mercator','lon',[-5.8 -2],'lat',[49.8 51.9]);
m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-6:0.5:-1.5],'ytick',[50:0.25:51.75],'fontsize',22)
k=m_line(coast(:,1),coast(:,2));
set(k,'color','k')
L(1)=m_line([-4.3 -4.05],[50.235 50.235],'color','k','linewidth',3);
m_line([-4.3 -4.05],[50.50 50.50],'color','k','linewidth',3);
m_line([-4.3 -4.3],[50.235 50.50],'color','k','linewidth',3);
m_line([-4.05 -4.05],[50.235 50.50],'color','k','linewidth',3);
% m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
% m_gshhs_f('patch',[.5 .5 .5]);
% m_ruler([0.75 0.90],0.1,[0 1000 2000])
hold on
for j = 1:length(fList(:,1));
    %plot all data on journeys back from L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
    else
        hold on
        load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_ret','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','Salinity_interp','undind1','undind2','undind3','undind4','ind''Salinity_interp','Licor_Datetime','dist');
        m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),5,Licor_Datetime(undind3:subsample:undind4))
        hold on
%         m_text(Longitude_interp(undind4)*-1,Latitude_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
        hold off
    end
    %plot all data on journeys out to L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')
    else
        hold on
        load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_ret','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','undind1','undind2','undind3','undind4','ind''Salinity_interp','Licor_Datetime','dist');
        m_scatter(Longitude_interp((undind1:subsample:undind2))*-1,Latitude_interp((undind1:subsample:undind2)),5,Licor_Datetime(undind1:subsample:undind2))
        hold on
%         m_text(Longitude_interp(undind2)*-1,Latitude_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
        hold off
    end
end
set(gca,'fontsize',22)
y=ylabel('Latitude');
set(y,'Units','Normalized','Position',[-0.12,0.5,0]);
set(gca,'fontsize',22);
x=xlabel('Longitude');
set(x,'Units','Normalized','Position',[0.5,-0.07,0]);
set(gca,'fontsize',22);
set(gca,'fontsize',22); hold on;

[C,D]=m_ll2xy(-Tamar_LongW,Tamar_LatN);    %add bott;e locations
m_scatter(-Tamar_LongW,Tamar_LatN,8,'rd','filled');
[CC,D]=m_ll2xy(-Tamar_LongW,Tamar_LatN);
% for V=1:length(CC)-3
% text(CC(V),D(V),Tamarlabels(V),'HorizontalAlignment', 'left','VerticalAlignment', 'top','fontsize',12);
% end




%plot of all lat/long coloured by date
acaba=figure(1163)

%plot of all lat/long coloured by date
subplot(1,2,1)
subsample=30;
m_proj('mercator','lon',[-5.8 -2],'lat',[49.8 51.9]);
m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-6:0.5:-1.5],'ytick',[50:0.25:51.75],'fontsize',22)
k=m_line(coast(:,1),coast(:,2));
set(k,'color','k')
L(1)=m_line([-4.3 -4.05],[50.235 50.235],'color','k','linewidth',3);
m_line([-4.3 -4.05],[50.50 50.50],'color','k','linewidth',3);
m_line([-4.3 -4.3],[50.235 50.50],'color','k','linewidth',3);
m_line([-4.05 -4.05],[50.235 50.50],'color','k','linewidth',3);
% m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
% m_gshhs_f('patch',[.5 .5 .5]);
% m_ruler([0.75 0.90],0.1,[0 1000 2000])
hold on
for j = 1:length(fList(:,1));
    %plot all data on journeys back from L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
    else
        hold on
        load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_ret','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','Salinity_interp','undind1','undind2','undind3','undind4','ind''Salinity_interp','Licor_Datetime','dist');
        m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),5,Licor_Datetime(undind3:subsample:undind4))
        hold on
%         m_text(Longitude_interp(undind4)*-1,Latitude_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
        hold off
    end
    %plot all data on journeys out to L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')
    else
        hold on
        load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_ret','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','undind1','undind2','undind3','undind4','ind''Salinity_interp','Licor_Datetime','dist');
        m_scatter(Longitude_interp((undind1:subsample:undind2))*-1,Latitude_interp((undind1:subsample:undind2)),5,Licor_Datetime(undind1:subsample:undind2))
        hold on
%         m_text(Longitude_interp(undind2)*-1,Latitude_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
        hold off
    end
end
set(gca,'fontsize',22)
y=ylabel('Latitude');
set(y,'Units','Normalized','Position',[-0.155,0.5,0]);
set(gca,'fontsize',22);
x=xlabel('Longitude');
set(x,'Units','Normalized','Position',[0.5,-0.085,0]);
set(gca,'fontsize',22);
set(gca,'fontsize',22); hold on;

[C,D]=m_ll2xy(-Tamar_LongW,Tamar_LatN);    %add bott;e locations
m_scatter(-Tamar_LongW,Tamar_LatN,8,'rd','filled');
[CC,D]=m_ll2xy(-Tamar_LongW,Tamar_LatN);
% for V=1:length(CC)-3
% text(CC(V),D(V),Tamarlabels(V),'HorizontalAlignment', 'left','VerticalAlignment', 'top','fontsize',12);
% end
text(0.02,1.07,'(a)','color','k','Fontsize',20,'Fontweight','bold','BackgroundColor','w','units','normalized'); 

subplot(1,2,2)
subsample=30;
m_proj('mercator','lon',[-4.3 -4.05],'lat',[50.235 50.50]);
m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-4.3:0.05:-4.05],'ytick',[50.225:0.025:50.55],'fontsize',22)
m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
% m_gshhs_f('patch',[.5 .5 .5]);
% m_ruler([0.75 0.90],0.1,[0 1000 2000])
hold on
for j = 1:length(fList(:,1));
    %plot all data on journeys back from L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
    else
        hold on
        load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_ret','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','Salinity_interp','undind1','undind2','undind3','undind4','ind''Salinity_interp','Licor_Datetime','dist');
        m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,Licor_Datetime(undind3:subsample:undind4))
        hold on
%         m_text(Longitude_interp(undind4)*-1,Latitude_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
        hold off
    end
    %plot all data on journeys out to L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')
    else
        hold on
        load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_ret','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','undind1','undind2','undind3','undind4','ind''Salinity_interp','Licor_Datetime','dist');
        m_scatter(Longitude_interp((undind1:subsample:undind2))*-1,Latitude_interp((undind1:subsample:undind2)),25,Licor_Datetime(undind1:subsample:undind2))
        hold on
%         m_text(Longitude_interp(undind2)*-1,Latitude_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
        hold off
    end
end
set(gca,'fontsize',22)
y=ylabel('Latitude');
set(y,'Units','Normalized','Position',[-0.28,0.5,0]);
set(gca,'fontsize',22);
x=xlabel('Longitude');
set(x,'Units','Normalized','Position',[0.5,-0.055,0]);
set(gca,'fontsize',22);
set(gca,'fontsize',22); hold on;
[C,D]=m_ll2xy(-4.217,50.25);    %add L4
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.219,50.25);
text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);

[C,D]=m_ll2xy(-4.188801,50.318280);    %add penlee
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.194801,50.318280);
text(CC,D,'PPAO','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% 
% [C,D]=m_ll2xy(-4.140932,50.334499);    %add breakwater
% line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% [CC,D]=m_ll2xy(-4.153932,50.334499);
% text(CC,D,'Breakwater','HorizontalAlignment', 'left','VerticalAlignment', 'bottom','fontsize',18);
% 
% [C,D]=m_ll2xy(-4.195,50.375);    %add breakwater
% line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% [CC,D]=m_ll2xy(-4.195,50.375);
% text(CC,D,{['River' char(10) 'Tamar']},'HorizontalAlignment', 'right','VerticalAlignment', 'top','fontsize',18);

annotation(acaba,'textarrow',[0.535937499999998 0.538541666666666],...
    [0.50434188034188 0.571581196581196],'String',{'Breakwater'},'FontSize',18);

annotation(acaba,'textarrow',[0.432812499999997 0.48125],...
    [0.792803418803418 0.725427350427351],'String',{'River','Tamar'},...
    'FontSize',18);

annotation(acaba,'textarrow',[0.432812499999997 0.48125],...
    [0.792803418803418 0.725427350427351],'String',{'Cawsand','Bay'},...
    'FontSize',18);

[C,D]=m_ll2xy(-4.2230,50.3124);    %add Rame
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.2180,50.3114);
text(CC,D,{['Rame' char(10) 'Head']},'HorizontalAlignment', 'right','VerticalAlignment', 'top','fontsize',18);

% [C,D]=m_ll2xy(-4.222,50.531);    %add Gunnislake gauge
% line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% [CC,D]=m_ll2xy(-4.222,50.531);
% text(CC,D,{['Gunnislake' char(10) 'flowrate gauge']},'HorizontalAlignment', 'left','VerticalAlignment', 'top','fontsize',18);

%this is off map
[C,D]=m_ll2xy(-4.185,50.368);    %add Gunnislake gauge
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.185,50.368);
text(CC,D,{['Devonport' char(10) 'tidal gauge']},'HorizontalAlignment', 'left','VerticalAlignment', 'bottom','fontsize',18);

[C,D]=m_ll2xy(-Tamar_LongW,Tamar_LatN);    %add bott;e locations
m_scatter(-Tamar_LongW,Tamar_LatN,25,'rd','filled');
[CC,D]=m_ll2xy(-Tamar_LongW,Tamar_LatN);
for V=1:length(CC)-3
text(CC(V),D(V),Tamarlabels(V),'HorizontalAlignment', 'left','VerticalAlignment', 'top','fontsize',12);
end
u=colorbar; 
set(u, 'YTick', 736482:7:736598)
cbdate(u,'mmm/dd');
ylabel(u,'Date(mmm/dd)','FontSize',22);    
text(0.02,1.07,'(b)','color','k','Fontsize',20,'Fontweight','bold','BackgroundColor','w','units','normalized'); 



%plot all salinity, temp and co2 tracks by lat and long
%%  
%plot of all salinity by lat/long
figure(117)
subplot(1,3,1)
subsample=30;
m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-4.25:0.05:-4.1],'ytick',[50.25:0.025:50.35],'fontsize',12)
m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
% m_ruler([0.75 0.90],0.1,[0 1000 2000])
hold on
text(0.02,1.05,'(a)','color','k','Fontsize',13,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
for j = 1:length(fList(:,1));
    %plot all data on journeys back from L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
    else
        hold on
        load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_ret','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','Salinity_interp','undind1','undind2','undind3','undind4','ind','Salinity_interp','Licor_Datetime','dist');
        m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,Salinity_interp(undind3:subsample:undind4))
        hold on
%         m_text(Longitude_interp(undind4)*-1,Latitude_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
        hold off
    end
    %plot all data on journeys out to L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')
    else
        hold on
        load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_ret','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','undind1','undind2','undind3','undind4','ind','Salinity_interp','Licor_Datetime','dist');
        m_scatter(Longitude_interp((undind1:subsample:undind2))*-1,Latitude_interp((undind1:subsample:undind2)),25,Salinity_interp(undind1:subsample:undind2))
        hold on
%         m_text(Longitude_interp(undind2)*-1,Latitude_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
        hold off
    end
end
set(gca,'fontsize',12)
y=ylabel('Latitude');
    set(y,'Units','Normalized','Position',[-0.2,0.5,0]);
    x=xlabel('Longitude');
    set(x,'Units','Normalized','Position',[0.5,-0.05,0]);
set(gca,'fontsize',12);
set(gca,'fontsize',12); hold on;
[C,D]=m_ll2xy(-4.217,50.25);    %add L4
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.217,50.243);
text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',12);
[C,D]=m_ll2xy(-4.188801,50.318280);    %add penlee
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.188801,50.318280);
text(CC,D,'PPAO','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',12);
[C,D]=m_ll2xy(-4.158932,50.334499);    %add breakwater
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.158932,50.334499);
text(CC,D,'Breakwater','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',12);
% ring1=m_range_ring(-4.217,50.25,0.5);
% set(ring1,'color','r')
% ring2=m_range_ring(-4.217,50.25,1);
% set(ring2,'color','b')
% ring3=m_range_ring(-4.217,50.25,2);
% set(ring3,'color','g');    
hj=colorbar;
title(hj,['Salinity'],'FontSize',12);
     
  
 
%plot of all temp by lat/long
subplot(1,3,2)
subsample=30;
m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-4.25:0.05:-4.1],'ytick',[50.25:0.025:50.35],'fontsize',12)
m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
% m_ruler([0.75 0.90],0.1,[0 1000 2000])
text(0.02,1.05,'(b)','color','k','Fontsize',13,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
hold on
for j = 1:length(fList(:,1));
    %plot all data on journeys back from L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
    else
        hold on
        load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_ret','SeaTemp_interp','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','Salinity_interp','undind1','undind2','undind3','undind4','ind','Salinity_interp','Licor_Datetime','dist');
        m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,SeaTemp_interp(undind3:subsample:undind4))
        hold on
%         m_text(Longitude_interp(undind4)*-1,Latitude_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
        hold off
    end
    %plot all data on journeys out to L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')
    else
        hold on
        load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_ret','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','undind1','undind2','undind3','undind4','ind','Salinity_interp','Licor_Datetime','dist');
        m_scatter(Longitude_interp((undind1:subsample:undind2))*-1,Latitude_interp((undind1:subsample:undind2)),25,SeaTemp_interp(undind1:subsample:undind2))
        hold on
%         m_text(Longitude_interp(undind2)*-1,Latitude_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
        hold off
    end
end
set(gca,'fontsize',12)
y=ylabel('Latitude');
    set(y,'Units','Normalized','Position',[-0.2,0.5,0]);
    x=xlabel('Longitude');
    set(x,'Units','Normalized','Position',[0.5,-0.05,0]);
set(gca,'fontsize',12);
set(gca,'fontsize',12); hold on;
[C,D]=m_ll2xy(-4.217,50.25);    %add L4
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.217,50.243);
text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',12);
[C,D]=m_ll2xy(-4.188801,50.318280);    %add penlee
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.188801,50.318280);
text(CC,D,'PPAO','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',12);
[C,D]=m_ll2xy(-4.158932,50.334499);    %add breakwater
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.158932,50.334499);
text(CC,D,'Breakwater','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',12);
% ring1=m_range_ring(-4.217,50.25,0.5);
% set(ring1,'color','r')
% ring2=m_range_ring(-4.217,50.25,1);
% set(ring2,'color','b')
% ring3=m_range_ring(-4.217,50.25,2);
% set(ring3,'color','g');    
    hj=colorbar;
    title(hj,['Temperature (',num2str(degree_symbol),'C)'],'FontSize',12);
     
  
%plot of all co2 by lat/long
subplot(1,3,3)
subsample=30;
m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-4.25:0.05:-4.1],'ytick',[50.25:0.025:50.35],'fontsize',12)
m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
% m_ruler([0.75 0.90],0.1,[0 1000 2000])
hold on
text(0.02,1.05,'(c)','color','k','Fontsize',13,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
for j = 1:length(fList(:,1));
    %plot all data on journeys back from L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
    else
        hold on
        load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_ret','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','undind1','undind2','undind3','undind4','ind''Salinity_interp','Licor_Datetime','dist');
        m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,Licor_fco2_combined(undind3:subsample:undind4))
        hold on
%         m_text(Longitude_interp(undind4)*-1,Latitude_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
        hold off
    end
    %plot all data on journeys out to L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')
    else
        hold on
        load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_ret','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','undind1','undind2','undind3','undind4','ind''Salinity_interp','Licor_Datetime','dist');
        m_scatter(Longitude_interp((undind1:subsample:undind2))*-1,Latitude_interp((undind1:subsample:undind2)),25,Licor_fco2_combined(undind1:subsample:undind2))
        hold on
%         m_text(Longitude_interp(undind2)*-1,Latitude_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
        hold off
    end
end
set(gca,'fontsize',12)
y=ylabel('Latitude');
set(y,'Units','Normalized','Position',[-0.2,0.5,0]);
set(gca,'fontsize',12);
    x=xlabel('Longitude');
    set(x,'Units','Normalized','Position',[0.5,-0.05,0]);set(gca,'fontsize',12);
set(gca,'fontsize',12); hold on;
[C,D]=m_ll2xy(-4.217,50.25);    %add L4
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.217,50.243);
text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',12);
[C,D]=m_ll2xy(-4.188801,50.318280);    %add penlee
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.188801,50.318280);
text(CC,D,'PPAO','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',12);
[C,D]=m_ll2xy(-4.158932,50.334499);    %add breakwater
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.158932,50.334499);
text(CC,D,'Breakwater','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',12);
% ring1=m_range_ring(-4.217,50.25,0.5);
% set(ring1,'color','r')
% ring2=m_range_ring(-4.217,50.25,1);
% set(ring2,'color','b')
% ring3=m_range_ring(-4.217,50.25,2);
% set(ring3,'color','g');    
hj=colorbar;
title(hj,['fCO_2_(_s_w_)(',num2str(micro_symbol),'atm)'],'FontSize',12);
%%
     

%use this code to plot transects by point along the tide rather than by
%order of date

    Sorttimes=[];Inout=[];fListtransects=[];
      %find all transects and their times
    for j = 1:length(fList(:,1));
        %find all data on journeys back from L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
        else
            load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_ret','SeaTemp_interp','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','Salinity_interp','undind1','undind2','undind3','undind4','ind','Salinity_interp','Licor_Datetime','dist','dt_since_lowtide_und_retHHMM');
            Sorttimes = [Sorttimes ; dt_since_lowtide_und_retHHMM];
            Inout = [Inout ; 1];
            fListtransects = [fListtransects ; fList(j,1:10)];
        end
        %find all data on journeys out to L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')
        else
            load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_outbound','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','dt_since_lowtide_und_outboundHHMM','undind1','undind2','undind3','undind4','ind','Salinity_interp','Licor_Datetime','dist','dt_since_lowtide_und_retHHMM');
            Sorttimes = [Sorttimes ; dt_since_lowtide_und_outboundHHMM];
            Inout = [Inout ; 2];
            fListtransects = [fListtransects ; fList(j,1:10)];
        end
    end


    numericalSorttimes=datenum(Sorttimes,'HH:MM');
    [p subplotnumorder]=sort(numericalSorttimes');
    [p  spnum]=sort(subplotnumorder);

    %plot of all temp by lat/long
    figure(120)
    subsample=30;counterF=1;
    hold on
    for j = 1:length(fList(:,1));
        %plot all data on journeys back from L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
        else
            hold on
            load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_ret','SeaTemp_interp','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','Salinity_interp','undind1','undind2','undind3','undind4','ind','Salinity_interp','Licor_Datetime','dist','dt_since_lowtide_und_retHHMM');
            Subplotcolumns=4;
            subplot(4,Subplotcolumns,spnum(counterF))%plot patch first underneath.
            m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
            m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-4.25:0.05:-4.1],'ytick',[50.25:0.025:50.35],'fontsize',12)
            m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
            m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,SeaTemp_interp(undind3:subsample:undind4),'filled')
            hold on
            %         m_text(Longitude_interp(undind4)*-1,Latitude_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
            set(gca,'fontsize',12)
            % y=ylabel('Latitude');
            % set(y,'Units','Normalized','Position',[0,0.5,0]);
            % x=xlabel('Longitude');
            % set(x,'Units','Normalized','Position',[0.5,0.065,0]);
            set(gca,'fontsize',12);
            set(gca,'fontsize',12); hold on;
            hj=colorbar;
            ylabel(hj,['Temperature (',num2str(degree_symbol),'C)'],'FontSize',12);
            title({[fList(j,9:10) '/' fList(j,6:7) 'LW+' dt_since_lowtide_und_retHHMM(1:5) 'hrs']},'FontSize',16);hold on
            hold off
            counterF=counterF+1;
        end
        %plot all data on journeys out to L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')
        else
            hold on
            load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_outbound','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','dt_since_lowtide_und_outboundHHMM','undind1','undind2','undind3','undind4','ind','Salinity_interp','Licor_Datetime','dist','dt_since_lowtide_und_retHHMM');
            Subplotcolumns=4;
            subplot(4,Subplotcolumns,spnum(counterF))%plot patch first underneath.
            m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
            m_grid('linestyle','none','tickdir','out', 'tickstyle','dd','xtick',[-4.25:0.05:-4.1],'ytick',[50.25:0.025:50.35],'fontsize',12)
            m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
            m_scatter(Longitude_interp((undind1:subsample:undind2))*-1,Latitude_interp((undind1:subsample:undind2)),25,SeaTemp_interp(undind1:subsample:undind2),'filled')
            hold on
            %         m_text(Longitude_interp(undind2)*-1,Latitude_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
            set(gca,'fontsize',12)
            % y=ylabel('Latitude');
            % set(y,'Units','Normalized','Position',[0,0.5,0]);
            % x=xlabel('Longitude');
            % set(x,'Units','Normalized','Position',[0.5,0.065,0]);
            set(gca,'fontsize',12);
            set(gca,'fontsize',12); hold on;
            hj=colorbar;
            ylabel(hj,['Temperature (',num2str(degree_symbol),'C)'],'FontSize',12);
            title({[fList(j,9:10) '/' fList(j,6:7) 'LW+' dt_since_lowtide_und_outboundHHMM(1:5) 'hrs']},'FontSize',16);hold on
            hold off
            counterF=counterF+1;
        end
    end
  
    %plot of all co2 by lat/long
    figure(121)
    subsample=30;counterF=1;
    hold on
    for j = 1:length(fList(:,1));
        %plot all data on journeys back from L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
        else
            hold on
            load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_ret','SeaTemp_interp','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','Salinity_interp','undind1','undind2','undind3','undind4','ind','Salinity_interp','Licor_Datetime','dist','dt_since_lowtide_und_retHHMM');
            Subplotcolumns=4;
            subplot(4,Subplotcolumns,spnum(counterF))%plot patch first underneath.
            m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
            m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-4.25:0.05:-4.1],'ytick',[50.25:0.025:50.35],'fontsize',12)
            m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
            m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,Licor_fco2_combined(undind3:subsample:undind4),'filled')
            hold on
            %         m_text(Longitude_interp(undind4)*-1,Latitude_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
            set(gca,'fontsize',12)
            % y=ylabel('Latitude');
            % set(y,'Units','Normalized','Position',[0,0.5,0]);
            % x=xlabel('Longitude');
            % set(x,'Units','Normalized','Position',[0.5,0.065,0]);
            set(gca,'fontsize',12);
            set(gca,'fontsize',12); hold on;
            hj=colorbar;
            ylabel(hj,['fCO_2_(_s_w_)' '(',num2str(micro_symbol),'atm)'],'FontSize',12);
            title({[fList(j,9:10) '/' fList(j,6:7) 'LW+' dt_since_lowtide_und_retHHMM(1:5) 'hrs']},'FontSize',16);hold on
            hold off
            counterF=counterF+1;
        end
        %plot all data on journeys out to L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')
        else
            hold on
            load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_outbound','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','dt_since_lowtide_und_outboundHHMM','undind1','undind2','undind3','undind4','ind','Salinity_interp','Licor_Datetime','dist','dt_since_lowtide_und_retHHMM');
            Subplotcolumns=4;
            subplot(4,Subplotcolumns,spnum(counterF))%plot patch first underneath.
            m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
            m_grid('linestyle','none','tickdir','out', 'tickstyle','dd','xtick',[-4.25:0.05:-4.1],'ytick',[50.25:0.025:50.35],'fontsize',12)
            m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
            m_scatter(Longitude_interp((undind1:subsample:undind2))*-1,Latitude_interp((undind1:subsample:undind2)),25,Licor_fco2_combined(undind1:subsample:undind2),'filled')
            hold on
            %         m_text(Longitude_interp(undind2)*-1,Latitude_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
            set(gca,'fontsize',12)
            % y=ylabel('Latitude');
            % set(y,'Units','Normalized','Position',[0,0.5,0]);
            % x=xlabel('Longitude');
            % set(x,'Units','Normalized','Position',[0.5,0.065,0]);
            set(gca,'fontsize',12);
            set(gca,'fontsize',12); hold on;
            hj=colorbar;
            ylabel(hj,['fCO_2_(_s_w_)' '(',num2str(micro_symbol),'atm)'],'FontSize',12);
            title({[fList(j,9:10) '/' fList(j,6:7) 'LW+' dt_since_lowtide_und_outboundHHMM(1:5) 'hrs']},'FontSize',16);hold on
            hold off
            counterF=counterF+1;
        end
    end

    %plot of all temp by lat/long
    figure(122)
    subsample=30;counterF=1;
    hold on
    for j = 1:length(fList(:,1));
        %plot all data on journeys back from L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
 else
            hold on
            load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_ret','SeaTemp_interp','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','Salinity_interp','undind1','undind2','undind3','undind4','ind','Salinity_interp','Licor_Datetime','dist','dt_since_lowtide_und_retHHMM');
            Subplotcolumns=4;
            subplot(4,Subplotcolumns,spnum(counterF))%plot patch first underneath.
            m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
            m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-4.25:0.05:-4.1],'ytick',[50.25:0.025:50.35],'fontsize',12)
            m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
            m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,Salinity_interp(undind3:subsample:undind4),'filled')
            hold on
            %         m_text(Longitude_interp(undind4)*-1,Latitude_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
            set(gca,'fontsize',12)
            % y=ylabel('Latitude');
            % set(y,'Units','Normalized','Position',[0,0.5,0]);
            % x=xlabel('Longitude');
            % set(x,'Units','Normalized','Position',[0.5,0.065,0]);
            set(gca,'fontsize',12);
            set(gca,'fontsize',12); hold on;
            hj=colorbar;
            ylabel(hj,'Salinity (PSU)','FontSize',12);
            title({[fList(j,9:10) '/' fList(j,6:7) 'LW+' dt_since_lowtide_und_retHHMM(1:5) 'hrs']},'FontSize',16);hold on
            hold off
            counterF=counterF+1;
        end
        %plot all data on journeys out to L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')
        else
            hold on
            load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_outbound','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','dt_since_lowtide_und_outboundHHMM','undind1','undind2','undind3','undind4','ind','Salinity_interp','Licor_Datetime','dist','dt_since_lowtide_und_retHHMM');
            Subplotcolumns=4;
            subplot(4,Subplotcolumns,spnum(counterF))%plot patch first underneath.
            m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
            m_grid('linestyle','none','tickdir','out', 'tickstyle','dd','xtick',[-4.25:0.05:-4.1],'ytick',[50.25:0.025:50.35],'fontsize',12)
            m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
            m_scatter(Longitude_interp((undind1:subsample:undind2))*-1,Latitude_interp((undind1:subsample:undind2)),25,Salinity_interp(undind1:subsample:undind2),'filled')
            hold on
            %         m_text(Longitude_interp(undind2)*-1,Latitude_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
            set(gca,'fontsize',12)
            % y=ylabel('Latitude');
            % set(y,'Units','Normalized','Position',[0,0.5,0]);
            % x=xlabel('Longitude');
            % set(x,'Units','Normalized','Position',[0.5,0.065,0]);
            set(gca,'fontsize',12);
            set(gca,'fontsize',12); hold on;
            hj=colorbar;
            ylabel(hj,'Salinity (PSU)','FontSize',12);
            title({[fList(j,9:10) '/' fList(j,6:7) 'LW+' dt_since_lowtide_und_outboundHHMM(1:5) 'hrs']},'FontSize',16);hold on
            hold off
            counterF=counterF+1;
        end
    end
 
%load in data for mooring figures
load ('C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\Seasonalstudyunderway.mat');
load ('C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\underwaypCO2Seasonalstudy.mat');
load ('C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\Seasonal_study_L4_mooring.mat');

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
deg = km2deg(1);
[lon,lat] = scircle1(-4.217,50.25,deg);
[in ~]= inpolygon(-1*Longitude,Latitude,lon,lat);
L4ind=find(in==1);

deg = km2deg(1);
[in2 ~]= inpolygon(underwayLong,underwayLat,lon,lat);
L4indco2=find(in2==1);

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


tamar_dtt=datenum(tamar_dt);
tamar_flowrate_mms3=str2double(tamar_flowrate_ms3);

[p o]=find(tamar_dtt>datenum('2016-04-01 00:00:00','yyyy-mm-dd HH:MM:SS') & tamar_dtt<datenum('2016-09-30 00:00:00','yyyy-mm-dd HH:MM:SS'));


% plot(-1*Longitude,Latitude,'r*')
% hold on
% plot(lon,lat)

[rows cols]=size(fList);
klp=datenum(2016,06,9,0,0,0)
jkl=datenum(2016,10,01,0,0,1)
[a ~]=find(Penlee_dt>klp & Penlee_dt<jkl);
L4_dt_PAR=L4_dt_PAR';
x1=ones(length(Latitude),1)*50.251;
x2=ones(length(Latitude),1)*-4.221;
for t=1:length(Longitude);
    dist2l4(t) = pos2dist(x1(t),x2(t),Latitude(t),-1*Longitude(t),2);
end
    
[bu ~]=find(dist2l4'<1);
[bx ~]=find(Und_DT>klp & Und_DT<jkl);
[by bt]=ismember(bu,bx)
bt(bt==0)=[]
btt=bu(by)


[b ~]=find(by==1)
[c ~]=find(L4_dt>klp & L4_dt<jkl);
[d ~]=find(L4_dt_PAR>klp & L4_dt_PAR<jkl);
[e ~]=find(tamar_dtt>klp & tamar_dtt<jkl);


%Wind and Par
figure(201)  
subplot(3,1,1)
set(gcf, 'Color','w','Position', get(0,'Screensize'));
for j=1:rows;
    path = 'C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\';
    fListfilename = [path fList(j,1:10) '.mat'];
    if datenum(fList(j,1:10),'yyyy-mm-dd')>klp
            %unique string sequence all days when want to compare NSOP data
    if strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-06-30')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-08-04')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-15')|strcmpi(fList(j,1:10),'2016-09-21')
     
    load(fListfilename,'inds','indt','Licor_Datetime')
    Licor_Datetime(inds(1))
    r=indt(end);
    x=[r inds(1) inds(1) r];
    y= [ max([L4Windspeed;Windspeedabs]) max([L4Windspeed;Windspeedabs]) 0 0];
    patch(Licor_Datetime(x),y,colour_greyshade); 
    clearvars indslag indtlag Licor_Datetime
    else
        end
    else
    end

end
hold on;
despiked_Windspeedabs=func_despike_phasespace3d(Windspeedabs);
despiked_Windspeedabs=func_despike_phasespace3d(L4Windspeed);
plot(Penlee_dt(a),Penlee_Wspd_m_s(a),'LineWidth',2);
scatter(Und_DT(btt),Windspeedabs(btt),4,'o','filled');
% plot(L4_dt(c),L4Windspeed(c),'LineWidth',2);
dynamicDateTicks([], [], 'mm/dd');
xlabel('Time (month/day)','FontSize',13);
ylabel({[' Wind speed {(ms^{-1})}']},'FontSize',13);
ylim([0 18.2])
% klp=datenum(2016,06,10,0,0,1)
% jkl=datenum(2016,09,30,0,0,1)
xlim([klp jkl])
set(gca,'FontSize',13);
set(gca,'FontSize',13);
text(0.02,0.95,'(a)','color','k','Fontsize',13,'Fontweight','bold','BackgroundColor','w','units','normalized'); 

subplot(3,1,2)
for i=1:rows;
    path = 'C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\';
    fListfilename = [path fList(i,1:10) '.mat'];
    if datenum(fList(i,1:10),'yyyy-mm-dd')>klp
                    %unique string sequence all days when want to compare NSOP data
    if strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-06-30')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-08-04')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-15')|strcmpi(fList(j,1:10),'2016-09-21')
     
    load(fListfilename,'inds','indt','Licor_Datetime')
    Licor_Datetime(inds(1))
    r=indt(end);
    x=[r inds(1) inds(1) r];
    y= [ max([Winddirabs;L4Winddir]) max([Winddirabs;L4Winddir]) 0 0];
    patch(Licor_Datetime(x),y,colour_greyshade); 
    clearvars indslag indtlag Licor_Datetime
    else
    end
    else 
    end
end
hold on;
plot(Penlee_dt(a),Penlee_Wdir_deg(a),'LineWidth',2);
scatter(Und_DT(btt),Winddirabs(btt),4,'o','filled');
plot(L4_dt(c),L4Winddir(c),'LineWidth',2);
dynamicDateTicks([], [], 'mm/dd');
xlabel('Time (month/day)','FontSize',13);
ylabel({[' Wind direction (',num2str(degree_symbol) ')']},'FontSize',13);
ylim([0 max(L4Winddir)])
set(gca,'FontSize',13);
set(gca,'FontSize',13);
text(0.02,0.95,'(b)','color','k','Fontsize',13,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
xlim([klp jkl])

subplot(3,1,3)
for i=1:rows;
    path = 'C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\';
    fListfilename = [path fList(i,1:10) '.mat'];
    if datenum(fList(i,1:10),'yyyy-mm-dd')>klp
                    %unique string sequence all days when want to compare NSOP data
    if strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-06-30')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-08-04')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-15')|strcmpi(fList(j,1:10),'2016-09-21')
     
    load(fListfilename,'inds','indt','Licor_Datetime')
    Licor_Datetime(inds(1))
    r=indt(end);
    x=[r inds(1) inds(1) r];
    y= [ 2640 2640 0 0];
    patch(Licor_Datetime(x),y,colour_greyshade); 
    clearvars indslag indtlag Licor_Datetime
    hold on
    else
    end
    else 
    end
end
% plot(L4_dt_PAR(d),L4PAR(d),'LineWidth',2);
hold on
plot(L4_dt_PAR(d),L4PAR(d),'LineWidth',2);
scatter(Und_DT(btt),Parport(btt),4,'o','filled');
dynamicDateTicks([], [], 'mm/dd');
xlabel('Time (month/day)','FontSize',13);
ylabel({['PAR (Wm^{-2})']},'FontSize',13);
% ylim([0 max(L4PAR)])
set(gca,'FontSize',13);
set(gca,'FontSize',13);
text(0.02,0.95,'(c)','color','k','Fontsize',13,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
xlim([klp jkl])


%%%%%%%%%%%%%%%%%Temperature and salinity and river flowrate
figure(202)    
subplot(3,1,1)
for j=1:rows;
    path = 'C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\';
    fListfilename = [path fList(j,1:10) '.mat'];
                        %unique string sequence all days when want to compare NSOP data
    if strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-06-30')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-08-04')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-15')|strcmpi(fList(j,1:10),'2016-09-21')
     
    load(fListfilename,'inds','indt','Licor_Datetime')
    Licor_Datetime(inds(1))
    r=indt(end);
    x=[r inds(1) inds(1) r];
    y= [ max(L4SST) max(L4SST) 0 0];
    patch(Licor_Datetime(x),y,colour_greyshade); 
    clearvars indslag indtlag Licor_Datetime
    else end
end
hold on;
plot(L4_dt(c),L4SST(c),'LineWidth',2);
scatter(Und_DT(btt),SeaTemp(btt),4,'o','filled');
ylabel({['Temperature (',num2str(degree_symbol),'C)']},'FontSize',13);
ylim([12 max(L4SST)])
dynamicDateTicks([], [], 'mm/dd');
xlim([klp jkl])
% setDateAxes(gca, 'XLim', [datenum('June 10, 2016') datenum('October 1, 2016')])
set(gca,'FontSize',13);
set(gca,'FontSize',13);
text(0.03,0.90,'(a)','color','k','Fontsize',13,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
xlabel('Time (month/day)','FontSize',13);


subplot(3,1,2)
for j=1:rows;
    path = 'C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\';
    fListfilename = [path fList(j,1:10) '.mat'];
                        %unique string sequence all days when want to compare NSOP data
    if strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-06-30')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-08-04')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-15')|strcmpi(fList(j,1:10),'2016-09-21')
     
    load(fListfilename,'inds','indt','Licor_Datetime')
    Licor_Datetime(inds(1))
    r=indt(end);
    x=[r inds(1) inds(1) r];
    y= [ max(L4Salinity) max(L4Salinity) 0 0];
    patch(Licor_Datetime(x),y,colour_greyshade); 
    clearvars indslag indtlag Licor_Datetime
    else 
    end
end
hold on;
plot(L4_dt(c),L4Salinity(c),'LineWidth',2);
scatter(Und_DT(btt),Salinity(btt),4,'o','filled');
dynamicDateTicks([], [], 'mm/dd');
xlabel('Time (month/day)','FontSize',13);
ylabel({['Salinity']},'FontSize',13);
ylim([34.5 max(L4Salinity)])
xlim([klp jkl])
set(gca,'FontSize',13);
set(gca,'FontSize',13);
text(0.03,0.90,'(b)','color','k','Fontsize',13,'Fontweight','bold','BackgroundColor','w','units','normalized'); 

[bv , bt]=find(underway_dist'<1);

subplot(3,1,3)
for j=1:rows;
    path = 'C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\';
    fListfilename = [path fList(j,1:10) '.mat'];
                        %unique string sequence all days when want to compare NSOP data
    if strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-06-30')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-08-04')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-15')|strcmpi(fList(j,1:10),'2016-09-21')
     
    load(fListfilename,'inds','indt','Licor_Datetime')
    Licor_Datetime(inds(1))
    r=indt(end);
    x=[r inds(1) inds(1) r];
    y= [ 435 435 0 0];
    patch(Licor_Datetime(x),y,colour_greyshade); 
    clearvars indslag indtlag Licor_Datetime
    else 
    end
end
hold on;
plot(0,0)
scatter(underway_DT(bt),underwayfCO2_sw(bt),20,'o','filled')
dynamicDateTicks([], [], 'mm/dd');
xlabel('Time (month/day)','FontSize',13);
ylabel({['fCO_2_(_s_w_)(',num2str(micro_symbol),'atm)']},'FontSize',13);
xlim([klp jkl])
ylim([330 435])
set(gca,'FontSize',13);
set(gca,'FontSize',13);
text(0.03,0.90,'(c)','color','k','Fontsize',13,'Fontweight','bold','BackgroundColor','w','units','normalized'); 



figure(203)
for j=1:rows;
    path = 'C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\';
    fListfilename = [path fList(j,1:10) '.mat'];
                        %unique string sequence all days when want to compare NSOP data
    if strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-06-30')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-08-04')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-15')|strcmpi(fList(j,1:10),'2016-09-21')
     
    load(fListfilename,'inds','indt','Licor_Datetime')
    Licor_Datetime(inds(1))
    r=indt(end);
    x=[r inds(1) inds(1) r];
    y= [ (max(tamar_flowrate_mms3(e))) (max(tamar_flowrate_mms3(e))) 0 0];
    patch(Licor_Datetime(x),y,colour_greyshade); 
    clearvars indslag indtlag Licor_Datetime
    else 
    end
end
hold on
plot(tamar_dtt(e),tamar_flowrate_mms3(e),'LineWidth',2);
ylabel([{'River Tamar flowrate (m^{3} s^{-1})'}],'FontSize',13);
xlabel('Time (month/day)','FontSize',13);
dynamicDateTicks([], [], 'mm/dd');
set(gca,'FontSize',13);
set(gca,'FontSize',13);
xlim([klp jkl])
ylim([0 max(tamar_flowrate_mms3(e))])

load('C:\Users\rps207\Documents\MATLAB\2016 - Temperature mooring output analysis\Data\NSTempmooringL4.mat')

[ a b]=find(T_DT>klp &T_DT<jkl);

figure(204)
degree_symbol= sprintf('%c', char(176));
micro_symbol= sprintf('%c', char(0181));
plot(T_DT(a),T_35m(a),'-m'); hold on;
plot(T_DT(a),T_star_29m(a),'-y'); hold on;
plot(T_DT(a),T_24m(a),'-k'); hold on;
plot(T_DT(a),T_15m(a),'-g'); hold on;
plot(T_DT(a),T_06m(a),'-b'); hold on;
plot(T_DT(a),T_03m(a),'-r'); hold on;
dynamicDateTicks([], [], 'dd/mm');
legend('3.5m','2.9m','2.4m','1.5m','0.6m','0.3m','FontSize',13);
xlabel('Time','FontSize',13);
ylabel({['Temperature (',num2str(degree_symbol),'C)']},'FontSize',13);
set(gca,'FontSize',13);
set(gca,'FontSize',13);
clearvars a  T_03m T_06m T_15m T_24m T_35m T_DT T_star_29m



load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\ph_L4.mat'],'SST_phpack','pH_1_phpack','dt_phpack','SAL_phpack');

[~ , Z]=unique(L4_dt(c));
L4DT_SS=L4_dt(c);
L4Sal_SS=L4Salinity(c)
%this fig shows that the two agree well enough could probably
% use either but use L4 one as this is what has been presented in the paper
% already
% figure(1234567) 
% plot(L4Sal_SS)
% hold on
% plot(SAL_phpack,'r')

%INTERP MOORING SALINITY
SAL_interp_ph= interp1(L4DT_SS(Z),L4Sal_SS(Z) ,dt_phpack); %interpolate both temperature probes so they are the same lemngth 
TA_L4_sal_calc=829 + SAL_interp_ph*42.9; %L4 surface relationship TA/SAL from KITIDIS 2012



[bv , bt]=find(underway_dist'<1);
A=[]%pre allocate matrix
for i=1:length(pH_1_phpack)
A(i,:)=CO2SYS(TA_L4_sal_calc(i),pH_1_phpack(i),1,3,SAL_phpack(i),20,SST_phpack(i),0,1,2,0.5,1,4,1);
end
fco2_calc_ta_ph=A(:,20);


figure(205)%plot Helens pH data and CO2 from the quest and CO2
subplot(1,2,1)
plot(dt_phpack,pH_1_phpack,'o','LineWidth',0.5);
hold on
addaxis(dt_phpack,TA_L4_sal_calc,'-o','LineWidth',2)
dynamicDateTicks([], [], ' dd/mm');
subplot(1,2,2)
plot(underway_DT(bt),underwayfCO2_sw(bt),'-o','LineWidth',2)
hold on
plot(dt_phpack,fco2_calc_ta_ph,'-r','LineWidth',2)
dynamicDateTicks([], [], ' dd/mm');


    %quest data within 1km of L4
    cHeader = {'year','month','day','hour','minute','second','pCO2 seawater','pCO2 atmosphere','fCO2 seawater','fCO2 atmosphere', 'Salinity','Temperature','Latitude','Longitude'}; %HEADERS WITH UNITS
    commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
    commaHeader = commaHeader(:)';
    textHeader = cell2mat(commaHeader); %cHeader in text with commas
    
    dt_str=(underway_DT(bt))
    dt_str_datvec=datevec(underway_DT(bt));
    fid = fopen(['Quest_2016_fCO2_L4_only.csv'],'w');
    fprintf(fid,'%s\n',textHeader);
    fclose(fid)
    %write data to end of file
    dlmwrite(['Quest_2016_fCO2_L4_only.csv'],[datevec(dt_str) underwaypCO2_sw(bt) underwaypCO2_atm(bt) underwayfCO2_sw(bt) underwayfCO2_atm(bt)  underwayS_sea(bt) underwayT_sea(bt) underwayLat(bt) underwayLong(bt)],'-append');


    %ALL Quest data
    cHeader = {'year','month','day','hour','minute','second','pCO2 seawater','pCO2 atmosphere','fCO2 seawater','fCO2 atmosphere', 'Salinity','Temperature','Latitude','Longitude'}; %HEADERS WITH UNITS
    commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
    commaHeader = commaHeader(:)';
    textHeader = cell2mat(commaHeader); %cHeader in text with commas
    
    dt_str=(underway_DT)
    fid = fopen(['Quest_2016_fCO2_all.csv'],'w');
    fprintf(fid,'%s\n',textHeader);
    fclose(fid)
    %write data to end of file
    dlmwrite(['Quest_2016_fCO2_all.csv'],[datevec(dt_str) underwaypCO2_sw underwaypCO2_atm underwayfCO2_sw underwayfCO2_atm  underwayS_sea underwayT_sea underwayLat underwayLong],'-append');


    figure(114)
    counterF=1;subsample=30;Subplotcolumns=4;
    for j = 1:length(fList(:,1));
        subplot(3,Subplotcolumns,counterF)
        %skip these
        %unique string sequence all days when want to compare NSOP data
        if strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-06-30')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-08-04')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-15')|strcmpi(fList(j,1:10),'2016-09-21')
            
            load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'underwayfCO2_sw_interp','dt_since_lowtide_und_ret','indt','Licor_fco2_surface_lag','inds','SeaTemp_interp','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','Salinity_interp','undind1','undind2','undind3','undind4','ind','Salinity_interp','Licor_Datetime','dist','dt_since_lowtide_und_retHHMM','avgdepth');
%             for ds=1:length(inds)-1
%                 x=[indt(ds) inds(ds) inds(ds) indt(ds)];
%                 y= [ max(Licor_fco2_surface_lag) max(Licor_fco2_surface_lag) 0 0];
%                 patch(Licor_Datetime(x),y,colour_greyshade,'LineStyle','none');
%             end

                x=[indt(end) inds(1) inds(1) indt(end)];
                y= [ max(Licor_fco2_surface_lag) max(Licor_fco2_surface_lag) 0 0];
                patch(Licor_Datetime(x),y,colour_peachback,'LineStyle','none');

            hold on
            
            %return journey
            %         if strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
            if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
            else
                %Underway return
                x=[undind4 undind3 undind3 undind4];
                y= [ max(Licor_fco2_surface_lag) max(Licor_fco2_surface_lag) 0 0];
                patch(Licor_Datetime(x),y,colour_verylightblue);
                %plot NSOP system return
                plot(Licor_Datetime(undind3:undind4),Licor_fco2_combined(undind3:undind4),'k','LineWidth',2);

            end
            %outward journey
            %         if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-21')
            if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')
                
            else
                %underway out
                x=[undind2 undind1 undind1 undind2];
                y= [ max(Licor_fco2_surface_lag) max(Licor_fco2_surface_lag) 0 0];
                patch(Licor_Datetime(x),y,colour_verylightblue);
                %plot NSOP system out
                plot(Licor_Datetime(undind1:undind2),Licor_fco2_combined(undind1:undind2),'k','LineWidth',2);

            end
            counterF=counterF+1;
            set(gca,'FontSize',12)
            dynamicDateTicks([], [], 'HH:MM');
            %add text labels
            day=str2double(fList(j,9:10));
            month=str2double(fList(j,6:7));
            transectdate=datenum(2016,month,day,1,1,1)
            dateoutformatfullmonth = 'mmmm ';
            monthfull=datestr(transectdate,dateoutformatfullmonth)
            dayordinal=num2ordinal(day)
            title({[monthfull dayordinal]});
            ylabel(['fCO_2' '(',num2str(micro_symbol),'atm)'],'FontSize',12);
            set(gca,'FontSize',12);
            
            %         datetick('x','HH:MM', 'keepticks')
            %plot showerhead data
%             plot(Licor_Datetime,underwayfCO2_sw_interp,'g','LineWidth',1);
            

%             %plot line as grey
%             for p=1:length(inds)
%                 plot(Licor_Datetime(inds(p):indt(p)),Licor_fco2_combined(inds(p):indt(p)),'color',colour_greyshade,'LineWidth',2);
%             end
            
            %black data that is 2.5-3.5m
            for p=1:length(inds)
                if avgdepth(p)>2.5  & avgdepth(p)<3.5 
                plot(Licor_Datetime(inds(p):indt(p)),Licor_fco2_combined(inds(p):indt(p)),'k','LineWidth',2);
                else
                end
            end
             
            plot(underway_DT,underwayfCO2_sw,'*','color',colour_crimson,'MarkerSize',14,'LineWidth',2);

            if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')  
            else
                x=[inds(1) undind2 undind2 inds(1)];
                y= [ max(Licor_fco2_surface_lag) max(Licor_fco2_surface_lag) 0 0];
                patch(Licor_Datetime(x),y,'w','LineStyle','none');
            end
            
            if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
            else
                x=[undind3 indt(end) indt(end) undind3];
                y= [ max(Licor_fco2_surface_lag) max(Licor_fco2_surface_lag) 0 0];
                patch(Licor_Datetime(x),y,'w','LineStyle','none');
            
            end
            
            
            %set y limits for each subplot
            if strcmpi(fList(j,1:10),'2016-06-10');
                ylim([343 365]);
                set(gca, 'XTick', datenum(2016,06,10,12,00,0):1/48:datenum(2016,06,10,13,00,0))
                rotateXLabels( gca(), 90 )
                datetick('x','HH:MM', 'keepticks')
                                xlim([datenum(2016,06,10,11,45,0) datenum(2016,06,10,13,15,0) ])

            elseif strcmpi(fList(j,1:10),'2016-06-15');
                ylim([330 361]);
                set(gca, 'XTick', datenum(2016,06,15,11,30,0):1/48:datenum(2016,06,15,12,30,0))
                rotateXLabels( gca(), 90 )
                datetick('x','HH:MM', 'keepticks')
                                xlim([datenum(2016,06,15,11,30,0) datenum(2016,06,15,12,37,0) ])

            elseif strcmpi(fList(j,1:10),'2016-06-22');
                ylim([340 386]);
                set(gca, 'XTick', datenum(2016,06,22,10,30,0):1/48:datenum(2016,06,22,14,00,0))
                rotateXLabels( gca(), 90 )
                datetick('x','HH:MM', 'keepticks')
                                xlim([datenum(2016,06,22,10,20,0) datenum(2016,06,22,14,00,0) ])

            elseif strcmpi(fList(j,1:10),'2016-06-30');
                ylim([340 380]);
                set(gca, 'XTick', datenum(2016,06,30,10,30,0):1/48:datenum(2016,06,30,14,00,0))
                rotateXLabels( gca(), 90 )
                datetick('x','HH:MM', 'keepticks')
                                xlim([datenum(2016,06,30,10,29,0) datenum(2016,06,30,14,15,0) ])

            elseif strcmpi(fList(j,1:10),'2016-07-07');
                ylim([340 450]);
                set(gca, 'XTick', datenum(2016,07,07,10,30,30):1/48:datenum(2016,07,07,15,30,30))
                rotateXLabels( gca(), 90 )
                datetick('x','HH:MM', 'keepticks')
                                xlim([datenum(2016,07,07,10,26,30) datenum(2016,07,07,15,31,30) ])

            elseif strcmpi(fList(j,1:10),'2016-07-13');
                ylim([350 400]);
                set(gca, 'XTick', datenum(2016,07,13,11,30,00):1/48:datenum(2016,07,13,14,30,00))
                rotateXLabels( gca(), 90 )
                datetick('x','HH:MM', 'keepticks')
                                xlim([datenum(2016,07,13,11,26,00) datenum(2016,07,13,14,39,00) ])

            elseif strcmpi(fList(j,1:10),'2016-07-20');
                ylim([315 370]);
                set(gca, 'XTick', datenum(2016,07,20,12,00,00):1/48:datenum(2016,07,20,14,30,00))
                rotateXLabels( gca(), 90 )
                datetick('x','HH:MM', 'keepticks')
                                xlim([datenum(2016,07,20,11,49,00) datenum(2016,07,20,14,35,00) ])

            elseif strcmpi(fList(j,1:10),'2016-08-04');
                ylim([335 445]);
                set(gca, 'XTick', datenum(2016,08,04,11,00,00):1/48:datenum(2016,08,04,14,30,00))
                rotateXLabels( gca(), 90 )
                datetick('x','HH:MM', 'keepticks')
                                xlim([datenum(2016,08,04,10,54,00) datenum(2016,08,04,14,39,00) ])

            elseif strcmpi(fList(j,1:10),'2016-08-10');
                ylim([330 415]);
                
            elseif strcmpi(fList(j,1:10),'2016-08-17');
                ylim([285 415]);
                set(gca, 'XTick', datenum(2016,08,17,11,30,00):1/48:datenum(2016,08,17,13,00,00))
                rotateXLabels( gca(), 90 )
                datetick('x','HH:MM', 'keepticks')
                                xlim([datenum(2016,08,17,11,04,00) datenum(2016,08,17,13,28,00) ])

            elseif strcmpi(fList(j,1:10),'2016-08-24');
                ylim([375 412]);
                set(gca, 'XTick', datenum(2016,08,24,10,00,00):1/48:datenum(2016,08,24,11,00,00))
                rotateXLabels( gca(), 90 )
                datetick('x','HH:MM', 'keepticks')
                                xlim([datenum(2016,08,24,09,55,00) datenum(2016,08,24,11,10,00) ])

            elseif strcmpi(fList(j,1:10),'2016-09-02');
                ylim([380 415]);               
                
            elseif strcmpi(fList(j,1:10),'2016-09-15');
                ylim([368 445]);
                set(gca, 'XTick', datenum(2016,09,15,10,00,00):1/48:datenum(2016,09,15,14,00,00))
                rotateXLabels( gca(), 90 )
                datetick('x','HH:MM', 'keepticks')
                                xlim([datenum(2016,09,15,09,57,00) datenum(2016,09,15,14,09,00) ])

            elseif strcmpi(fList(j,1:10),'2016-09-21');
                ylim([382 460])
                set(gca, 'XTick', datenum(2016,09,21,09,00,00):1/48:datenum(2016,09,21,12,30,00))
                rotateXLabels( gca(), 90 )
                datetick('x','HH:MM', 'keepticks')
                                xlim([datenum(2016,09,21,08,31,00) datenum(2016,09,21,12,54,00) ])

                
            else
            end
        else
        end
    end
    
        
    
 %repeat code as above sort fewer transects
    
    Sorttimes=[];Inout=[];fListtransects=[];
      %find all transects and their times
    for j = 1:length(fList(:,1));
        %find all data on journeys back from L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-06-30')|strcmpi(fList(j,1:10),'2016-07-07')
       
            load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'dt_since_lowtide_und_ret','SeaTemp_interp','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','Salinity_interp','undind1','undind2','undind3','undind4','ind','Salinity_interp','Licor_Datetime','dist','dt_since_lowtide_und_retHHMM');
            Sorttimes = [Sorttimes ; dt_since_lowtide_und_retHHMM];
            Inout = [Inout ; 1];
            fListtransects = [fListtransects ; fList(j,1:10)];
        else
        end
    end

    numericalSorttimes=datenum(Sorttimes,'HH:MM');
    [p subplotnumorder]=sort(numericalSorttimes');
    [p  spnum2]=sort(subplotnumorder);
    
    
    %plot fluxes for 4 timeframe
    figure(805)
    subsample=30;counterF=1;
    hold on
    for j = 1:length(fList(:,1));
        if strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-06-30')|strcmpi(fList(j,1:10),'2016-07-07')
            
            load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'Longitude_interp','Latitude_interp','undind1','undind2','undind3','undind4','Licor_fco2_combined','SeaTemp_interp','Salinity_interp','Licor_Datetime','dist','Fmolhrx','dt_since_lowtide_und_retHHMM','Flux_showerhead','showerlong_overlap_und_only','showerlat_overlap_und_only');
            subsample=30;
            Subplotcolumns=2;
            subplot(2,Subplotcolumns,spnum2(counterF))
            m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
            m_grid('linestyle','none','tickdir','out',  'tickstyle','dd','xtick',[-4.25:0.05:-4.1],'ytick',[50.25:0.02:50.35],'fontsize',13,'%.3f')
            m_patch(coast(:,1),coast(:,2),[0.7412    0.7176    0.4196],'EdgeColor','none'); hold on;
%             text(0.02,1.08,'(a)','color','k','Fontsize',13,'Fontweight','bold','BackgroundColor','w','units','normalized');
            m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,Fmolhrx(undind3:subsample:undind4),'filled');%-Flux_showerhead
            
            m_scatter(showerlong_overlap_und_only,showerlat_overlap_und_only,50,Flux_showerhead,'d','filled','MarkerEdgeColor','m','LineWidth',1.5);%-Flux_showerhead

            y=ylabel('Latitude');
            set(y,'Units','Normalized','Position',[-0.2,0.5,0]);
            x=xlabel('Longitude');
            set(x,'Units','Normalized','Position',[0.5,-0.05,0]);
            set(gca,'fontsize',13);
            set(gca,'fontsize',13);
            [C,D]=m_ll2xy(-4.217,50.25);    %add L4
            line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
            [CC,D]=m_ll2xy(-4.217,50.25);
            text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',12);
            [C,D]=m_ll2xy(-4.188801,50.318280);       %add penlee
            line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
            [CC,D]=m_ll2xy(-4.188801,50.318280);
            text(CC,D,'PPAO','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',12);
            [C,D]=m_ll2xy(-4.158932,50.334499);    %add breakwater
            line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
            [CC,D]=m_ll2xy(-4.158932,50.334499);
            text(CC,D,'Breakwater','HorizontalAlignment', 'left','VerticalAlignment', 'top','fontsize',12);
            %     ring1=m_range_ring(-4.217,50.25,0.5); set(ring1,'color','r')ring2=m_range_ring(-4.217,50.25,1); set(ring2,'color','b');
            hj=colorbar;
            ylabel(hj,['Sea' char(8212) 'Air flux of CO_2 (mmol m^2 hr^-1)'],'FontSize',13);
            title({[fList(j,9:10) '/' fList(j,6:7) ' LW+' dt_since_lowtide_und_retHHMM(1:5) 'hrs']},'FontSize',16);hold on
            hold off
            counterF=counterF+1;
        else
        end
    end


    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%stats for paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% stats for paper
%difference between the two systems!
mean(Shower_overlap-NSOP_overlap_) %2.77
std(Shower_overlap([1:11 13:end])-NSOP_overlap_([1:11 13:end])) %6.25
RMSE = sqrt(mean(((Shower_overlap-NSOP_overlap_).^2))) 


mean(Shower_overlap_und([1:7 9:10 12:end])-NSOP_overlap_und_([1:7 9:10 12:end])) %2.68
std(Shower_overlap_und([1:7 9:10 12:end])-NSOP_overlap_und_([1:7 9:10 12:end])) %28.15
RMSE = sqrt(mean(((Shower_overlap_und([1:7 9:10 12:end])-NSOP_overlap_und_([1:7 9:10 12:end])).^2))) %27.12


%paper stat L4 salinity-  'this is supported by the record of surface salinity
%from the L4 mooring which had a value of 35.15 ±0.08  during this period '
[v,k]=find(L4_dt<klp)
ts_start=max(v)
[v,k]=find(L4_dt>jkl)
ts_end=min(v)
mean(L4Salinity(ts_start:ts_end))
std(L4Salinity(ts_start:ts_end))

%paper stat range of salinity peaks - 'also present in the salinity record for L4 were intermittent short lived pulses in salinity of -0.3 – -0.4 .
    sF = 3600; %sampling frequency in Hz
    spec = hspec_ss(L4Salinity(ts_start:ts_end),L4Salinity(ts_start:ts_end),sF);
    figure(99); loglog(spec(:,1),spec(:,1).*spec(:,2),'.'); %plot the frequency response and use it to set cutoff frequency
    N=2; Wn = 1500/(sF/2);  % Normalized cutoff frequency
    [B,A] = butter(N,Wn);   %butterworth filter
    h=fvtool(B,A);          %plot magnitude response
    smoothL4Salinity = filtfilt(B,A,L4Salinity(ts_start:ts_end)); %smooth data using filtfilt
    diffL4sal=diff(smoothL4Salinity); %find differential
%     plot(diffL4sal)
    threshold=0.1;
    [pk, ish] = findpeaks(diffL4sal,'MINPEAKHEIGHT',threshold); %set threshhold and find peaks
%     [pk2, ind2] = findpeaks(-diffL4sal,'MINPEAKHEIGHT',threshold);%set threshhold and find peaks
%     ind=([ind1;ind2]); %concatinate index matrices
    salind=([ish]); %concatinate index matrices
    L4_dtind=L4_dt(ts_start:ts_end);
    
    %check this has worked with plot
    figure(999);
    plot(L4_dt(c),L4Salinity(c),'LineWidth',2);
    hold on
    scatter(Und_DT(b),Salinity(b),4,'o','filled');
    dynamicDateTicks([], [], 'mm/dd');
    xlabel('Time (month/day)','FontSize',13);
    ylabel({['Salinity']},'FontSize',13);
    ylim([34 max(L4Salinity)])
    klp=datenum(2016,06,10,0,0,1)
    jkl=datenum(2016,09,30,0,0,1)
    xlim([klp jkl])
    hold on
    plot(L4_dtind(salind),35*ones(length(salind),1),'r*')
    dynamicDateTicks([], [], 'mm/dd');
    
    %found peaks
    rowsToDelete = salind < 12 | salind > length(L4Salinity)-12;
    salind(rowsToDelete) = [];%remove first point nuisance.
    L4Salinitypeaks=L4Salinity(ts_start:ts_end+2); %data subset
    for k=1:length(salind)
        dailysalavgbeforefresh(k)=max(L4Salinitypeaks(salind(k)-12:salind(k)+12));%look for the maximum that day
    end
    %work out difference in salinity
    salinitydrops=dailysalavgbeforefresh'-L4Salinitypeaks(salind);
    figure(999)
    plot(L4_dtind(salind),L4Salinitypeaks(salind),'k*')
    plot(L4_dtind(salind),L4Salinitypeaks(salind)+salinitydrops,'g*')
    %this totally worked!
    %stat for paper
    mean(salinitydrops)% 0.3103
    max(salinitydrops)% 0.68
    min(salinitydrops)% 0.04
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %get all the data from all 15 transects!
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Salinity_trans=[];
    DistL4_trans = [];%
    Datetime_trans = [];%
    Licor_fco2_trans = [];%
    Licor_fco2_atm_trans = [];%
    SeaTemp_trans = [];%
    Latitude_trans = [];%
    Longitude_trans = [];%
    del_Salinity_trans=[];
    del_co2_trans=[];
    for j = 1:length(fList(:,1)); 
 %plot all data on journeys back from L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
else    hold on
load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'undind1','undind2','undind3','undind4','Latitude_interp','Longitude_interp','Salinity_interp','SeaTemp_interp','Licor_fco2_combined','Licor_Datetime','dist','inds','indt');
        Salinity_trans = [Salinity_trans ; Salinity_interp(undind3:undind4)];%
        DistL4_trans = [DistL4_trans ; dist(undind3:undind4)];%
        Datetime_trans = [Datetime_trans ; Licor_Datetime(undind3:undind4)];%
         Licor_fco2_trans = [Licor_fco2_trans ; Licor_fco2_combined(undind3:undind4)];%
         Licor_fco2_atm_trans = [Licor_fco2_atm_trans ; Licor_fco2_combined(undind3:undind4)];%
         SeaTemp_trans = [SeaTemp_trans ; SeaTemp_interp(undind3:undind4)];%
         Latitude_trans = [Latitude_trans ; Latitude_interp(undind3:undind4)];%
         Longitude_trans = [Longitude_trans ; Longitude_interp(undind3:undind4)];%

         s = struct(char({['temp' fList(j,[1:4 6:7 9:10])]}),SeaTemp_interp(undind3:undind4));
        del_Salinity_trans = [del_Salinity_trans ; Salinity_interp(undind3:undind4)-Salinity_interp(undind3)];%
        del_co2_trans = [del_co2_trans ; Licor_fco2_combined(undind3:undind4)-Licor_fco2_combined(undind3)];%


hold off
end
    end
for j = 1:length(fList(:,1)); 
%plot all data on journeys out to L4 , if listed below ignore, wrong or didnt work for reasons above
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')
else    hold on
load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'undind1','undind2','undind3','undind4','Salinity_interp','Salinity_interp','SeaTemp_interp','Latitude_interp','Longitude_interp','Licor_Datetime','dist','inds','indt');
%skip if underway indexes are empty
        Salinity_trans = [Salinity_trans ; Salinity_interp(undind1:undind2)];%
        DistL4_trans = [DistL4_trans ; dist(undind1:undind2)];%
        Datetime_trans = [Datetime_trans ; Licor_Datetime(undind1:undind2)];%
        Licor_fco2_trans = [Licor_fco2_trans ; Licor_fco2_combined(undind1:undind2)];%
        SeaTemp_trans = [SeaTemp_trans ; SeaTemp_interp(undind1:undind2)];%
        Latitude_trans = [Latitude_trans ; Latitude_interp(undind1:undind2)];%
        Longitude_trans = [Longitude_trans ; Longitude_interp(undind1:undind2)];%hold off
        del_Salinity_trans = [del_Salinity_trans ; Salinity_interp(undind1:undind2)-Salinity_interp(undind2)];%
        del_co2_trans = [del_co2_trans ; Licor_fco2_combined(undind1:undind2)-Licor_fco2_combined(undind2)];%

end
end

%export high frequency data

[B,I] = sort(Datetime_trans,'ascend') 


    %write header to file
    cHeader = {'year','month','day','hour','minute','second','fCO2 seawater (micro atm)', 'Salinity (PSU)','Temperature(degrees C)','Latitude(deg N)','Longitude (deg W)'}; %HEADERS WITH UNITS
    commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
    commaHeader = commaHeader(:)';
    textHeader = cell2mat(commaHeader); %cHeader in text with commas
    
    dt_str=(Datetime_trans(I))
    fid = fopen(['Sims_2016_highfrequency_WCO_fCO2.csv'],'w');
    fprintf(fid,'%s\n',textHeader);
    fclose(fid)
    %write data to end of file
    dlmwrite(['Sims_2016_highfrequency_WCO_fCO2.csv'],[datevec(dt_str) Salinity_trans(I) Licor_fco2_trans(I) SeaTemp_trans(I) Latitude_trans(I) Longitude_trans(I)],'-append');

    
    

    
    
    

max(Salinity_trans) %35.23
min(Salinity_trans) %34.17

dist_ppao=[]

%find range of salinity close to PPAO
     for t=1:length(Latitude_trans)
        dist_ppao(t,:)= pos2dist(50.318, -4.189,Latitude_trans(t),-1*Longitude_trans(t),2);
     end
   [c , ]=find(dist_ppao<0.5);
   

dist_break=[]
%find range of salinity close to PPAO
     for t=1:length(Latitude_trans)
        dist_break(t,:)= pos2dist(50.334499, -4.158932,Latitude_trans(t),-1*Longitude_trans(t),2);
     end   
   
   
      [d , ]=find(dist_break<4);


   min(Salinity_trans(c)) %34.51
   mean(Salinity_trans(c)) % 34.84
   
   mean(Licor_fco2_trans(d)) % 385.56
   min(Licor_fco2_trans(d)) % 338.03
   max(Licor_fco2_trans(d)) % 440.47

   
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %get differences between start and end of profiles
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    Salinitydel_trans =[];
    date_trans=[];
    SeaTempdel_trans=[];
    SeaTemprange_trans=[]; 
    avgtempl4_trans=[];
    fCO2del_trans=[];
    
    for j = 1:length(fList(:,1));
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
        else    hold on
            load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'undind1','undind2','undind3','undind4','Latitude_interp','Longitude_interp','Salinity_interp','SeaTemp_interp','Licor_fco2_combined','Licor_Datetime','dist','inds','indt');
            Salinitydel_trans = [Salinitydel_trans ; Salinity_interp(undind3)-Salinity_interp(undind4)];%
            SeaTempdel_trans = [SeaTempdel_trans ; SeaTemp_interp(undind3)-SeaTemp_interp(undind4)];%
            SeaTemprange_trans = [SeaTemprange_trans ; max(SeaTemp_interp(undind3:undind4))-min(SeaTemp_interp(undind3:undind4))];%
            date_trans = [date_trans ; fList(j,1:10)];%
            fCO2del_trans = [fCO2del_trans ; Licor_fco2_combined(undind3)-Licor_fco2_combined(undind4)];%

        end
    end
    for j = 1:length(fList(:,1));
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')
        else
            load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'undind1','undind2','undind3','undind4','Salinity_interp','Salinity_interp','SeaTemp_interp','Latitude_interp','Longitude_interp','Licor_Datetime','dist','inds','indt');
            Salinitydel_trans = [Salinitydel_trans ; Salinity_interp(undind2)-Salinity_interp(undind1)];%
            SeaTempdel_trans = [SeaTempdel_trans ; SeaTemp_interp(undind2)-SeaTemp_interp(undind1)];%
            SeaTemprange_trans = [SeaTemprange_trans ; max(SeaTemp_interp(undind3:undind4))-min(SeaTemp_interp(undind3:undind4))];%
            date_trans = [date_trans ; fList(j,1:10)];%
            fCO2del_trans = [fCO2del_trans ; Licor_fco2_combined(undind2)-Licor_fco2_combined(undind1)];%

        end   

    end
    mean(Salinitydel_trans) %0.36
    mean(SeaTempdel_trans) %-0.15 this isn't useful because of warming 
    mean(SeaTemprange_trans) %0.58 degrees
    datevec(date_trans)
        mean(fCO2del_trans) %20.33

    %max change in co2 along transect
    fCO2del_trans(8) %61.08
    
% for this check SeaTempdel_trans and date_trans
% from paper - (0.002 , 0.039 and 0.025°C) between L4 and the breakwater on July 20th  August 10th
% and September 15th  respectively whereas other temperate differentials were
% much larger (1.18 , 0.64, 0.59 and 0.57°C) on  June 10th, 
% June 15th, June 22nd and August 4th  
    
%small temperature gradients
%july 20th - 
datevec(date_trans(6,:))
SeaTempdel_trans(6)
%sept 15 10th - 
datevec(date_trans(9,:))
SeaTempdel_trans(9)
    
%large gradients
%june10
datevec(date_trans(1,:))
SeaTempdel_trans(1)
%june 15th
datevec(date_trans(2,:))
SeaTempdel_trans(2)
%june30
datevec(date_trans(3,:))
SeaTempdel_trans(3)
%august 4th   
datevec(date_trans(7,:))
SeaTempdel_trans(7)
    
datevec(date_trans(11,:))
SeaTempdel_trans(11)

%difference in temperature from first to last transect
%4th profile june 10th
 load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(4,1:10) '.mat'],'undind1','undind2','undind3','undind4','Latitude_interp','Longitude_interp','Salinity_interp','SeaTemp_interp','Licor_fco2_combined','Licor_Datetime','dist','inds','indt');
 SeaTemp_trans_beg=SeaTemp_interp(undind3:undind4);
 mean(SeaTemp_trans_beg)% 13.90
 
 %18th profile - september 21st
 load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(18,1:10) '.mat'],'undind1','undind2','undind3','undind4','Latitude_interp','Longitude_interp','Salinity_interp','SeaTemp_interp','Licor_fco2_combined','Licor_Datetime','dist','inds','indt');
 SeaTemp_trans_end=SeaTemp_interp(undind3:undind4);
 mean(SeaTemp_trans_end)% 16.70
 
 
%weekly warming 
%plot all temperature by time?
figure(456)
plot(Datetime_trans-Datetime_trans(1),SeaTemp_trans,'*')
 c = polyfit(Datetime_trans-Datetime_trans(1),SeaTemp_trans,1);
% Display evaluated equation y = m*x + b
disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))]) %0.0365 per day
 

    
%difference in CO2 from first to last transect
%4th profile june 10th
 load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(4,1:10) '.mat'],'undind1','undind2','undind3','undind4','Latitude_interp','Longitude_interp','Salinity_interp','SeaTemp_interp','Licor_fco2_combined','Licor_Datetime','dist','inds','indt');
 Licor_fco2_combined_trans_beg=Licor_fco2_combined(undind3:undind4);
 mean(Licor_fco2_combined_trans_beg)% 354.46
 
 %18th profile - september 21st
 load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(18,1:10) '.mat'],'undind1','undind2','undind3','undind4','Latitude_interp','Longitude_interp','Salinity_interp','SeaTemp_interp','Licor_fco2_combined','Licor_Datetime','dist','inds','indt');
 Licor_fco2_combined_trans_end=Licor_fco2_combined(undind3:undind4);
 mean(Licor_fco2_combined_trans_end)% 420.75

 
 %co2 variability around penlee on certain profiles
%4th profile july 7th
 load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(8,1:10) '.mat'],'undind1','undind2','undind3','undind4','Latitude_interp','Longitude_interp','Salinity_interp','SeaTemp_interp','Licor_fco2_combined','Licor_Datetime','dist','inds','indt','ATMSCO2meantransmatrix');
 dist_ppao=[]
 x=Latitude_interp(undind3:undind4)
 y=Longitude_interp(undind3:undind4)
%find range of salinity close to PPAO
     for t=1:length(x)
        dist_ppao(t,:)= pos2dist(50.318, -4.189,x(t),-1*y(t),2);
     end
     Licor_fco2_combined_trans=Licor_fco2_combined(undind3:undind4);
    Sal_trans=Salinity_interp(undind3:undind4);
    
    [c , ]=find(dist_ppao<2);
   min(Licor_fco2_combined_trans(c))
      max(Licor_fco2_combined_trans(c))
      
       [d , ]=find(dist_ppao<1);
       min(Sal_trans(d)) %34.5

  max(Salinity_interp(undind3:undind4))
  min(Salinity_interp(undind3:undind4))

  datevec(date_trans(4,:))
fCO2del_trans(4) %26.43
  
     min(Licor_fco2_combined_trans)%360.98
      max(Licor_fco2_combined_trans)%434.42

      max(Licor_fco2_combined_trans)-min(Licor_fco2_combined_trans) %73.43
ATMSCO2meantransmatrix(1)



%co2 variability around penlee on certain profiles
%5th profile june 15th
 load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(5,1:10) '.mat'],'undind1','undind2','undind3','undind4','Latitude_interp','Longitude_interp','Salinity_interp','SeaTemp_interp','Licor_fco2_combined','Licor_Datetime','dist','inds','indt','ATMSCO2meantransmatrix');
  max(Salinity_interp(undind3:undind4)) %35.09
  min(Salinity_interp(undind3:undind4)) %34.61

dist_break=[]
 x=Latitude_interp(undind3:undind4)
 y=Longitude_interp(undind3:undind4)
      Licor_fco2_combined_trans=Licor_fco2_combined(undind3:undind4);
    Sal_trans=Salinity_interp(undind3:undind4);
%find range of salinity close to PPAO
     for t=1:length(x)
        dist_break(t,:)= pos2dist(50.334499, -4.158932,x(t),-1*y(t),2);
     end

   [c , ]=find(dist_break<2);
   min(Sal_trans(c))
      max(Sal_trans(c))

%find range of salinity close to PPAO
     for t=1:length(x)
        dist_ppao(t,:)= pos2dist(50.318, -4.189,x(t),-1*y(t),2);
     end
   [d , ]=find(dist_ppao<1);
      mean(Sal_trans(d))
      mean(Licor_fco2_combined_trans(d))%353.02 

  datevec(date_trans(2,:))
fCO2del_trans(2) %18.02




%co2 variability around penlee on certain profiles
%7th profile june 30th
 load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(7,1:10) '.mat'],'undind1','undind2','undind3','undind4','Latitude_interp','Longitude_interp','Salinity_interp','SeaTemp_interp','Licor_fco2_combined','Licor_Datetime','dist','inds','indt','ATMSCO2meantransmatrix');
  max(Salinity_interp(undind3:undind4)) %35.09
  min(Salinity_interp(undind3:undind4)) %34.17

  max(SeaTemp_interp(undind3:undind4)) %14.83
  max(Licor_fco2_combined(undind3:undind4)) %385.51

  %gradient just to penleee
x=(Licor_fco2_combined(undind3:undind4));
z=Latitude_interp(undind3:undind4);
[h n]=  find(z<50.305);
v=x(h);
  v(1)-v(end) %12.18
  
  %over whole transect large 37.78
  datevec(date_trans(3,:))
fCO2del_trans(3) %18.02
  
  
  
  
  
  
  
  
%salinity vs co2 relationship tamar bottles
x=Tamar_fco2_calc_ta_dic(3:end-5)-429.93;
y=Tamar_Salinity(3:end-5)-35.16;
c = polyfit(y,x,1);%slope is 59.74

 
%salinity vs co2 relationship tamar bottles T3-T7
x=Tamar_fco2_calc_ta_dic(3:7);
y=Tamar_Salinity(3:7);
c = polyfit(y,x,1);%slope is 59.74
 
 
 
%same code as used for figure above, save matrix together to derive
%equation, also exclude 07/07
del_Salinity_trans = []
    del_co2_trans=[]

for j = 1:length(fList(:,1));
        %plot all data on journeys back from L4 , if listed below ignore, wrong or didnt work for reasons above
% if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02');
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
    
else
    hold on
    load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'Licor_fco2_combined','undind1','undind2','undind3','undind4','Salinity_interp','Licor_Datetime','dist');
    hold off
        del_Salinity_trans = [del_Salinity_trans ; Salinity_interp(undind3:undind4)-Salinity_interp(undind3)];%
    del_co2_trans = [del_co2_trans ; Licor_fco2_combined(undind3:undind4)-Licor_fco2_combined(undind3)];%
        end
    %plot all data on journeys out to L4 , if listed below ignore, wrong or didnt work for reasons above
% if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-04')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')|strcmpi(fList(j,1:10),'2016-07-07');
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')
else
    hold on
    load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'Licor_fco2_combined','undind1','undind2','undind3','undind4','Salinity_interp','Licor_Datetime','dist');
   
    hold off
    del_Salinity_trans = [del_Salinity_trans ; Salinity_interp(undind1:undind2)-Salinity_interp(undind2)];%
    del_co2_trans = [del_co2_trans ; Licor_fco2_combined(undind1:undind2)-Licor_fco2_combined(undind2)];%
        end    
end

 %for transect data get equation for sal vs fco2
c = polyfit(del_Salinity_trans,del_co2_trans,1)%-39.8273    5.4980



%get drop in co2 on station from sept 21st


        load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(18,1:10) '.mat'],'underwayfCO2_sw_interp','dt_since_lowtide_und_ret','indt','Licor_fco2_surface_lag','inds','SeaTemp_interp','Longitude_interp','Latitude_interp','dt_since_lowtide_und_ret','Licor_fco2_combined','Salinity_interp','undind1','undind2','undind3','undind4','ind','Salinity_interp','Licor_Datetime','dist','dt_since_lowtide_und_retHHMM');
allco2nsop=[]
meanbinco2nsop=[]
        for h=1:length(inds)
       allco2nsop=[ allco2nsop;(Licor_fco2_combined(inds(h):indt(h)))];
        meanbinco2nsop=[ meanbinco2nsop;mean(Licor_fco2_combined(inds(h):indt(h)))];

        end
max(allco2nsop) %458.98
min(allco2nsop) %383.32
        458.98-383.32
        
        meanbinco2nsop(16)%431.75
                meanbinco2nsop(26) %393.84
431.75-393.84

                
                
     %ALL underway data
    cHeader = {'year','month','day','hour','minute','second', 'Salinity','Temperature','Chlorophyll','Pressure','Windspeed (m/s)','Latitude','Longitude'}; %HEADERS WITH UNITS
    commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
    commaHeader = commaHeader(:)';
    textHeader = cell2mat(commaHeader); %cHeader in text with commas
    
    dt_str=(Und_DT)
    fid = fopen(['Quest_underway_2016.csv'],'w');
    fprintf(fid,'%s\n',textHeader);
    fclose(fid)
    %write data to end of file
    dlmwrite(['Quest_underway_2016.csv'],[datevec(dt_str) Salinity SeaTemp Chla Pressure Windspeedabs Latitude Longitude],'-append');
               
                
    %l4 mooring data as 1 text file
    cHeader = {'year','month','day','hour','minute','second', 'Salinity','Temperature','Windspeed (m/s)'}; %HEADERS WITH UNITS
    commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
    commaHeader = commaHeader(:)';
    textHeader = cell2mat(commaHeader); %cHeader in text with commas
    
    dt_str=(L4_dt)
    fid = fopen(['L4_mooring_2016.csv'],'w');
    fprintf(fid,'%s\n',textHeader);
    fclose(fid)
    %write data to end of file
    dlmwrite(['L4_mooring_2016.csv'],[datevec(dt_str) L4Salinity L4SST L4Windspeed ],'-append');
               
c = polyfit(del_Salinity_trans,del_co2_trans,1);%-39.8273    5.4980
[correlation_coeff,pvalue_rel] = corrcoef (del_Salinity_trans(:,1),del_co2_trans(:,1));
r_sqr = power(correlation_coeff,2);




%salinity vs co2 relationship tamar bottles T3-T9
x=Tamar_fco2_calc_ta_dic(3:9);
y=Tamar_Salinity(3:9);
cc = polyfit(x,y,1);%slope is 42.24

[correlation_coeff_2,pval] = corrcoef (x,y);
r_sqr_bottles = power(correlation_coeff_2,2);



figure(1212)
%     offset_T7=Tamar_fco2_calc_ta_dic([7])-410.93;
%     offset_T7_SAL=Tamar_Salinity([7])-35.16;
    y=(10:-0.1:-25);
     x=(y*c(1))+ c(2); %+offset_T7;
     plot(x,y,'--k','LineWidth',1)
%     plot(x,y+offset_T7_SAL,'--k','LineWidth',1)
    text(400,-12,'R^{2} = 0.2165');
    text(400,-13.5,'\xifCO_2 _(_s_w_) = -39.83 \xiS + 5.498');
    hold on
    
%     y=(10:-0.1:-25);
%     x=(y*-cc(2))+ cc(1)+offset_T7;
%     plot(x,y+offset_T7_SAL,'--r','LineWidth',1)
%     xlim([ 0 700])
%     text(400,-1,'R^{2} = 0.7796','Color','r');
%     text(400,-2.5,'\xifCO_2 _(_s_w_) = -42.24 \xiS + 0.01','Color','r');%plot river Tamar bottle dat
%     
load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\Tamar_bottles.mat'],'Tamar_Date','Tamar_DICumolkg','Tamar_fco2_calc_ta_dic','Tamar_LatN','Tamar_LongW','Tamar_Salinity','Tamar_TAumolkg','Tamar_Temp1','Salinity_interp','Licor_Datetime','dist');
%on october 1st salinity was 35.15, march 11 35.18 , this is 12 days
%after(best i can do though)unless i can get it from quest underway?. 
        Tamar_distL4=[];
     for t=1:length(Tamar_LongW)
        Tamar_distL4(t,:)= pos2dist(50.25, -4.217,Tamar_LatN(t),-Tamar_LongW(t),2);
     end
     scatter(Tamar_fco2_calc_ta_dic([1:9])-410.93,Tamar_Salinity([1:9])-35.16,30,Tamar_distL4([1:9]),'d','filled')
hold on
        
    labs=text(Tamar_fco2_calc_ta_dic([1:9])-410.93+5,Tamar_Salinity([1:9])-35.16, Tamarlabels([1:9]),'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
     set(labs,'BackgroundColor', 'none');
     
    scatter(410.93-410.93,35.16-35.16,30,0,'d','filled')
    labs_l4=text(410.93-410.93+5,35.16-35.16, 'L4','color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
     set(labs_l4,'BackgroundColor', 'none');

     
    xlabel(['\xi fCO_2 _(_s_w_)(',num2str(micro_symbol),'atm)'],'fontsize',22);
    set(gca,'fontsize',22);
    set(gca,'fontsize',22);
    ylabel('\xi S (PSU)','fontsize',22);
    hj=colorbar;
    ylabel(hj,'Distance from L4 (km)','FontSize',16);
    set(gca,'fontsize',16);
    set(gca,'fontsize',16);
    hold on
    ylim([-22 1])





    
figure(1211)
subsample=30;
hold on
for j = 1:length(fList(:,1));
        %plot all data on journeys back from L4 , if listed below ignore, wrong or didnt work for reasons above
% if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02');
if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-02')
    
elseif strcmpi(fList(j,1:10),'2016-07-07')
        hold on
    load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'Licor_fco2_combined','undind1','undind2','undind3','undind4','Salinity_interp','Licor_Datetime','dist');
    scatter(Salinity_interp(undind3:undind4)-Salinity_interp(undind3),Licor_fco2_combined(undind3:undind4)-Licor_fco2_combined(undind3),8,dist(undind3:undind4),'MarkerEdgeColor',colour_greyshade)
    ylabel(['\xi fCO_2 _(_s_w_)(',num2str(micro_symbol),'atm)'],'fontsize',22)
    set(gca,'fontsize',22)
    set(gca,'fontsize',22)
    xlabel('\xi S (PSU)','fontsize',22)
    hj=colorbar;
    ylabel(hj,'Distance from L4 (km)','FontSize',16);
    set(gca,'fontsize',16)
    set(gca,'fontsize',16)
    hold off
    
else
    hold on
    load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'Licor_fco2_combined','undind1','undind2','undind3','undind4','Salinity_interp','Licor_Datetime','dist');
    scatter(Salinity_interp(undind3:undind4)-Salinity_interp(undind3),Licor_fco2_combined(undind3:undind4)-Licor_fco2_combined(undind3),8,dist(undind3:undind4))
    ylabel(['\xi fCO_2 _(_s_w_)(',num2str(micro_symbol),'atm)'],'fontsize',22)
    set(gca,'fontsize',22)
    set(gca,'fontsize',22)
    xlabel('\xi S (PSU)','fontsize',22)
    hj=colorbar;
    ylabel(hj,'Distance from L4 (km)','FontSize',16);
    set(gca,'fontsize',16)
    set(gca,'fontsize',16)
    hold off
end
     
    %plot all data on journeys out to L4 , if listed below ignore, wrong or didnt work for reasons above
% if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-04')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')|strcmpi(fList(j,1:10),'2016-07-07');
        if strcmpi(fList(j,1:10),'2016-04-27')|strcmpi(fList(j,1:10),'2016-05-11')|strcmpi(fList(j,1:10),'2016-05-26')|strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-07-27')|strcmpi(fList(j,1:10),'2016-08-10')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-09-02')|strcmpi(fList(j,1:10),'2016-09-21')

else
    hold on
    load(['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\' fList(j,1:10) '.mat'],'Licor_fco2_combined','undind1','undind2','undind3','undind4','Salinity_interp','Licor_Datetime','dist');
    scatter(Salinity_interp(undind1:undind2)-Salinity_interp(undind2),Licor_fco2_combined(undind1:undind2)-Licor_fco2_combined(undind2),8,dist(undind1:undind2))
    ylabel(['\xi fCO_2_(_s_w_)(',num2str(micro_symbol),'atm)'],'fontsize',22)
    set(gca,'fontsize',22)
    set(gca,'fontsize',22)
    xlabel('\xi Salinity (PSU)','fontsize',22)
    hj=colorbar;
    ylabel(hj,'Distance from L4 (km)','FontSize',16);
   
    set(gca,'fontsize',16)
    set(gca,'fontsize',16)
    hold off
        end  
%   plot(Tamar_fco2_calc_ta_dic([6:7,9])-Tamar_fco2_calc_ta_dic(8),Tamar_Salinity([6:7,9])-Tamar_Salinity(8),'d','MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor','r') 
%   labs=text(20+(Tamar_fco2_calc_ta_dic([6:7,9])-Tamar_fco2_calc_ta_dic(8)),(Tamar_Salinity([6:7,9])-Tamar_Salinity(8)), Tamarlabels([6:7,9]),'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w');   
end
    hold on
    x=0.2:-0.01:-1;
    y=(x*c(1))+ c(2);
    plot(x,y,'--k','LineWidth',3);
    text(-0.8,-45,'R^{2} = 0.2165');
    text(-0.9,-45,'\xifCO_2 _(_s_w_) = -39.83 \xiS + 5.50');
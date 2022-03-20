%WCO map 
clc;  clear all ; close all; %reset workspace
addpath('c:\Users\rps207\Documents\Matlab\Functions');
addpath('c:\Users\rps207\Documents\Matlab\Functions\Add_Axis');
addpath('c:\Users\rps207\Documents\Matlab\Functions\cbdate');
addpath('c:\Users\rps207\Documents\Matlab\Functions\m_map');
addpath('c:\Users\rps207\Documents\Matlab\Functions\mixing_library');
addpath('c:\Users\rps207\Documents\Matlab\Functions\cm_and_cb_utilities');
mfileDir = 'C:\Users\rps207\Documents\Matlab\CO2 NSOP output analysis\'; %path for main matlab analysis

load coast.dat %high resollution coastline for south west of uk

load('Data\bathymetry.mat')
u=longitudeceltic(longcelticsea)
i=(latitudeceltic(latcelticsea))'

Jan1_serial = datenum([2015, 1, 0, 0, 0, 0]);


fname='C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data';
cd(fname); %change current directory to folder with files want to load
fList = ls('*.mat'); 
fList=fList(30:45,:)% ste custom cruise range, 31-46 =sesonal studyl, 1:30 for cruises


%Load in underway data seasonal study

% 
% %preallocate variables
% Latitude=[];Longitude=[];
% flistsize=size(fList);
% for p=1:flistsize(1);
%     % load variable and redefine
%     S=load(fList(p,:)); 
%     Latitudex=S.Latitude_interp;
%     Longitudex=S.Longitude_interp;
%     Latitudex=mean(Latitudex);
%     Longitudex=mean(Longitudex);
%     Latitude=[Latitude,Latitudex]
%     Longitude=[Longitude,Longitudex]
% end



figure(1)
m_proj('mercator','lon',[-4.6 -4],'lat',[50 54]);  
clf
m_grid('linestyle','none','tickdir','out','fontsize',22)
% k=m_line(coast(:,1),coast(:,2));
% set(k,'color','b')
m_patch(coast(:,1),coast(:,2),'b','EdgeColor','none'); hold on;
% m_contourf(u,i,bathymetrycelticsea',[0 0]); hold on;
p=m_contour(u,i,bathymetrycelticsea',[-80:5:-15],'ShowText','on');
set(gca,'fontsize',22)
y=ylabel('Latitude');
set(y,'Units','Normalized','Position',[-0.1,0.5,0]);
set(gca,'fontsize',22)
xlabel('Longitude');
set(gca,'fontsize',22)
set(gca,'fontsize',22)
hold on
%add L4
[C,D]=m_ll2xy(-4.217,50.25);
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.217,50.25);
text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',24);
%add E1
[C,D]=m_ll2xy(-4.368,50.035);
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.368,50.035);
text(CC,D,'E1','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',24);
%add Penlee
[C,D]=m_ll2xy(-4.193056,50.318889);
line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
[CC,D]=m_ll2xy(-4.200,50.323);
text(CC,D,'Penlee Point','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',24);
%Add plymouth to  map
annotation('textbox',[0.598916666666667 0.855659397715473 0.0792083333333333 0.0477673935617861],'String',{'Plymouth'},'FontSize',24,'FitBoxToText','off','LineStyle','none');
hcb = colorbar;
cblabel('Depth','fontsize',22)
set(hcb,'fontsize',22)



figure(2)
m_proj('mercator','lon',[-5 -3.5],'lat',[49.5 50.5]);  
clf
m_contourf(u,i,bathymetrycelticsea',[0 0]); hold on;
p=m_contour(u,i,bathymetrycelticsea',[-50:5:0],'ShowText','on');
m_grid('linestyle','none','tickdir','out','fontsize',22)
% k=m_line(coast(:,1),coast(:,2));
% .set(k,'color','k')
set(gca,'fontsize',22)
y=ylabel('Latitude');
set(y,'Units','Normalized','Position',[-0.13,0.5,0]);
set(gca,'fontsize',22)
x=xlabel('Longitude');
set(gca,'fontsize',22)
set(x, 'Units', 'Normalized', 'Position', [0.5 , -0.07 , 0])
set(gca,'fontsize',22)
hold on
% %add deployment Sites to map
% h=m_scatter(Longitude([1:7 14:29]),Latitude([1:7 14:29]),22,'k');
% set(h,'marker','square','MarkerFaceColor','k');
% hcb = colorbar;
% cblabel('Depth','fontsize',22)
% set(hcb,'fontsize',22)



%Create folder and save figure
       if ~exist('Figures/Maps', 'dir')
             mkdir('Figures/Maps');
        end
%export_fig('C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Figures\Maps\DY030CruiseTrack.png'); 

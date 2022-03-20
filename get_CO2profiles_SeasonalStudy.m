%Seasonal Study workup of CO2 data
opengl('save', 'software')
clc;  clear all ; close all; %reset workspace
tstart = tic;
addpath('c:\Users\rps207\Documents\Matlab\Functions');
addpath('c:\Users\rps207\Documents\Matlab\Functions\Add_Axis');
addpath('c:\Users\rps207\Documents\Matlab\Functions\cbdate');
addpath('c:\Users\rps207\Documents\Matlab\Functions\mixing_library');
addpath('c:\Users\rps207\Documents\Matlab\Functions\despiking_tooblox');
addpath('c:\Users\rps207\Documents\Matlab\Functions\cm_and_cb_utilities');
addpath('c:\Users\rps207\Documents\Matlab\Functions\tsplot');
addpath('c:\Users\rps207\Documents\Matlab\Functions\m_map');
set(groot,'DefaultFigureColormap',jet)

addpath('C:\Users\rps207\Documents\MATLAB\VINDTA');
mfileDir = 'C:\Users\rps207\Documents\Matlab\CO2 NSOP output analysis\'; %path for main matlab analysis

degree_symbol= sprintf('%c', char(176));
micro_symbol= sprintf('%c', char(0181));

%loop through Licor files and 
path = 'C:\Users\rps207\Documents\Seasonal Study Data\LicorData\';%identify file path
cd(path); %change current directory to folder with files want to load
fList = ls('*.txt');  % whilst in files folder, generates list of LogData files in directory , ls('*LogData*.lvm');     
iRow=[];
% for i=1:size(fList,1);    iRow = [iRow ; i];  end;
%   LINE TO EXCLUDE PROFILES 
% for i=1:size(fList,1);  
%     if strfind(fList(i,:),'2016-07-20.txt') > 0; 
%     elseif strfind(fList(i,:),'2016-07-27.txt') > 0;
%     else iRow = [iRow ; i]; end;    end;
% fList = fList(iRow,:);

cd(mfileDir);   %change currect directory original matlab folder

% fList = fList(4,:);

randomcolor=distinguishable_colors(length(fList));

pFigs = 0; %1 = display, 0 = no display display figures or not
pequT = 1; %1 = plot equ temp, 0 = no plot equilibrator temperature and temp probes
ptxtnum=1; %1 plot numbers of bins on figure 2 , 0 = dont plot bins
axT=0; %1 = Plot time since start of the profile, 0= plot time hhmmss
randomcolor=distinguishable_colors(length(fList));

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
counter=0;counterF=0;counterG=0;dt_since_lowtide_und_return=[];dt_since_lowtide_und_outbound=[]; counterH=0;
NSOP_overlap=[];  shower_overlap=[];und_dist_L4=[];co2indall=[];


%work through 16 nsop profiles
for j = 1:length(fList(:,1));               % loops through each file in logdata list
%    close all;
    %%%%%%%%%%%%%%%%%%%   pt100s data inport %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %importpt100 text docs into workspace be selecting folder location
    path = 'C:\Users\rps207\Documents\Seasonal Study Data\PT100Data\';
    fullpathPT100a = [path fList(j,1:10) 'a.txt'];
    fullpathPT100b =[path fList(j,1:10) 'b.txt'];
    [Temperature_probe_a , Temperature_probe_b , Datetime_probe_a, Datetime_probe_b]= getPT100seasonalstudy(fullpathPT100a,fullpathPT100b);

    %%%%%%%%%%%%%%%%%%%   Licor data import %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %import Licor text docs into workspace be selecting folder location
    path = 'C:\Users\rps207\Documents\Seasonal Study Data\LicorData\';%identify file path
    fullpathLicor =[path fList(j,1:10) '.txt'];
    [LicorH20B,LicorPressure,LicorTemperature,Licor_Datetime,Licor_Datetime_doy,LicorxCO2B ] = getLicorSeasonalStudy(fullpathLicor,j,fList);
    
    %%%%%%%%%%%%%%%%%%%   CO2 calibration import %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %import Licor text docs into workspace be selecting folder location
    path = 'C:\Users\rps207\Documents\Seasonal Study Data\Licor CO2 calibrations\';%identify file path
    fullpathLicorcal =[path fList(j,1:10) '_co2cal.txt'];
    [Std_nom_values, Std_runtimes] = getCO2calSeasonalStudy(fullpathLicorcal,j,fList);
    
    %%%%%%%%%%%%%%%%%%%   microCTD data import %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %import microctd text doc into workspace be selecting folder location
    path = 'C:\Users\rps207\Documents\Seasonal Study Data\CTD\';
    fullpathCTD = [path fList(j,1:10) 'CTD.txt'];
    [CTDtemp, CTDdepthm,CTDsalinity ,CTDDT,CTDdepth_interp,CTDTemperature_interp,CTDsalinity_interp,CTDchl,CTDchl_interp] = getCTDseasonalstudy(fullpathCTD,Licor_Datetime,j,fList);  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%   flow rate data import, lag time calc, response time   %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %import flow rate sensor info
    path = 'C:\Users\rps207\Documents\Seasonal Study Data\Flow Sensor\';
    fullpathflow = [path fList(j,1:10) 'flowatest.txt'];
    [flowDT, MeasuredFlow_interp, lagsecs,respsecs,smoothflow] = getflowSeasonalStudy(fullpathflow, Licor_Datetime,j,fList);
     
    %import pump on/off time
    path = 'C:\Users\rps207\Documents\Seasonal Study Data\Pump on-off\';
    fullpathpump = [path fList(j,1:10) 'pump.txt'];
    [pumponind, pumpoffind] = getpumpSeasonalStudy(fullpathpump,Licor_Datetime);
    
    
    %import connected udnerway on/off times
    path = 'C:\Users\rps207\Documents\Seasonal Study Data\Connected underway\';
    fullpathund = [path fList(j,1:10) '_und.txt'];
   [undind1, undind2,undind3, undind4] = getunderwaySeasonalStudy(fullpathund,Licor_Datetime,fList,j);
    
    %     %import TADIC time
    path = 'C:\Users\rps207\Documents\SSBjulyCruiseData\Bottlefill data\';
    fullpathTADIC = [path fList(j,1:10) 'TADIC.txt'];
    if exist(fullpathTADIC)
    [Bottno,Fillstr, Fillend,TA,DIC,samplefco2,sampledepth,meanfilltime,semsampledepth] = getTADICjuly(fullpathTADIC,lagsecs,CTDdepth_interp,CTDTemperature_interp,CTDsalinity_interp,Licor_Datetime,j,fList);
    else
    sprintf('NO TA/DIC file')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%   identify depth bins %%%%%%%%%%%%%%%%%%%%%
    sF = 1; %sampling frequency in Hz
    spec = hspec_ss(CTDdepth_interp,CTDdepth_interp,sF);
    %figure(99); loglog(spec(:,1),spec(:,1).*spec(:,2),'.'); %plot the frequency response and use it to set cutoff frequency
    N=2; Wn = 0.008/(sF/2);  % Normalized cutoff frequency
    [B,A] = butter(N,Wn);   %butterworth filter
    %h=fvtool(B,A);          %plot magnitude response
    smoothCTD = filtfilt(B,A,CTDdepth_interp); %smooth data using filtfilt
    diffCTD=diff(smoothCTD); %find differential
    threshold=0.008;
    [pk, ind1] = findpeaks(diffCTD,'MINPEAKHEIGHT',threshold); %set threshhold and find peaks
    [pk2, ind2] = findpeaks(-diffCTD,'MINPEAKHEIGHT',threshold);%set threshhold and find peaks
    ind=([ind1;ind2]); %concatinate index matrices
    ind=([ind;pumponind;pumpoffind]);%add pump on and in as indexes
    ind=sort(ind);%sort rows so they are ordered by time rather than up and down.
    ind(ind>pumpoffind)=[];%remove inds greater then when pump went off
    ind(ind<pumponind)=[];%remove inds greater then when pump went on
    
    % work out average lag for each bin
    avglag=[];
    for w=1:length(ind)-1 ;%loop through the indexes to determine the average depths between indexes
        avglag(w)=mean(lagsecs(ind(w):ind(w+1)));  
    end
    avglag(end+1)=avglag(end); %make lag time longer so it doesn't cause problems in loop below, replaced by pump off anyway
    
    
    %%%%%%%%%%%%%%%%%%%%%%% remove false midpoints, account for transition,response and lagtime   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bincutoff=60; % set cutoffs to remove false transitions
    diffind=diff(ind);
    transsec=15; %remove the transition between depths, index points midway through transition, geustimate, transition both ways
    ind(1) = pumponind;
    indr=[];inds=[];  indt=[]; 
    for m = 1:length(diffind);
        if (diffind(m)>bincutoff);
            startInd = round(ind(m)+transsec+avglag(m)+(respsecs*2));
            stopInd = round(ind(m+1)-transsec+avglag(m+1));
            indr = [indr ; ind(m)]; % used only for plotting visualisation of unused midpoints
            inds = [inds ; startInd];
            indt = [indt ; stopInd];
        end
    end
    indt(end)=pumpoffind;%stop counting CO2 as soon as pump is turned off
    avglag=avglag(1:end-1); %undo lagtime increase extension above
    avglag=round(avglag(find(diffind>bincutoff)));
 
    
    indexrt=[];
    for k=1:length(inds);
        indexrt(k)=(abs(indt(k)-inds(k))<60);
    end

   %exception needed for temperature drop out on the 19th July
   if  fList(j,1:10)==('2016-06-30');  
   %find the indexes for where temperature cuts out
   temp_cutoutsrt= datenum('2016-06-30 13:02:00','yyyy-mm-dd HH:MM:SS');
   temp_cutoutend= datenum('2016-06-30 13:23:00','yyyy-mm-dd HH:MM:SS'); 
   tmp = abs(Licor_Datetime-temp_cutoutsrt);
   [icf tempoutind] = min(tmp); %index of closest value
   tmp2 = abs(Licor_Datetime-temp_cutoutend);
   [icf tempbackind] = min(tmp2); %index of closest value
   data=(tempoutind-indt); [value,indbeg] = min(data(data>=0)); %find the last bin end index before temp cuts out
   data=(tempbackind-inds);ipos = find(data<=0); [value,imin]=max(data(ipos)); imin=ipos(imin) ;%find first bin start index after temperature ccomes back
   inds=[inds(1:indbeg);inds(imin:end)];
   indt=[indt(1:indbeg);indt(imin:end)];
   avglag=[avglag(1:indbeg),avglag(imin:end)];
   end
    
   
   %exception needed for temperature drop out on the 19th July
   if  fList(j,1:10)==('2016-08-24');  
   %find the indexes for where temperature cuts out
   temp_cutoutsrt= datenum('2016-08-24 09:50:30','yyyy-mm-dd HH:MM:SS');
   temp_cutoutend= datenum('2016-08-24 09:55:00','yyyy-mm-dd HH:MM:SS'); 
   tmp = abs(Licor_Datetime-temp_cutoutsrt);
   [icf tempoutind] = min(tmp); %index of closest value
   tmp2 = abs(Licor_Datetime-temp_cutoutend);
   [icf tempbackind] = min(tmp2); %index of closest value
   data=(tempoutind-indt); [value,indbeg] = min(data(data>=0)); %find the last bin end index before temp cuts out
   data=(tempbackind-inds);ipos = find(data<=0); [value,imin]=max(data(ipos)); imin=ipos(imin) ;%find first bin start index after temperature ccomes back
   inds=[inds(1:indbeg);inds(imin:end)];
   indt=[indt(1:indbeg);indt(imin:end)];
   avglag=[avglag(1:indbeg),avglag(imin:end)];
   end
   
   
   if  fList(j,1:10)==('2016-07-13');  
   %find the indexes for where pump is turned off
   temp_cutoutsrt= datenum('2016-07-13 11:31:00','yyyy-mm-dd HH:MM:SS');
   temp_cutoutend= datenum('2016-07-13 11:37:00','yyyy-mm-dd HH:MM:SS'); 
   tmp = abs(Licor_Datetime-temp_cutoutsrt);
   [icf tempoutind] = min(tmp); %index of closest value
   tmp2 = abs(Licor_Datetime-temp_cutoutend);
   [icf tempbackind] = min(tmp2); %index of closest value
   data=(tempoutind-indt); [value,indbeg] = min(data(data>=0)); %find the last bin end index before temp cuts out
   data=(tempbackind-inds);ipos = find(data<=0); [value,imin]=max(data(ipos)); imin=ipos(imin) ;%find first bin start index after temperature ccomes back
   inds=[inds(1:indbeg);inds(imin:end)];
   indt=[indt(1:indbeg);indt(imin:end)];
   avglag=[avglag(1:indbeg),avglag(imin:end)]; 
   end
   
   if  fList(j,1:10)==('2016-08-10');  
   %find the indexes for where pump is turned off
   temp_cutoutsrt= datenum('2016-08-10 12:45:30','yyyy-mm-dd HH:MM:SS');
   temp_cutoutend= datenum('2016-08-10 12:46:30','yyyy-mm-dd HH:MM:SS'); 
   tmp = abs(Licor_Datetime-temp_cutoutsrt);
   [icf tempoutind] = min(tmp); %index of closest value
   tmp2 = abs(Licor_Datetime-temp_cutoutend);
   [icf tempbackind] = min(tmp2); %index of closest value
   data=(tempoutind-indt); [value,indbeg] = min(data(data>=0)); %find the last bin end index before temp cuts out
   data=(tempbackind-inds);ipos = find(data<=0); [value,imin]=max(data(ipos)); imin=ipos(imin) ;%find first bin start index after temperature ccomes back
   inds=[inds(1:indbeg);inds(imin:end)];
   indt=[indt(1:indbeg);indt(imin:end)];
   avglag=[avglag(1:indbeg),avglag(imin:end)]; 
   end
  
   
   %weird spike -remove
   if  fList(j,1:10)==('2016-09-21');  
   %find the indexes for where temperature cuts out
   temp_cutoutsrt= datenum('2016-09-21 10:24:00','yyyy-mm-dd HH:MM:SS');
   temp_cutoutend= datenum('2016-09-21 10:30:00','yyyy-mm-dd HH:MM:SS'); 
   tmp = abs(Licor_Datetime-temp_cutoutsrt);
   [icf tempoutind] = min(tmp); %index of closest value
   tmp2 = abs(Licor_Datetime-temp_cutoutend);
   [icf tempbackind] = min(tmp2); %index of closest value
   data=(tempoutind-indt); [value,indbeg] = min(data(data>=0)); %find the last bin end index before temp cuts out
   data=(tempbackind-inds);ipos = find(data<=0); [value,imin]=max(data(ipos)); imin=ipos(imin) ;%find first bin start index after temperature ccomes back
   inds=[inds(1:indbeg);inds(imin:end)];
   indt=[indt(1:indbeg);indt(imin:end)];
   avglag=[avglag(1:indbeg),avglag(imin:end)];
   end
   
   
      %water in syysteme -remove
   if  fList(j,1:10)==('2016-06-15');  
   %find the indexes for where temperature cuts out
   temp_cutoutsrt= datenum('2016-06-15 10:00:00','yyyy-mm-dd HH:MM:SS');
   temp_cutoutend= datenum('2016-06-15 11:30:00','yyyy-mm-dd HH:MM:SS'); 
   tmp = abs(Licor_Datetime-temp_cutoutsrt);
   [icf tempoutind] = min(tmp); %index of closest value
   tmp2 = abs(Licor_Datetime-temp_cutoutend);
   [icf tempbackind] = min(tmp2); %index of closest value
   data=(tempoutind-indt); [value,indbeg] = min(data(data>=0)); %find the last bin end index before temp cuts out
   data=(tempbackind-inds);ipos = find(data<=0); [value,imin]=max(data(ipos)); imin=ipos(imin) ;%find first bin start index after temperature ccomes back
   inds=[inds(1:indbeg);inds(imin:end)];
   indt=[indt(1:indbeg);indt(imin:end)];
   avglag=[avglag(1:indbeg),avglag(imin:end)];
   end
   
   
   
   
   
        %weird spike -remove
   if  fList(j,1:10)==('2016-08-10');  
   %find the indexes for where temperature cuts out
   temp_cutoutsrt= datenum('2016-08-10 14:19:00','yyyy-mm-dd HH:MM:SS');
   temp_cutoutend= datenum('2016-08-10 14:22:00','yyyy-mm-dd HH:MM:SS'); 
   tmp = abs(Licor_Datetime-temp_cutoutsrt);
   [icf tempoutind] = min(tmp); %index of closest value
   tmp2 = abs(Licor_Datetime-temp_cutoutend);
   [icf tempbackind] = min(tmp2); %index of closest value
   replacementvec=linspace(LicorxCO2B(tempoutind),LicorxCO2B(tempbackind),tempbackind-tempoutind+1);
   LicorxCO2B(tempoutind:tempbackind)=replacementvec;
   end
   
   
           %weird spike -remove
   if  fList(j,1:10)==('2016-06-30');  
   %find the indexes for where temperature cuts out
   temp_cutoutsrt= datenum('2016-06-30 10:34:00','yyyy-mm-dd HH:MM:SS');
   temp_cutoutend= datenum('2016-06-30 10:36:30','yyyy-mm-dd HH:MM:SS'); 
   tmp = abs(Licor_Datetime-temp_cutoutsrt);
   [icf tempoutind] = min(tmp); %index of closest value
   tmp2 = abs(Licor_Datetime-temp_cutoutend);
   [icf tempbackind] = min(tmp2); %index of closest value
   replacementvec=linspace(LicorxCO2B(tempoutind),LicorxCO2B(tempbackind),tempbackind-tempoutind+1);
   LicorxCO2B(tempoutind:tempbackind)=replacementvec;
   end
   

   
   %remove weird negative spike
   if  fList(j,1:10)==('2016-06-30');  
   %find the indexes for where temperature cuts out
   temp_cutoutsrt= datenum('2016-06-30 10:34:00','yyyy-mm-dd HH:MM:SS');
   temp_cutoutend= datenum('2016-06-30 11:36:30','yyyy-mm-dd HH:MM:SS'); 
   tmp = abs(Licor_Datetime-temp_cutoutsrt);
   [icf tempoutind] = min(tmp); %index of closest value
   tmp2 = abs(Licor_Datetime-temp_cutoutend);
   [icf tempbackind] = min(tmp2); %index of closest value
   data=(tempoutind-indt); [value,indbeg] = min(data(data>=0)); %find the last bin end index before temp cuts out
   data=(tempbackind-inds);ipos = find(data<=0); [value,imin]=max(data(ipos)); imin=ipos(imin) ;%find first bin start index after temperature ccomes back
   inds=[inds(1:indbeg);inds(imin:end)];
   indt=[indt(1:indbeg);indt(imin:end)];
   avglag=[avglag(1:indbeg),avglag(imin:end)];
   end
   
   
   
   
   
   %correct indexes to remove point when pump was turned off
   if  fList(j,1:10)==('2016-08-10');  
       avglag(4)=200;
       inds(9)=inds(9)+20;
       indt(8)=[];inds(8)=[];avglag(8)=[];
       indt(4)=[];inds(4)=[];avglag(4)=[];
       indt(2)=[];inds(2)=[];avglag(2)=[];
   end
   
      %correct indexes to remove point when pump was turned off
   if  fList(j,1:10)==('2016-08-04'); 
       indt(24)=[];inds(24)=[];avglag(24)=[];
       indt(23)=[];inds(23)=[];avglag(23)=[];
       indt(20)=[];inds(20)=[];avglag(20)=[];
       indt(7)=[];inds(7)=[];avglag(7)=[];
       indt(6)=[];inds(6)=[];avglag(6)=[];
       indt(5)=[];inds(5)=[];avglag(5)=[];
   end
   
   %correct indexes to remove point when pump was turned off
   if  fList(j,1:10)==('2016-05-11'); 
       indt(20)=[];inds(20)=[];avglag(20)=[];
   end

%       %correct indexes to remove point when pump was turned off
%    if  fList(j,1:10)==('2016-06-30'); 
%        indt(39)=[];inds(39)=[];avglag(39)=[];
%    end
   
    %remove wrongly identified false transitions
   if  fList(j,1:10)==('2016-06-22'); 
       indt(8)=[];inds(8)=[];avglag(8)=[];
   end
   if  fList(j,1:10)==('2016-07-07'); 
       indt(27)=[];inds(27)=[];avglag(27)=[];
   end
     if  fList(j,1:10)==('2016-09-21'); 
       indt(31)=[];inds(31)=[];avglag(31)=[];
   end
   if  fList(j,1:10)==('2016-07-13');
       indt(28)=[];inds(28)=[];avglag(28)=[];
       indt(12)=[];inds(12)=[];avglag(12)=[];
       indt(6)=[];inds(6)=[];avglag(6)=[];
       indt(5)=[];inds(5)=[];avglag(5)=[];
   end
%    if  fList(j,1:10)==('2016-08-24'); 
%        indt(25)=[];inds(25)=[];avglag(25)=[];
%    end
   if  fList(j,1:10)==('2016-09-02'); 
       indt(24)=[];inds(24)=[];avglag(24)=[];
       indt(14)=[];inds(14)=[];avglag(14)=[];
   end
   if  fList(j,1:10)==('2016-09-15'); 
       indt(9)=[];inds(9)=[];avglag(9)=[];
       indt(6)=[];inds(6)=[];avglag(6)=[];
   end
   if  fList(j,1:10)==('2016-08-24'); 
%        indt(18)=[];inds(18)=[];avglag(18)=[];
   end
   if  fList(j,1:10)==('2016-09-21'); 
       indt(34)=[];inds(34)=[];avglag(34)=[];
   end
   
   
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%     CO2 calculations  Part 1  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    Temperature_probe_a_interp= interp1(Datetime_probe_a,Temperature_probe_a ,Datetime_probe_b); %interpolate both temperature probes so they are the same lemngth 
    Equilibrator_Temperature = (Temperature_probe_a_interp + Temperature_probe_b)/2; %define equilibrator temperature as the average between inlet and outlet probes
    Equilibrator_Temperature_kelvin=Equilibrator_Temperature +273.15; %convert  into kelvin
    Equilibrator_Datetime=Datetime_probe_b; % use probe b time as equilibrator time
    Equilibrator_Temperature_kelvin_interp= interp1(Equilibrator_Datetime,Equilibrator_Temperature_kelvin ,Licor_Datetime);%interp equilibrator temperature to the length of the licor data
    LicorTemperaturekelvin=LicorTemperature +273.15; %convert to kelvin
    %use Wagner and garb 2002 to get water vapour
    temp_mod = 1-Equilibrator_Temperature_kelvin_interp./647.096;
    vapor_0sal_kPa=(22.064e3)*(exp((647.076./(Equilibrator_Temperature_kelvin_interp)).*((-7.85951783*temp_mod)+((1.84408259)*(temp_mod.^(3/2)))+(-11.7866497*(temp_mod.^3))+(22.6807411*(temp_mod.^3.5))+(-15.9618719*(temp_mod.^4))+(1.80122502*(temp_mod.^7.5)))));
    %Correct vapor pressure for salinity
    molality = 31.998 * CTDsalinity_interp ./(1e3-1.005*CTDsalinity_interp);
    osmotic_coef = 0.90799 -0.08992*(0.5*molality) +0.18458*(0.5*molality).^2 -0.07395*(0.5*molality).^3 -0.00221*(0.5*molality).^4;
    vapor_press_kPa = vapor_0sal_kPa .* exp(-0.018 * osmotic_coef .* molality);
    Vapour_pressure_mbar = 10*(vapor_press_kPa/101.32501);%Convert to mbar
    %Vapour_pressure_mbar= 0.981 * exp(14.32602 - (5306.83./(Equilibrator_Temperature_kelvin_interp))) * 1013.25;% temperature only not salinity(Cooper 1998)
    Vapour_pressure_atm=0.000986923267 .*Vapour_pressure_mbar;%convert to atm 1 millibar = 0.000986923267 atmosphere
    Atmospheric_pressure_kpa= LicorPressure; %assume that cell B pressure is equivalent to ambient pressure 
    Atmospheric_pressure_pa= Atmospheric_pressure_kpa *1000; %convert from kpa into pa
    Atmospheric_pressure_atm= Atmospheric_pressure_kpa *0.00986923267;  %convert from kpa into atm
    Equilibrator_pressure_atm= Atmospheric_pressure_kpa *0.00986923267; % assume that the pressure in the licor is representative of that in the equilibrator
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%     CO2 standards calculation  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if fList(j,1:10)==('2016-04-27');
        % get licor index point where calibration gas is venting through Licor in datenum format
        %add 7 secs for pressure to level out
        [~, precal250] = (min(abs(Licor_Datetime - (datenum(Std_runtimes(2),'yyyy-mm-dd HH:MM:SS')))));
        [~, precal380] = (min(abs(Licor_Datetime - (datenum(Std_runtimes(3),'yyyy-mm-dd HH:MM:SS')))));
        [~, precal450] = (min(abs(Licor_Datetime - (datenum(Std_runtimes(4),'yyyy-mm-dd HH:MM:SS')))));
        [~, postcal250] = (min(abs(Licor_Datetime - (datenum(Std_runtimes(6),'yyyy-mm-dd HH:MM:SS')))));
        [~, postcal380] = (min(abs(Licor_Datetime - (datenum(Std_runtimes(7),'yyyy-mm-dd HH:MM:SS')))));
        [~, postcal450] = (min(abs(Licor_Datetime - (datenum(Std_runtimes(8),'yyyy-mm-dd HH:MM:SS')))));
        precal250=precal250+7; precal380=precal380+7; precal450=precal450+7;     postcal250=postcal250+7; postcal380=postcal380+7; postcal450=postcal450+7;
        % calibrated xco2 for standards are
        true250std=244.46;  true380std= 376.27; true450std=447.44; truestd=[true250std,true380std,true450std];
        %interpolated xco2 values for length of licor file
        veclength=(1:length(LicorxCO2B))';
        interpolated250x=LicorxCO2B(precal250)+(((LicorxCO2B(postcal250)-LicorxCO2B(precal250)).*(( Licor_Datetime(veclength) -(datenum(Std_runtimes(2),'yyyy-mm-dd HH:MM:SS')))./((datenum(Std_runtimes(6),'yyyy-mm-dd HH:MM:SS'))-(datenum(Std_runtimes(2),'yyyy-mm-dd HH:MM:SS'))))));
        interpolated380x= LicorxCO2B(precal380)+(((LicorxCO2B(postcal380)-LicorxCO2B(precal380))*(( Licor_Datetime(veclength) -(datenum(Std_runtimes(3),'yyyy-mm-dd HH:MM:SS')))./((datenum(Std_runtimes(7),'yyyy-mm-dd HH:MM:SS'))-(datenum(Std_runtimes(3),'yyyy-mm-dd HH:MM:SS'))))));
        interpolated450x= LicorxCO2B(precal450)+(((LicorxCO2B(postcal450)-LicorxCO2B(precal450))*(( Licor_Datetime(veclength) -(datenum(Std_runtimes(4),'yyyy-mm-dd HH:MM:SS')))./((datenum(Std_runtimes(8),'yyyy-mm-dd HH:MM:SS'))-(datenum(Std_runtimes(4),'yyyy-mm-dd HH:MM:SS'))))));
        interpolatedvec=[interpolated250x,interpolated380x,interpolated450x];
        %calculate the linear fit for each data point
        A=[truestd(1:size(interpolatedvec, 2));ones(1,size(interpolatedvec, 2))]'; AA=kron(speye(length(LicorxCO2B)), A); bb=interpolatedvec'; bb=bb(:);
        z=AA\bb; z=reshape(z,2,[]); slopes=z(1,:)'; intercept=z(2,:)';    clearvars AA bb z
       %check calibration numbers are reasonable
       LicorxCO2B(precal250);
       LicorxCO2B(postcal250);
       LicorxCO2B(precal380);
       LicorxCO2B(postcal380);
       LicorxCO2B(precal450);
       LicorxCO2B(postcal450);
    else
      %LICOR WASNT LOGGING DURING PRE CAL- thankfully wrote down the
   if  fList(j,1:10)==('2016-09-21'); 
        % get licor index point where calibration gas is venting through Licor in datenum format
        %6 toc hange valve and for pressure to leevl out, add vent time in
        %seconds as indexes are seconds 
        [~, postcal250] = (min(abs(Licor_Datetime - (datenum(Std_runtimes(3),'yyyy-mm-dd HH:MM:SS')))));
        [~, postcal450] = (min(abs(Licor_Datetime - (datenum(Std_runtimes(4),'yyyy-mm-dd HH:MM:SS')))));
         postcal250=postcal250+6; postcal450=postcal450+6;
        % calibrated xco2 for standards are
        true250std=263.04; true450std=483.36; truestd=[true250std,true450std];
        %interpolated xco2 values for length of licor file
        veclength=(1:length(LicorxCO2B))';
        interpolated250x=264.35+(((LicorxCO2B(postcal250)-264.35).*(( Licor_Datetime(veclength) -(datenum(Std_runtimes(1),'yyyy-mm-dd HH:MM:SS')))./((datenum(Std_runtimes(3),'yyyy-mm-dd HH:MM:SS'))-(datenum(Std_runtimes(1),'yyyy-mm-dd HH:MM:SS'))))));
        interpolated450x=469.80+(((LicorxCO2B(postcal450)-469.80).*(( Licor_Datetime(veclength) -(datenum(Std_runtimes(2),'yyyy-mm-dd HH:MM:SS')))./((datenum(Std_runtimes(4),'yyyy-mm-dd HH:MM:SS'))-(datenum(Std_runtimes(2),'yyyy-mm-dd HH:MM:SS'))))));
        interpolatedvec=[interpolated250x,interpolated450x];
        %calculate the linear fit for each data point
        A=[truestd(1:size(interpolatedvec, 2));ones(1,size(interpolatedvec, 2))]'; AA=kron(speye(length(LicorxCO2B)), A); bb=interpolatedvec'; bb=bb(:);
        z=AA\bb; z=reshape(z,2,[]); slopes=z(1,:)'; intercept=z(2,:)';    clearvars AA bb z
   else
         % get licor index point where calibration gas is venting through Licor in datenum format
        %6 toc hange valve and for pressure to leevl out, add vent time in
        %seconds as indexes are seconds 
        [~, precal250] = (min(abs(Licor_Datetime - (datenum(Std_runtimes(1),'yyyy-mm-dd HH:MM:SS')))));
        [~, precal450] = (min(abs(Licor_Datetime - (datenum(Std_runtimes(2),'yyyy-mm-dd HH:MM:SS')))));
        [~, postcal250] = (min(abs(Licor_Datetime - (datenum(Std_runtimes(3),'yyyy-mm-dd HH:MM:SS')))));
        [~, postcal450] = (min(abs(Licor_Datetime - (datenum(Std_runtimes(4),'yyyy-mm-dd HH:MM:SS')))));
        precal250=precal250+6;  precal450=precal450+6; postcal250=postcal250+6; postcal450=postcal450+6;
        % calibrated xco2 for standards are
        true250std=263.04; true450std=483.36; truestd=[true250std,true450std];
        %interpolated xco2 values for length of licor file
        veclength=(1:length(LicorxCO2B))';
        interpolated250x=LicorxCO2B(precal250)+(((LicorxCO2B(postcal250)-LicorxCO2B(precal250)).*(( Licor_Datetime(veclength) -(datenum(Std_runtimes(1),'yyyy-mm-dd HH:MM:SS')))./((datenum(Std_runtimes(3),'yyyy-mm-dd HH:MM:SS'))-(datenum(Std_runtimes(1),'yyyy-mm-dd HH:MM:SS'))))));
        interpolated450x=LicorxCO2B(precal450)+(((LicorxCO2B(postcal450)-LicorxCO2B(precal450)).*(( Licor_Datetime(veclength) -(datenum(Std_runtimes(2),'yyyy-mm-dd HH:MM:SS')))./((datenum(Std_runtimes(4),'yyyy-mm-dd HH:MM:SS'))-(datenum(Std_runtimes(2),'yyyy-mm-dd HH:MM:SS'))))));
        interpolatedvec=[interpolated250x,interpolated450x];
        %calculate the linear fit for each data point
        A=[truestd(1:size(interpolatedvec, 2));ones(1,size(interpolatedvec, 2))]'; AA=kron(speye(length(LicorxCO2B)), A); bb=interpolatedvec'; bb=bb(:);
        z=AA\bb; z=reshape(z,2,[]); slopes=z(1,:)'; intercept=z(2,:)';    clearvars AA bb z
        %check calibration numbers are reasonable
        LicorxCO2B(precal250);
        LicorxCO2B(postcal250);
        LicorxCO2B(precal450);
        LicorxCO2B(postcal450);
    end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%     CO2 calculations  Part 2  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    LicorxCO2Bcal=(LicorxCO2B.*slopes)+intercept;
    
%         efficency correction
    if  fList(j,1:10)==('2016-04-27');
           LicorxCO2Bcal=LicorxCO2Bcal*(1/0.9638515);   
    elseif fList(j,1:10)==('2016-05-11');
          LicorxCO2Bcal=LicorxCO2Bcal*(1/0.9638515);   
    elseif fList(j,1:10)==('2016-05-26');
          LicorxCO2Bcal=LicorxCO2Bcal*(1/0.9638515);   
    elseif fList(j,1:10)==('2016-06-10');
           LicorxCO2Bcal=LicorxCO2Bcal*(1/0.9638515);   
    elseif fList(j,1:10)==('2016-06-15');
           LicorxCO2Bcal=LicorxCO2Bcal*(1/0.9712889);   
    elseif fList(j,1:10)==('2016-06-22');
          LicorxCO2Bcal=LicorxCO2Bcal*(1/0.976612 );   
    elseif fList(j,1:10)==('2016-06-30');
          LicorxCO2Bcal=LicorxCO2Bcal*(1/0.979082);   
    elseif fList(j,1:10)==('2016-07-07');
           LicorxCO2Bcal=LicorxCO2Bcal*(1/0.9629);  
    elseif fList(j,1:10)==('2016-07-13');
            LicorxCO2Bcal=LicorxCO2Bcal*(1/0.88209); 
    elseif fList(j,1:10)==('2016-07-20');
        
    elseif fList(j,1:10)==('2016-07-27');
        
    elseif fList(j,1:10)==('2016-08-04');
          LicorxCO2Bcal=LicorxCO2Bcal*(1/0.9241);
   elseif fList(j,1:10)==('2016-08-10'); %SAME AS 08/04
          LicorxCO2Bcal=LicorxCO2Bcal*(1/0.9241);      
   elseif fList(j,1:10)==('2016-08-17');
           LicorxCO2Bcal=LicorxCO2Bcal*(1/0.89808);   
   elseif fList(j,1:10)==('2016-08-24');
          LicorxCO2Bcal=LicorxCO2Bcal*(1/0.950349);  
   elseif fList(j,1:10)==('2016-09-02');
          LicorxCO2Bcal=LicorxCO2Bcal*(1/0.954459);
   elseif fList(j,1:10)==('2016-09-15');
         LicorxCO2Bcal=LicorxCO2Bcal*(1/0.781612); 
   elseif fList(j,1:10)==('2016-09-15');
        LicorxCO2Bcal=LicorxCO2Bcal*(1/ 0.8425);     
    else
         LicorxCO2Bcal=LicorxCO2Bcal*(1/0.843506);       
    end

    LicorpCO2= LicorxCO2Bcal.*(Equilibrator_pressure_atm - Vapour_pressure_atm);%Correction for water vapour pressure as described in (dickson(2007)-section 8.5.3 TO GET PCO2
    Temperature_probe_a_interp= interp1(Datetime_probe_a,Temperature_probe_a ,Datetime_probe_b); %interpolate both temperature probes so they are the same lemngth 
    Equilibrator_Temperature = (Temperature_probe_a_interp + Temperature_probe_b)/2; %define equilibrator temperature as the average between inlet and outlet probes
    Equilibrator_Temperature_kelvin=Equilibrator_Temperature +273.15; %convert temperature into kelvin
    Equilibrator_Datetime=Datetime_probe_b; % use probe b time as equilibrator time
    Equilibrator_Temperature_kelvin_interp= interp1(Equilibrator_Datetime,Equilibrator_Temperature_kelvin ,Licor_Datetime);%interp equilibrator temperature to the length of the licor data
    BCO2te= -1636.75 + (12.0408*(Equilibrator_Temperature_kelvin_interp)) - (3.27957*0.01*(Equilibrator_Temperature_kelvin_interp).^2) + (3.16528*0.00001*(Equilibrator_Temperature_kelvin_interp).^3);%units of cm^3 mol^-1 ,determine the two viral coefficents of co2, (Dickson 2007) SOP 24 and p98 section 8.3 use equilibrator_pressure and equilibrator_temperature
    deltaCO2te= 57.7 - 0.118*( Equilibrator_Temperature_kelvin_interp); %units of cm^3 mol^-1
    R=8.31447; %specific gas constant J/Mol*K
    Licor_fco2= LicorpCO2.*exp(((BCO2te+(2*deltaCO2te))*0.000001.*(Atmospheric_pressure_pa))./(R *Equilibrator_Temperature_kelvin_interp));%use viral coefficents to determine fco2
    Equilibrator_Temperature_interp= interp1(Equilibrator_Datetime, Equilibrator_Temperature,Licor_Datetime); %interpolate equilibrator to licor date time
    
% old code that is used for doing CO2 analysis point by point, keep in for
% purpose of plotting below
    Y=Licor_Datetime; %  %Create 2 random matrices with the same value for n
    g=(Licor_Datetime-(lagsecs/(24*60*60))); %Create 2 random matrices with the same value for n
    IDX = knnsearch(Y(:),g(:)) ;%For the number of values of length of A  find the nearest in the first column of B , designate this as an index either 1 or 2
    CTDTemperature_interp_lag=CTDTemperature_interp(IDX);%use this index to obtain ancillary information from matrix B
    Licor_fco2_surface=Licor_fco2.*exp(0.0423*(CTDTemperature_interp_lag - Equilibrator_Temperature_interp));%correction of co2 to sea surface temperature , (Dickson 2007) p98 8.4
   
    
    load (['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\Underwayinterp\' fList(j,1:10) '.mat'],'AtmTemp_interp','new_names2','AtmTemp_interplongtime','underwayDOY_interp','Atten_interp','underwayDOY_interplongtime','Atten_interplongtime','underwayDfCO2_interp','Cdom_interp','underwayDfCO2_interplongtime','Cdom_interplongtime','underwayDpCO2_interp','Chla_interp','underwayDpCO2_interplongtime','Chla_interplongtime','underwayLat_interp','underwayLat_interplongtime','Humidity_interp','underwayLong_interp','Humidity_interplongtime','underwayLong_interplongtime','Latitude_interp','underwayO2_optode_interp','Latitude_interplongtime','underwayO2_optode_interplongtime','Longitude_interp','underwayP_atm_interp','Longitude_interplongtime','underwayP_atm_interplongtime','O2saturation_interp','underwayRel_Wind_interp','O2saturation_interplongtime','underwayRel_Wind_interplongtime','Parport_interp','underwayS_sea_interp','Parport_interplongtime','underwayS_sea_interplongtime','Pressure_interp','underwayT_sea_interp','Pressure_interplongtime','underwayT_sea_interplongtime','Salinity_interp','underwayWind_Dir_interp','Salinity_interplongtime','underwayWind_Dir_interplongtime','SeaTemp_interp','underway_DT','SeaTemp_interplongtime','underway_DT_interp','Time_interp','underway_DT_interplongtime','Turbidity_interp','underwaychlorophyll_interp','Turbidity_interplongtime','underwaychlorophyll_interplongtime','Und_DT_interp','underwayfCO2_atm_interp','Und_DT_interplongtime','underwayfCO2_atm_interplongtime','Winddirabs_interp','underwayfCO2_sw_interp','Winddirabs_interplongtime','underwayfCO2_sw_interplongtime','Windspeedabs_interp','underwaypCO2_atm_interp','Windspeedabs_interplongtime','underwaypCO2_atm_interplongtime','underwaypCO2_sw_interp','underwaypCO2_sw_interplongtime','names2')
    
    
    
    %NOTE for underway sections when not connected to underway use the
    %correct salinity and temperature from the underway not the CTD
    
    %note need to apply 79 second lag to underway measurements, simpler
    %than the avglag stuff i did for the NSOP code. Firstly use the
    %adjusted sea temp to get fco2 then apply the lag back to get the co2
    %on the correct timeframe
    %lag part 1
    %create matrix for seatemp_interp that can be shifted.
    underwaylag=79;
    a=nan(underwaylag,1);
    SeaTemp_interp_shifted=[a ; SeaTemp_interp(1:end-underwaylag)];
    Salinity_interp_shifted=[a ; Salinity_interp(1:end-underwaylag)];
    %now temperature used is correct
    
    molalityund = 31.998 * SeaTemp_interp_shifted ./(1e3-1.005*Salinity_interp_shifted);
    osmotic_coefund = 0.90799 -0.08992*(0.5*molalityund) +0.18458*(0.5*molalityund).^2 -0.07395*(0.5*molalityund).^3 -0.00221*(0.5*molalityund).^4;
    vapor_press_kPaund = vapor_0sal_kPa .* exp(-0.018 * osmotic_coefund .* molalityund);
    Vapour_pressure_mbarund = 10*(vapor_press_kPaund/101.32501);%Convert to mbar
    %Vapour_pressure_mbar= 0.981 * exp(14.32602 - (5306.83./(Equilibrator_Temperature_kelvin_interp))) * 1013.25;% temperature only not salinity(Cooper 1998)
    Vapour_pressure_atmund=0.000986923267 .*Vapour_pressure_mbarund;%convert to atm 1 millibar = 0.000986923267 atmosphere
    LicorpCO2und= LicorxCO2Bcal.*(Equilibrator_pressure_atm - Vapour_pressure_atmund);%Correction for water vapour pressure as described in (dickson(2007)-section 8.5.3 TO GET PCO2
    Licor_fco2und= LicorpCO2und.*exp(((BCO2te+(2*deltaCO2te))*0.000001.*(Atmospheric_pressure_pa))./(R *Equilibrator_Temperature_kelvin_interp));%use viral coefficents to determine fco2
    Licor_fco2_surfaceund=Licor_fco2und.*exp(0.0423*(SeaTemp_interp_shifted - Equilibrator_Temperature_interp));
    

%%%%%%% Use this code to move Licor Timeseries backwards for plotting
    Zp=(Licor_Datetime-(lagsecs/(24*60*60))); %Create 2 random matrices with the same value for n
    IDX2 = knnsearch(Y(:),Zp(:)) ;%For the number of values of length of A  find the nearest in the first column of B , designate this as an index either 1 or 2
    Licor_fco2_lag=Licor_fco2(IDX2);
    LicorxCO2B_lag=LicorxCO2B(IDX2);
    LicorpCO2_lag=LicorpCO2(IDX2);
    LicorxCO2Bcal_lag=LicorxCO2Bcal(IDX2);
    Licor_fco2_surface_lag=Licor_fco2_surface(IDX2);
    
   Temperature_probe_a_interpx= interp1(Datetime_probe_a,Temperature_probe_a ,Licor_Datetime);%interp equilibrator temperature to the length of the licor data
   Temperature_probe_a_interp_lag=Temperature_probe_a_interpx(IDX2);
   Temperature_probe_b_interpx= interp1(Datetime_probe_b,Temperature_probe_b ,Licor_Datetime);%interp equilibrator temperature to the length of the licor data
   Temperature_probe_b_interp_lag=Temperature_probe_b_interpx(IDX2);
   Equilibrator_Temperature_interp_lag=Equilibrator_Temperature_interp(IDX2);
    
%for the average binned co2 first bin ctd and equ temperature and fco2
    for z=1:length(inds); %loop through the indexes to determine the average depths between indexes
        avgEquilibrator_Temperature_interp(z)=mean(Equilibrator_Temperature_interp((inds(z)):(indt(z))));
        stdEquilibrator_Temperature_interp(z)=std(Equilibrator_Temperature_interp(inds(z):indt(z))); 
        semEquilibrator_Temperature_interp(z)=sem(Equilibrator_Temperature_interp(inds(z):indt(z))); 
        avgfco2(z)=mean(Licor_fco2(inds(z):indt(z)));%loop through the indexes to determine the average temperature between indexes
        avgCTDTemperature_interp_lag(z)=mean(CTDTemperature_interp((inds(z)-avglag(z)):(indt(z)-avglag(z))));
        avgLicorxCO2B(z)=mean(Equilibrator_Temperature_interp(inds(z):indt(z)));  
    end

    % funny business going on here try
    avgEquilibrator_Temperature_interp=avgEquilibrator_Temperature_interp(1:length(inds));
    stdEquilibrator_Temperature_interp=stdEquilibrator_Temperature_interp(1:length(inds));
    avgCTDTemperature_interp_lag=avgCTDTemperature_interp_lag(1:length(inds));
    avgfco2=avgfco2(1:length(inds));

    
    avgfco2surf=avgfco2.*exp(0.0423.*(avgCTDTemperature_interp_lag - avgEquilibrator_Temperature_interp));%correction of co2 to sea surface temperature , (Dickson 2007) p98 8.4
    %clearvars R Equilibrator_Temperature_kelvin  deltaCO2te Vapour_pressure_atm Vapour_pressure_mbar Equilibrator_Datetime Datetime_probe_b Datetime_probe_a BCO2te  Temperature_probe_a Temperature_probe_a_interp Temperature_probe_b
                
    
    %apply lag 2 here
   
    undind4=undind4-underwaylag;
    %this is underway
    Licor1_preshift=Licor_fco2_surfaceund(1:inds(1)-1);
    Licor1=[Licor1_preshift(underwaylag+1:end);a];
    
    %this is nsop
    Licor2=Licor_fco2_surface_lag(inds(1):indt(end));
        
    %this is underway
    Licor3_preshift=Licor_fco2_surfaceund(indt(end)+1:end);
    Licor3=[Licor3_preshift(underwaylag+1:end);a]; 
    
    %combined
    Licor_fco2_combined=[Licor1;Licor2;Licor3];
        
    %exception for weird co2 spike
      if  fList(j,1:10)==('2016-08-17');  
   %find the indexes for where pump is turned off
   temp_cutoutsrt= datenum('2016-08-17 12:59:30','yyyy-mm-dd HH:MM:SS');
   temp_cutoutend= datenum('2016-08-17 13:03:00','yyyy-mm-dd HH:MM:SS'); 
   [tmp tempoutind] = find(Licor_Datetime>temp_cutoutsrt & Licor_Datetime<temp_cutoutend);
   Licor_fco2_combined(tmp)=323.60;
      end
   
    %load in underway variables interpolated onto same length as Licor
    load (['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\Underwayinterp\' fList(j,1:10) '.mat'],'AtmTemp_interp','new_names2','AtmTemp_interplongtime','underwayDOY_interp','Atten_interp','underwayDOY_interplongtime','Atten_interplongtime','underwayDfCO2_interp','Cdom_interp','underwayDfCO2_interplongtime','Cdom_interplongtime','underwayDpCO2_interp','Chla_interp','underwayDpCO2_interplongtime','Chla_interplongtime','underwayLat_interp','underwayLat_interplongtime','Humidity_interp','underwayLong_interp','Humidity_interplongtime','underwayLong_interplongtime','Latitude_interp','underwayO2_optode_interp','Latitude_interplongtime','underwayO2_optode_interplongtime','Longitude_interp','underwayP_atm_interp','Longitude_interplongtime','underwayP_atm_interplongtime','O2saturation_interp','underwayRel_Wind_interp','O2saturation_interplongtime','underwayRel_Wind_interplongtime','Parport_interp','underwayS_sea_interp','Parport_interplongtime','underwayS_sea_interplongtime','Pressure_interp','underwayT_sea_interp','Pressure_interplongtime','underwayT_sea_interplongtime','Salinity_interp','underwayWind_Dir_interp','Salinity_interplongtime','underwayWind_Dir_interplongtime','SeaTemp_interp','underway_DT','SeaTemp_interplongtime','underway_DT_interp','Time_interp','underway_DT_interplongtime','Turbidity_interp','underwaychlorophyll_interp','Turbidity_interplongtime','underwaychlorophyll_interplongtime','Und_DT_interp','underwayfCO2_atm_interp','Und_DT_interplongtime','underwayfCO2_atm_interplongtime','Winddirabs_interp','underwayfCO2_sw_interp','Winddirabs_interplongtime','underwayfCO2_sw_interplongtime','Windspeedabs_interp','underwaypCO2_atm_interp','Windspeedabs_interplongtime','underwaypCO2_atm_interplongtime','underwaypCO2_sw_interp','underwaypCO2_sw_interplongtime','names2');
    %calculate distance form starting point of profile
    Lat1=ones(length(Latitude_interp),1).*Latitude_interp(ind(1));
    Long1=ones(length(Longitude_interp),1).*Longitude_interp(ind(1)); dist=[];
    for t=1:length(Longitude_interp);
        dist(t) = pos2dist(Lat1(t),Long1(t),Latitude_interp(t),Longitude_interp(t),2);
    end
    
    dist=dist';
    
    load('C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\Plymouth_tides.mat');
    Plymouth_tidal_height_interp=interp1(Plymouth_tidal_height_dt,Plymouth_tidal_height,Licor_Datetime);

    
    
    %calculate density of profile
    avgpresdb1 = sw_pres(CTDdepth_interp,max(Latitude_interp(1)));
    CTDdensity_interp=sw_dens(CTDsalinity_interp,CTDTemperature_interp,avgpresdb1);
    
    %%%%%%%%%%%%%%%%%%%%%%% bin CTD data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    avgtime=[];     avgdepth=[]; stddepth=[];  semdepth=[];    avgsal=[]; stdsal=[];  semsal=[];
    avgtemp=[]; stdtemp=[];  semtemp=[];    avgfco2=[]; stdfco2=[];  semfco2=[];   avgflow=[]; avgtide=[];
    avgWindspeed=[];  avgden=[];stdden=[],semden=[];avgLongitude=[];avgLatitude=[]; avgchl=[];stdchl=[];semchl=[];
    stdEqu_temp=[];semEqu_temp=[]; avgfco2und=[];stdfco2und=[];semfco2und=[]; avgEqu_temp=[];   stdEqu_temp=[];semEqu_temp=[];
    
    for v=1:length(inds) %loop through the indexes to determine the average depths between indexes

        avgtime(v)=mean(Licor_Datetime((inds(v)-avglag(v)):(indt(v)-avglag(v)))); % make this CTD times as this is what CTD variables and CO2(INDS/INDT) actually is at this time
        avgflow(v)=mean(MeasuredFlow_interp(inds(v):indt(v)));

        avgdepth(v)=mean(CTDdepth_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        stddepth(v)=std(CTDdepth_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        semdepth(v)=sem(CTDdepth_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
      
        avgsal(v)=mean(CTDsalinity_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));%loop through the indexes to determine the average salinity between indexes
        stdsal(v)=std(CTDsalinity_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        semsal(v)=sem(CTDsalinity_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));

        avgtemp(v)=mean(CTDTemperature_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));%loop through the indexes to determine the average temperature between indexes
        stdtemp(v)=std(CTDTemperature_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        semtemp(v)=sem(CTDTemperature_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        
        avgden(v)=mean(CTDdensity_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));%loop through the indexes to determine the average temperature between indexes
        stdden(v)=std(CTDdensity_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        semden(v)=semjuly(CTDdensity_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        
        avgchl(v)=mean(CTDchl_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));%loop through the indexes to determine the average temperature between indexes
        stdchl(v)=std(CTDchl_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        semchl(v)=semjuly(CTDchl_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));

        avgfco2und(v)=mean(underwayfCO2_sw_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));%loop through the indexes to determine the average temperature between indexes
        stdfco2und(v)=std(underwayfCO2_sw_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        semfco2und(v)=semjuly(underwayfCO2_sw_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        
        avgfco2(v)=mean(Licor_fco2(inds(v):indt(v)));%loop through the indexes to determine the average temperature between indexes
        stdfco2(v)=std(Licor_fco2(inds(v):indt(v)));
        semfco2(v)=sem(Licor_fco2(inds(v):indt(v)));
        
        indslagrow(v)=inds(v)-avglag(v); % give the indexes for bin start and end points for plotting below 
        indtlagrow(v)=indt(v)-avglag(v);
        
        avgWindspeed(v)=mean(Windspeedabs_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));%loop through the indexes to determine the average temperature between indexes
        avgLatitude(v)=mean(Latitude_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));%loop through the indexes to determine the average temperature between indexes
        avgLongitude(v)=mean(Longitude_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));%loop through the indexes to determine the average temperature between indexes
        avgtide(v)=mean(Plymouth_tidal_height_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));%loop through the indexes to determine the average temperature between indexes

        avgEqu_temp(v)=mean(Equilibrator_Temperature_interp((inds(v)):(indt(v))));
        stdEqu_temp(v)=std(Equilibrator_Temperature_interp(inds(v):indt(v))); 
        semEqu_temp(v)=sem(Equilibrator_Temperature_interp(inds(v):indt(v))); 
        
    end
 
    %mooring data
    C=load('C:\Users\rps207\Documents\MATLAB\2016 - Temperature mooring output analysis\Data\NSTempmooringl4.mat');
    load('C:\Users\rps207\Documents\MATLAB\2016 - Temperature mooring output analysis\Data\NSTempmooringl4.mat','T_DT');
    names=({'T_03m','T_06m','T_15m','T_24m','T_35m','T_star_29m'})';
    %create new names for variables with _interp suffix
    
        %create new names for variables with _interp suffix
    for t=1:numel(names)
        rty=[names(t) '_interp'];x=horzcat(rty{:}) ;new_names(:,t)=cellstr(x);
    end
    % Interpolate variables in a loop using new variable names
    for i=1:numel(names)
        interpvars.(new_names{i})=interp1(T_DT,C.(names{i}),Licor_Datetime);
    end
    v2struct(interpvars)

    for v=1:length(inds) %loop through the indexes to determine the average depths between indexes
        avgT_03m(v)=mean(T_03m_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        stdT_03m(v)=std(T_03m_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        semT_03m(v)=semmoor(T_03m_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
      
        avgT_06m(v)=mean(T_06m_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        stdT_06m(v)=std(T_06m_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        semT_06m(v)=semmoor(T_06m_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        
        avgT_15m(v)=mean(T_15m_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        stdT_15m(v)=std(T_15m_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        semT_15m(v)=semmoor(T_15m_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        
        avgT_35m(v)=mean(T_35m_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        stdT_35m(v)=std(T_35m_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        semT_35m(v)=semmoor(T_35m_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        
        avgT_24m(v)=mean(T_24m_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        stdT_24m(v)=std(T_24m_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        semT_24m(v)=semmoor(T_24m_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));

        avgT_29m(v)=mean(T_star_29m_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        stdT_29m(v)=std(T_star_29m_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));
        semT_29m(v)=semmoor(T_star_29m_interp((inds(v)-avglag(v)):(indt(v)-avglag(v))));

        
    end

    
    depth03m=ones(length(avgT_03m),1)*0.3;
    depth06m=ones(length(avgT_06m),1)*0.6;
    depth15m=ones(length(avgT_15m),1)*1.5;
    depth35m=ones(length(avgT_35m),1)*3.5;
    depth24m=ones(length(avgT_24m),1)*2.4;
    depth29m=ones(length(avgT_29m),1)*2.9;
    
    indslag=indslagrow';
    indtlag=indtlagrow';
    
    profdt=Licor_Datetime(inds(1)-avglag(1)) ;
    avgTstr=avgtime-profdt;
        
    avgpresdb = sw_pres(avgdepth,max(Latitude_interp(1)));
    avgdensity=sw_dens(avgsal,avgtemp,avgpresdb);

    

    
    %%%%%%%%%%%%%%%%%%%%%% error propogation exercise %%%%%%%%%%%
    
    %For FCO2 to FCO2 insitu
    %error for subtractions p=A-B 
    %std_p=sqrt(A_sigma^2 + B_sigma^2)
    errminustemps=sqrt((stdtemp.^2) + ((stdEqu_temp).^2)); 
    
    %error for multiplication with constant c=kp
    %std_c=k*std_p
    errmultconstant=0.0423*errminustemps;
    
    % error for exponential t=exp(c)
    % std_t = t * c_std
    t=exp(0.0423.*(avgtemp - avgEqu_temp));
    errexp=(t.*errmultconstant);
    
    % error for multiplication z=F*t
    %std_z = sqrt((F_sigma/F_mu)^2 + (t_sigma/t_mu)^2)*F_mu*t_mu  %%noteF_mu*t_mu = z
    stdavgfco2surf = sqrt((stdfco2/avgfco2)^2 + (errexp/t)^2)*avgfco2surf;


    %need to propogate error for this
    %For FCO2 to FCO2 insitu
    %error for subtractions p=A-B 
    %std_p=sqrt(A_sigma^2 + B_sigma^2)
    
    errminustemps=sqrt((semtemp.^2) + ((semEqu_temp).^2)); 
    
    %error for multiplication with constant c=kp
    %std_c=k*std_p
    errmultconstant=0.0423*errminustemps;
    
    % error for exponential t=exp(c)
    % std_t = t * c_std
    t=exp(0.0423.*(avgtemp - avgEqu_temp));
    errexp=(t.*errmultconstant);
    
    % error for multiplication z=F*t
    %std_z = sqrt((F_sigma/F_mu)^2 + (t_sigma/t_mu)^2)*F_mu*t_mu  %%noteF_mu*t_mu = z
    semavgfco2surf = sqrt((semfco2/avgfco2)^2 + (errexp/t)^2)*avgfco2surf;
    

    
    
    ty=[avgdepth ];
    
    
    S=regstats([avgfco2surf],ty,'linear');
    
    load('Data/underwaypCO2Seasonalstudy.mat','underwayfCO2_sw','underway_DT','underwayLat','underwayLong','underwayS_sea','underwayT_sea','underwaypCO2_atm')
    load('C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\Plymouth_tides.mat')
    
    %underway_DT is in UTC everything is in local time , convert utc to
    %local time
    
     for t=1:length(underwayLat)
        underway_dist(t,:)= pos2dist(50.25, -4.217,underwayLat(t),underwayLong(t),2);
     end
     
     
                 undtimeavg=mean(Licor_Datetime(undind3:undind4));

      undtimeavgret=mean(Licor_Datetime(undind3:undind4));
       [ ind_temp_lowtide_und,~  ]=find(Plymouth_tidal_height_dt(ind_low_tide)<undtimeavgret);
       nearest_lowtide_before_und=Plymouth_tidal_height_dt(ind_low_tide(max(ind_temp_lowtide_und)));
       dt_since_lowtide_und_ret=undtimeavgret-nearest_lowtide_before_und;
       dt_since_lowtide_und_retHHMM=datestr(dt_since_lowtide_und_ret,'HH:MM');


       undtimeavgoutbound=mean(Licor_Datetime(undind1:undind2));
       [ ind_temp_lowtide_und,~  ]=find(Plymouth_tidal_height_dt(ind_low_tide)<undtimeavgoutbound);
       nearest_lowtide_before_und=Plymouth_tidal_height_dt(ind_low_tide(max(ind_temp_lowtide_und)));
       dt_since_lowtide_und_outbound=undtimeavgoutbound-nearest_lowtide_before_und;
       dt_since_lowtide_und_outboundHHMM=datestr(dt_since_lowtide_und_outbound,'HH:MM');

       
       
       
       
       
       
       
  
   

           %find points of CO2 system overlap
    
    %points of overlap on station
    NSOP_overlap=[];  shower_overlap=[];und_dist_L4=[];co2indall=[];
    %unique string sequence all days when want to compare NSOP data
    if strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-06-30')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-08-04')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-15')|strcmpi(fList(j,1:10),'2016-09-21')
        [x y ]=find(underway_DT>Licor_Datetime(inds(1)) & underway_DT<Licor_Datetime(indt(end)));
        yyyyyy=[];
        for l=1:length(x);
            [r s]=min(abs(underway_DT(x(l))-Licor_Datetime));
            yyyyyy=[yyyyyy,s];
        end
        
        
        
        for jkl=1:length(yyyyyy); 
        if CTDdepth_interp(yyyyyy(jkl))>2.5 && CTDdepth_interp(yyyyyy(jkl))<3.5
            NSOP_overlap=[NSOP_overlap ;Licor_fco2_combined(yyyyyy(jkl))];
            und_dist_L4=[und_dist_L4;dist(yyyyyy(jkl))];
            f=x(jkl)
            shower_overlap=[shower_overlap;underwayfCO2_sw(f)];
        end
        end
    else
    end
    
    
    
    %points of overlap during underway voyages
     NSOP_overlap_und=[];  shower_overlap_und=[];und_dist_L4_und=[];Licor_Datetime_und_dist_L4_und=[];

    %outward voyages
    d=[];e=[];w=[];z=[];
    emp=isempty(undind1);
    if emp==1;
    else
    [d e ]=find(underway_DT>Licor_Datetime(undind1) & underway_DT<Licor_Datetime(undind2));
    end
    emp3=isempty(d);
     %find out if the overlap is empty, then find the index in licor
     %datetime for overlap
     if emp3==1
     else
     yyyyyy=[];
    for l=1:length(d);
    [r s]=min(abs(underway_DT(d(l))-Licor_Datetime));
    yyyyyy=[yyyyyy,s];
    end
     end
         
    %return voyages
    emp2=isempty(undind3);
    if emp2==1;
    else
    [w z ]=find(underway_DT>Licor_Datetime(undind3) & underway_DT<Licor_Datetime(undind4+underwaylag));
    end
    emp4=isempty(w);
     if emp4==1
     else
    xxxxx=[];
    for l=1:length(w);
    [r s]=min(abs(underway_DT(w(l))-Licor_Datetime));
    xxxxx=[xxxxx,s];
    end
     end
     
    
    %unique string sequence all days when want to compare NSOP data
    if strcmpi(fList(j,1:10),'2016-06-10')|strcmpi(fList(j,1:10),'2016-06-15')|strcmpi(fList(j,1:10),'2016-06-22')|strcmpi(fList(j,1:10),'2016-06-30')|strcmpi(fList(j,1:10),'2016-07-07')|strcmpi(fList(j,1:10),'2016-07-13')|strcmpi(fList(j,1:10),'2016-07-20')|strcmpi(fList(j,1:10),'2016-08-04')|strcmpi(fList(j,1:10),'2016-08-17')|strcmpi(fList(j,1:10),'2016-08-24')|strcmpi(fList(j,1:10),'2016-09-15')|strcmpi(fList(j,1:10),'2016-09-21')
     
    if emp3==1 & emp3==1;
        %both empty and do nothing
    elseif emp3==1 & emp4==0;
    shower_overlap_und=[underwayfCO2_sw(d)];
    NSOP_overlap_und=[Licor_fco2_combined(yyyyyy)];  
    und_dist_L4_und=[dist(yyyyyy)];
    Licor_Datetime_und_dist_L4_und=[Licor_Datetime(yyyyyy)];
    elseif emp3==0 & emp4==1;
    shower_overlap_und=[underwayfCO2_sw(w)];
    NSOP_overlap_und=[Licor_fco2_combined(xxxxx)];  
     und_dist_L4_und=[dist(xxxxx)];
    Licor_Datetime_und_dist_L4_und=[Licor_Datetime(xxxxx)]
    else
    shower_overlap_und=[underwayfCO2_sw(d) ; underwayfCO2_sw(w)];
    NSOP_overlap_und=[Licor_fco2_combined(yyyyyy) ; Licor_fco2_combined(xxxxx)];  
    und_dist_L4_und=[dist(yyyyyy) ; dist(xxxxx)];
    Licor_Datetime_und_dist_L4_und=[Licor_Datetime(yyyyyy) ; Licor_Datetime(xxxxx)];
    end
    else
    end
     
     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %Calculate underway fluxes point by point
    for q=1:length(Windspeedabs_interp);
    %equations
    
    %%this gives conversation factor in takhashi 2009
    % (10^-2)*24*(365/12) * (12*(10^3)*(10^-6)) *0.26 *(660^0.5)
    
    %this does point by point wind speed - don't use
    %k660=(0.222*((Windspeedabs_interp(k))^2) +0.333*(Windspeedabs_interp(k)));%from nightingale 2000 with units of cm hr-1
   
    %average wind over return transect and make whole transects
    windspeedmeantransect=mean(Windspeedabs_interp(undind3:undind4));
    windspeedmeantransmatrix=ones(length(Windspeedabs_interp),1)*windspeedmeantransect;
    k660(q)=(0.222*((windspeedmeantransmatrix(q))^2) +0.333*(windspeedmeantransmatrix(q)));%from nightingale 2000 with units of cm hr-1
    %same for ATM pCO2
    ARMCO2meantransect=mean(underwaypCO2_atm_interp(undind3:undind4));
    ATMSCO2meantransmatrix=ones(length(Windspeedabs_interp),1)*ARMCO2meantransect;

    Schcw(q)=(2073.1 -(125.62*SeaTemp_interp(q)) + (3.6276*((SeaTemp_interp(q))^2)) - (0.04321*((SeaTemp_interp(q))^3)));%from wanniqhoif 1992 - for seawater
    Schdep(q)=(Schcw(q)/660)^-0.5; %johnson 2010
    Kw(q)=Schdep(q)*k660(q);%k with units of cm hr-1
    k0(q)=exp(-60.2409 + 93.4517*(100/(((SeaTemp_interp(q))+273.15))) + 23.3585*log((((SeaTemp_interp(q))+273.15))/100) +(Salinity_interp(q))*(0.023517 -0.023656*((((SeaTemp_interp(q))+273.15))/100) + 0.0047036*(((((SeaTemp_interp(q))+273.15))/100)^2)));%K0 with units(mol l-1 atm-1) weiss 1974
    scal=(12*1000)/(1e6);% scaling K0- sclaing factor to convert from mol L-1 atm-1 to gCm_3 micro atm-1
    scal2=(1/100); %scaling Kw from cm hr-1 to m hr-1
    %transfer coefficent
    TR(q)=Kw(q)*k0(q)*scal*scal2; %Units of gC m-2 hr-1 atm-1
    TRmo(q)=TR(q)*24*365/12;%Units of gC m-2 month-1 atm-1
    %calculate flux using delta co2f
    DELPCO2_TRANS(q)=(Licor_fco2_surfaceund(q)-ATMSCO2meantransmatrix(q));
    Fhr(q)=TR(q)*DELPCO2_TRANS(q); %Units of gC m-2 hr-1  - NOTE WHEN multiplied by 365*24 its about 10g which is very reasonable
    Fmo(q)=TRmo(q)*(Licor_fco2_combined(q)-ATMSCO2meantransmatrix(q)); %Units of gC m-2 mo-1
    %flux in mmol m2 hr-1
    Fmolhrx(q)=Fhr(q)*(1/12)*1e4;%units in mmol C m2 hr-1;
    end
    
    
    
    
    
   %SHOWERHEAD FLUX all points
  if emp3==1;
          for k=1:length(underwayfCO2_sw(w));

    %if multiple pco2 points during underway journey back average them
    shower_overlap_und_only=[underwayfCO2_sw(w)];
    showerSAL_overlap_und_only=[underwayS_sea(w)];
    showertemp_overlap_und_only=[underwayT_sea(w)];
    showeratm_overlap_und_only=[underwaypCO2_atm(w)];       
        
    showerlat_overlap_und_only=[underwayLat(w)];
    showerlong_overlap_und_only=[underwayLong(w)];    
    showerdt_overlap_und_only=[underway_DT(w)];
    
    %average wind over return transect
    windspeedmeantransect=mean(Windspeedabs_interp(undind3:undind4));
    k660=(0.222*((windspeedmeantransect)^2) +0.333*(windspeedmeantransect));%from nightingale 2000 with units of cm hr-1
              

 Schcw=(2073.1 -(125.62*showertemp_overlap_und_only(k)) + (3.6276*((showertemp_overlap_und_only(k))^2)) - (0.04321*((showertemp_overlap_und_only(k))^3)));%from wannikhoif 1992 - for seawater
    Schdep=(Schcw/660)^-0.5; %johnson 2010
    Kw=Schdep*k660;%k with units of cm hr-1
    k0=exp(-60.2409 + 93.4517*(100/(((showertemp_overlap_und_only(k))+273.15))) + 23.3585*log((((showertemp_overlap_und_only(k))+273.15))/100) +(showerSAL_overlap_und_only(k))*(0.023517 -0.023656*((((showertemp_overlap_und_only(k))+273.15))/100) + 0.0047036*(((((showertemp_overlap_und_only(k))+273.15))/100)^2)));%K0 with units(mol l-1 atm-1) weiss 1974
    scal=(12*1000)/(1e6);% scaling K0- sclaing factor to convert from mol L-1 atm-1 to gCm_3 micro atm-1
    scal2=(1/100); %scaling Kw from cm hr-1 to m hr-1
    %transfer coefficent
    TR=Kw*k0*scal*scal2; %Units of gC m-2 hr-1 atm-1
    TRmo=TR*24*365/12;%Units of gC m-2 month-1 atm-1
    %calculate flux using delta co2f
    Fhr=TR*(shower_overlap_und_only(k)-showeratm_overlap_und_only(k)); %Units of gC m-2 hr-1  - NOTE WHEN multiplied by 365*24 its about 10g which is very reasonable
    Fmo=TRmo*(shower_overlap_und_only(k)-showeratm_overlap_und_only(k)); %Units of gC m-2 mo-1

    %flux in mmol m2 hr-1
    Flux_showerhead(k)=Fhr*(1/12)*1e4;%units in mmol C m2 hr-1;
    
          end
  else
  end
    
   
    
    
    
    
    
    
    
    
     
    if pFigs == 1
        
        %Figures
        %List of figures
        %Figure 1 - Depth profile with indexed points plotted and smoothed differential plot
        %Figure 2 - Timeseries subplots of data(depth,temp,co2,salinity), cropped to the length of the profile
        %Figure 3 - Subplots of profiles. temp,sal.fco2 and fco2 corrected
        %Figure 4 -
        %Figure 5 -        
        
        Temperature_probe_a_interp = interp1(Datetime_probe_a,Temperature_probe_a,Licor_Datetime);l = find(isnan(Temperature_probe_a_interp ))'; Temperature_probe_a_interp (l) = 0;
        Temperature_probe_b_interp = interp1(Datetime_probe_b,Temperature_probe_b,Licor_Datetime);m = find(isnan(Temperature_probe_b_interp ))'; Temperature_probe_b_interp (m) = 0;
        profdt=Licor_Datetime(inds(1)-avglag(1)) ;
        avgTstr=avgtime-profdt;
        
        
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
        

        
        %Plot depth profile and index points
       %Plot depth profile and index points
        figure(1);
        set(gcf, 'Color','w','Position', get(0,'Screensize'));
        subplot(1,2,1); 
        %%
        plot(Licor_Datetime(ind(1):indt(end)),CTDdepth_interp(ind(1):indt(end)),'.');  hold on;
        plot(Licor_Datetime(ind(1):indt(end)),smoothCTD(ind(1):indt(end)),'.r');         
        plot(Licor_Datetime(ind),smoothCTD(ind),'.g','MarkerSize',32)
        plot(Licor_Datetime(indr),smoothCTD(indr),'.k','MarkerSize',32)
        plot(Licor_Datetime(indr),smoothCTD(indr),'.y','MarkerSize',20)
        dynamicDateTicks([], [], ' HH:MM');
        setDateAxes(gca,'Fontsize',24);
        set(gca,'YDir','reverse')
        indtextlabel1 = (1:length(ind))'; indtextlabel2 = num2str(indtextlabel1); indtextlabel3= cellstr(indtextlabel2);
        text(Licor_Datetime(ind+300), smoothCTD(ind), indtextlabel3,'color','k','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        ylabel('Depth (m)','Fontsize',24); xlabel('Time','Fontsize',24);
        legend('original','filtered','rejected midpoints','','midpoints used analysis','Location','SouthEast');
        hold off;
        %%
        subplot(1,2,2);
        %%
        h1=plot(Licor_Datetime(ind(1):indt(end)),diffCTD(ind(1):indt(end)),'.');
        h2=line([Licor_Datetime(ind(1)),Licor_Datetime(indt(end))],[threshold,threshold],'Color','r','LineWidth',2);
        h3=line([Licor_Datetime(ind(1)),Licor_Datetime(indt(end))],[-1*threshold,-1*threshold],'Color','r','LineWidth',2);
        legend(h2,'Bin identification threshold');
        dynamicDateTicks([], [], 'dd/mm HH:MM');
        setDateAxes(gca,'Fontsize',24);
        ylabel('Differential filtered data (m)'); xlabel('Time');
        set(gcf, 'Color','w','Position', get(0,'Screensize'));        
        if ~exist('Figures/Depthbins', 'dir')
                 mkdir('Figures','Depthbins');
        end
        
        mtit(['Depth Bin Identification' ' '  fList(j,1:10)],'fontsize',24,'xoff',-.025,'yoff',.025);
        if ~exist('Figures/Depthbins', 'dir')
             mkdir('Figures/Depthbins');
        end
        saveas(gca,[pwd '/Figures/Depthbins/' fList(j,1:10) '.png']);
        hold off;
        %%

       figure(2) %tidier depth profile
       %%
       %Plot depth profile and index points
       set(gcf, 'Color','w','Position', get(0,'Screensize'));
       plot(Licor_Datetime(ind(1):indt(end)),CTDdepth_interp(ind(1):indt(end)),'.');  hold on;
       %plot(Licor_Datetime(ind(1):indt(end)),smoothCTD(ind(1):indt(end)),'.r');
       plot(Licor_Datetime(ind),smoothCTD(ind),'.g','MarkerSize',32)
       dynamicDateTicks([], [], ' HH:MM');
       set(gca,'FontSize',16);
       set(gca,'FontSize',16);
       set(gca,'YDir','reverse')
       indtextlabel1 = (1:length(ind))'; indtextlabel2 = num2str(indtextlabel1); indtextlabel3= cellstr(indtextlabel2);
       text(Licor_Datetime(ind+300), smoothCTD(ind), indtextlabel3,'color','k','Fontsize',12,'Fontweight','bold','BackgroundColor','w');
       ylabel('Depth (m)','Fontsize',16); xlabel('Time','Fontsize',16);
       %%
       
        figure (3); %timseries of co2,temp etc
        set(gcf, 'Color','w','Position', get(0,'Screensize'));
        subplot(5,1,1) % Depth timeseries
        %%
        %plot patch first underneath.
        for r=1:length(inds)
        x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
        y= [ max(CTDdepth_interp) max(CTDdepth_interp) 0 0];
        patch(Licor_Datetime(x),y,colour_greyshade); 
        end
        hold on
        % plot main part of plot
        plot(Licor_Datetime(ind(1):indt(end)),CTDdepth_interp(ind(1):indt(end)),'LineWidth',2);
        title('Depth timeseires Plot','FontSize',12);
        xlabel('Time','FontSize',12); 
        ylabel('Depth (m)','FontSize',12);
        set(gca,'FontSize',12)
        dynamicDateTicks([], [], ' HH:MM');
        hold on
        %add text labels isntead of legend and before vertical lines drawn
        textbp('Depth','color','b','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        textbp('Bin Start','color','r','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        textbp('Bin End','color',colour_green,'Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        %plot bin start and end points and numerical labels
        plot(Licor_Datetime([indslag(1:end),indslag(1:end)]),[0,max(CTDdepth_interp)],'Color','r','LineWidth',2);
        plot(Licor_Datetime([indtlag(1:end),indtlag(1:end)]),[0,max(CTDdepth_interp)],'Color',colour_green,'LineWidth',2);
        indtextlabel1 = (1:length(inds))'; indtextlabel2 = num2str(indtextlabel1); indtextlabel3= cellstr(indtextlabel2);
        if ptxtnum== 1
        text(Licor_Datetime(indslag), (ones(length(inds),1)*max(CTDdepth_interp)), indtextlabel3,'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        text(Licor_Datetime(indtlag), (ones(length(inds),1)*max(CTDdepth_interp)), indtextlabel3,'color',colour_green,'Fontsize',12,'Fontweight','bold','BackgroundColor','w');         
        end; hold off;
        %%   
        subplot(5,1,2) %Insitu temperature timeseries
        %%
        %plot patch first underneath.
        for r=1:length(inds)
        x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
        y= [ max(Equilibrator_Temperature_interp) max(Equilibrator_Temperature_interp) 0 0];
        patch(Licor_Datetime(x),y,colour_greyshade); 
        end
        dynamicDateTicks([], [], ' HH:MM');
        title('Temperature timeseires plot','FontSize',12);
        xlabel('Time','FontSize',12); 
        ylabel({['CTD' char(10) 'Temperature (',num2str(degree_symbol),'C)']},'FontSize',12);       
        %plot bin start and end points and numerical labels
        mmequm2std=mean(Equilibrator_Temperature_interp(inds(1):indt(end)))-(2*std(Equilibrator_Temperature_interp(inds(1):indt(end))));
        mnequp2std=mean(Equilibrator_Temperature_interp(inds(1):indt(end)))+(2*std(Equilibrator_Temperature_interp(inds(1):indt(end))));
        %Plot ctd temp
        mnCTDtempm2std=mean(CTDTemperature_interp(ind(1):indt(end)))-(2*std(CTDTemperature_interp(ind(1):indt(end))));
        mnCTDtempp2std=mean(CTDTemperature_interp(ind(1):indt(end)))+(2*std(CTDTemperature_interp(ind(1):indt(end))));
%         set(gca, 'YLim', [mnCTDtempm2std,mnCTDtempp2std])
        set(gca,'FontSize',12); hold on;
        %add CTD temperature to 2nd axis
        plot(Licor_Datetime(ind(1):indt(end)),CTDTemperature_interp(ind(1):indt(end)),'m','LineWidth',2); 
        hold on
        plot(Licor_Datetime(ind(1):indt(end)),Temperature_probe_a_interp_lag(ind(1):indt(end)),'r','LineWidth',2); 
        plot(Licor_Datetime(ind(1):indt(end)),Temperature_probe_b_interp_lag(ind(1):indt(end)),'b','LineWidth',2); 

        %add text labels instead of a legend
        textbp( 'CTD','color','m','Fontsize',10,'Fontweight','bold','BackgroundColor','w');  hold on;
        if pequT== 1
        mmequm2std=mean(Equilibrator_Temperature_interp(inds(1):indt(end)))-(2*std(Equilibrator_Temperature_interp(inds(1):indt(end))));
        mnequp2std=mean(Equilibrator_Temperature_interp(inds(1):indt(end)))+(2*std(Equilibrator_Temperature_interp(inds(1):indt(end))));
        addaxis(Licor_Datetime(ind(1):indt(end)),Equilibrator_Temperature_interp_lag(ind(1):indt(end)),[mmequm2std mnequp2std],'Color',colour_green,'LineWidth',2); 
        hold on;
        addaxisplot(Licor_Datetime(ind(1):indt(end)),Temperature_probe_a_interp_lag(ind(1):indt(end)) ,2,'-yo','LineWidth',2); 
        addaxisplot(Licor_Datetime(ind(1):indt(end)),Temperature_probe_b_interp_lag(ind(1):indt(end)) ,2,'-bo','LineWidth',2);
        % add text labels instead of a legend
        textbp( 'Equilibrator' ,'Color',colour_green,'Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        textbp( 'Probe A','Color','y','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        textbp( 'Probe B','color','b','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        addaxislabel(2,['Equilibrator Temperature (',num2str(degree_symbol),'C)']);
        end
        %plot bin start and end points and numerical labels
        addaxisplot(Licor_Datetime([indslag(1:end),indslag(1:end)]),[0,max(Equilibrator_Temperature_interp)],1,'Color','r','LineWidth',2);
        addaxisplot(Licor_Datetime([indtlag(1:end),indtlag(1:end)]),[0,max(Equilibrator_Temperature_interp)],1,'Color',colour_green,'LineWidth',2);
        if ptxtnum== 1
        text(Licor_Datetime(indslag), (ones(length(inds),1)*mnCTDtempp2std), indtextlabel3,'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        text(Licor_Datetime(indtlag), (ones(length(inds),1)*mnCTDtempp2std), indtextlabel3,'color',colour_green,'Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        end;
        hold off;
              
%          load (['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\Underwayinterp\' fList(j,1:10) '.mat'],'ShipSST_interp')
%         set(gca,'FontSize',12) 
        %%
        subplot(5,1,3)% Licor timeseries
        %%
        %plot patch first underneath.
        for r=1:length(inds)
        x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
        y= [ max(Licor_fco2_surface) max(Licor_fco2_surface) 0 0];
        patch(Licor_Datetime(x),y,colour_greyshade); 
        end
        hold on
        %
        plot(Licor_Datetime(ind(1):indt(end)),LicorxCO2Bcal_lag(ind(1):indt(end)),'LineWidth',2); hold on;
        plot(Licor_Datetime(ind(1):indt(end)),Licor_fco2_lag(ind(1):indt(end)),'m','LineWidth',2); hold on;
        plot(Licor_Datetime(ind(1):indt(end)),Licor_fco2_surface_lag(ind(1):indt(end)),'k','LineWidth',2); 
        plot(Licor_Datetime(ind(1):indt(end)),underwaypCO2_sw_interp(ind(1):indt(end)),'g','LineWidth',2); 
        set(gca,'FontSize',12)
        dynamicDateTicks([], [], 'HH:MM');
        %add text labels
        textbp( 'XCO_2','color','b','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        textbp( 'F CO_2','color','m','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        textbp( 'F CO_2 underway','color','m','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        textbp( ['F CO_2' char(10) 'insitu corrected'],'color','k','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        title('CO_2 Timeseires Plot','FontSize',12);xlabel('Time','FontSize',12); ylabel(['CO_2' char(10) '(',num2str(micro_symbol),'atm)'],'FontSize',12);
        mnco2m2std=mean(Licor_fco2_surface(inds(1):indt(end)))-(2*std(Licor_fco2_surface(inds(1):indt(end))));
        mnco2p2std=mean(Licor_fco2_surface(inds(1):indt(end)))+(2*std(Licor_fco2_surface(inds(1):indt(end))));
%         set(gca, 'YLim', [mnco2m2std,mnco2p2std]);
        set(gca,'FontSize',12); 
        %plot bin start and end points and numerical labels
        plot(Licor_Datetime([indslag(1:end),indslag(1:end)]),[0,max(Licor_fco2_surface)],'Color','r','LineWidth',2);
        plot(Licor_Datetime([indtlag(1:end),indtlag(1:end)]),[0,max(Licor_fco2_surface)],'Color',colour_green,'LineWidth',2);
        if ptxtnum== 1
        text(Licor_Datetime(indslag), (ones(length(inds),1)*mnco2p2std), indtextlabel3,'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        text(Licor_Datetime(indtlag), (ones(length(inds),1)*mnco2p2std), indtextlabel3,'color',colour_green,'Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        end; hold off;
        %%
        subplot(5,1,4)%Salinity timseries
        %%
        %plot patch first underneath.
        for r=1:length(inds)
        x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
        y= [ max(CTDsalinity_interp) max(CTDsalinity_interp) 0 0];
        patch(Licor_Datetime(x),y,colour_greyshade); 
        end
        hold on
        plot(Licor_Datetime(ind(1):indt(end)),CTDsalinity_interp(ind(1):indt(end)),'LineWidth',2);
        dynamicDateTicks([], [], ' HH:MM');
        title('Salinity timeseires plot','FontSize',12);xlabel('Time','FontSize',12); ylabel(['Salinity' char(10) '(PSU)'],'FontSize',12); set(gca,'FontSize',12); hold on;
        %set axis limits as 2 stnd devs
        mnsalm2std=mean(CTDsalinity_interp(inds(1):indt(end)))-(2*std(CTDsalinity_interp(inds(1):indt(end))));
        mnsalp2std=mean(CTDsalinity_interp(inds(1):indt(end)))+(2*std(CTDsalinity_interp(inds(1):indt(end))));
        set(gca, 'YLim', [mnsalm2std,mnsalp2std])
        textbp( 'Salinity','color','b','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        %PLot bin divisers and labels
        plot(Licor_Datetime([indslag(1:end-1),indslag(1:end-1)]),[0,max(CTDsalinity_interp)],'Color','r','LineWidth',2);
        plot(Licor_Datetime([indtlag(2:end),indtlag(2:end)]),[0,max(CTDsalinity_interp)],'Color',colour_green,'LineWidth',2);
        if ptxtnum== 1
        text(Licor_Datetime(indslag), (ones(length(inds),1)*mnsalp2std), indtextlabel3,'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        text(Licor_Datetime(indtlag), (ones(length(inds),1)*mnsalp2std), indtextlabel3,'color',colour_green,'Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        end; hold off;
        %%
        subplot(5,1,5)%Flow rate timseries
        %%
        %plot patch first underneath.
        for r=1:length(inds)
        x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
        y= [ max(MeasuredFlow_interp) max(MeasuredFlow_interp) 0 0];
        patch(Licor_Datetime(x),y,colour_greyshade); 
        end
        hold on
        plot(Licor_Datetime(ind(1):indt(end)),MeasuredFlow_interp(ind(1):indt(end)),'LineWidth',2); hold on;
        plot(Licor_Datetime(ind(1):indt(end)),smoothflow(ind(1):indt(end)),'m','LineWidth',2);
        dynamicDateTicks([], [], ' HH:MM');
        title('Flowrate timseries plot','FontSize',12);xlabel('Time','FontSize',12); ylabel(['Flowrate' char(10) 'L/min'],'FontSize',12);set(gca,'FontSize',12);
        %set axis limits as 2 stnd devs
        mnflowm2std=mean(smoothflow(inds(1):indt(end)))-(2*std(smoothflow(inds(1):indt(end))));
        mnflowp2std=mean(smoothflow(inds(1):indt(end)))+(2*std(smoothflow(inds(1):indt(end))));
        set(gca, 'YLim', [mnflowm2std,mnflowp2std])
        textbp( 'Smoothed Flowrate','color','m','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        %PLot bin divisers and labels
        plot(Licor_Datetime([indslag(1:end-1),indslag(1:end-1)]),[0,max(CTDsalinity_interp)],'Color','r','LineWidth',2);
        plot(Licor_Datetime([indtlag(2:end),indtlag(2:end)]),[0,max(CTDsalinity_interp)],'Color',colour_green,'LineWidth',2);
        if ptxtnum== 1
        text(Licor_Datetime(indslag), (ones(length(inds),1)*mnflowp2std), indtextlabel3,'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        text(Licor_Datetime(indtlag), (ones(length(inds),1)*mnflowp2std), indtextlabel3,'color',colour_green,'Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        end

        mtit(['Timeseries for' ' ' fList(j,1:10) ],'fontsize',24,'xoff',-.025,'yoff',.025);

        if ~exist('Figures/Timeseries', 'dir')
             mkdir('Figures/Timeseries');
        end
        
        saveas(gca,[pwd '/Figures/Timeseries/' fList(j,1:10) '.png']);
        hold off;
        %%
        
        figure (300001); %timseries of co2,temp etc
        set(gcf, 'Color','w','Position', get(0,'Screensize'));
        subplot(5,1,1) % Depth timeseries
        %%
        %plot patch first underneath.
        for r=1:length(inds)
        x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
        y= [ max(CTDdepth_interp) max(CTDdepth_interp) 0 0];
        patch(Licor_Datetime(x),y,colour_greyshade); 
        end
        hold on
        % plot main part of plot
        plot(Licor_Datetime,CTDdepth_interp,'LineWidth',2);
        title('Depth timeseires Plot','FontSize',12);
        xlabel('Time','FontSize',12); 
        ylabel('Depth (m)','FontSize',12);
        set(gca,'FontSize',12)
        dynamicDateTicks([], [], ' HH:MM');
        hold on
        %add text labels isntead of legend and before vertical lines drawn
        textbp('Depth','color','b','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        textbp('Bin Start','color','r','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        textbp('Bin End','color',colour_green,'Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        %plot bin start and end points and numerical labels
        plot(Licor_Datetime([indslag(1:end),indslag(1:end)]),[0,max(CTDdepth_interp)],'Color','r','LineWidth',2);
        plot(Licor_Datetime([indtlag(1:end),indtlag(1:end)]),[0,max(CTDdepth_interp)],'Color',colour_green,'LineWidth',2);
        indtextlabel1 = (1:length(inds))'; indtextlabel2 = num2str(indtextlabel1); indtextlabel3= cellstr(indtextlabel2);
        if ptxtnum== 1
        text(Licor_Datetime(indslag), (ones(length(inds),1)*max(CTDdepth_interp)), indtextlabel3,'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        text(Licor_Datetime(indtlag), (ones(length(inds),1)*max(CTDdepth_interp)), indtextlabel3,'color',colour_green,'Fontsize',12,'Fontweight','bold','BackgroundColor','w');         
        end; hold off;
        %%   
        subplot(5,1,2) %Insitu temperature timeseries
        %%
        %plot patch first underneath.
        for r=1:length(inds)
        x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
        y= [ max(Equilibrator_Temperature_interp) max(Equilibrator_Temperature_interp) 0 0];
        patch(Licor_Datetime(x),y,colour_greyshade); 
        end
        dynamicDateTicks([], [], ' HH:MM');
        title('Temperature timeseires plot','FontSize',12);
        xlabel('Time','FontSize',12); 
        ylabel({['CTD' char(10) 'Temperature (',num2str(degree_symbol),'C)']},'FontSize',12);       
        %plot bin start and end points and numerical labels
        mmequm2std=mean(Equilibrator_Temperature_interp(inds(1):indt(end)))-(2*std(Equilibrator_Temperature_interp(inds(1):indt(end))));
        mnequp2std=mean(Equilibrator_Temperature_interp(inds(1):indt(end)))+(2*std(Equilibrator_Temperature_interp(inds(1):indt(end))));
        %Plot ctd temp
        mnCTDtempm2std=mean(CTDTemperature_interp(ind(1):indt(end)))-(2*std(CTDTemperature_interp(ind(1):indt(end))));
        mnCTDtempp2std=mean(CTDTemperature_interp(ind(1):indt(end)))+(2*std(CTDTemperature_interp(ind(1):indt(end))));
%         set(gca, 'YLim', [mnCTDtempm2std,mnCTDtempp2std])
        set(gca,'FontSize',12); hold on;
        %add CTD temperature to 2nd axis
        plot(Licor_Datetime,CTDTemperature_interp,'m','LineWidth',2); 
        hold on
        plot(Licor_Datetime,Temperature_probe_a_interp_lag,'r','LineWidth',2); 
        plot(Licor_Datetime,Temperature_probe_b_interp_lag,'b','LineWidth',2); 
        plot(Licor_Datetime,underwayT_sea_interp,'g','LineWidth',2); 
        plot(Licor_Datetime,SeaTemp_interp,'k','LineWidth',2); 

        %add text labels instead of a legend
        textbp( 'CTD','color','m','Fontsize',10,'Fontweight','bold','BackgroundColor','w');  hold on;

        %plot bin start and end points and numerical labels
        addaxisplot(Licor_Datetime([indslag(1:end),indslag(1:end)]),[0,max(Equilibrator_Temperature_interp)],1,'Color','r','LineWidth',2);
        addaxisplot(Licor_Datetime([indtlag(1:end),indtlag(1:end)]),[0,max(Equilibrator_Temperature_interp)],1,'Color',colour_green,'LineWidth',2);
        if ptxtnum== 1
        text(Licor_Datetime(indslag), (ones(length(inds),1)*mnCTDtempp2std), indtextlabel3,'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        text(Licor_Datetime(indtlag), (ones(length(inds),1)*mnCTDtempp2std), indtextlabel3,'color',colour_green,'Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        end;
        hold off;
              
%          load (['C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\Underwayinterp\' fList(j,1:10) '.mat'],'ShipSST_interp')
%         set(gca,'FontSize',12) 
        %%
        subplot(5,1,3)% Licor timeseries
        %%
        %plot patch first underneath.
        for r=1:length(inds)
        x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
        y= [ max(Licor_fco2_surface) max(Licor_fco2_surface) 0 0];
        patch(Licor_Datetime(x),y,colour_greyshade); 
        end
        hold on
        %
        plot(Licor_Datetime,LicorxCO2Bcal_lag,'LineWidth',2); hold on;
        plot(Licor_Datetime,Licor_fco2_lag,'m','LineWidth',2); hold on;
        plot(Licor_Datetime,Licor_fco2_surface_lag,'k','LineWidth',2); 
        plot(Licor_Datetime,underwaypCO2_sw_interp,'g','LineWidth',2); 
        set(gca,'FontSize',12)
        dynamicDateTicks([], [], 'HH:MM');
        %add text labels
        textbp( 'XCO_2','color','b','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        textbp( 'F CO_2','color','m','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        textbp( 'F CO_2 underway','color','m','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        textbp( ['F CO_2' char(10) 'insitu corrected'],'color','k','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        title('CO_2 Timeseires Plot','FontSize',12);xlabel('Time','FontSize',12); ylabel(['CO_2' char(10) '(',num2str(micro_symbol),'atm)'],'FontSize',12);
        mnco2m2std=mean(Licor_fco2_surface(inds(1):indt(end)))-(2*std(Licor_fco2_surface(inds(1):indt(end))));
        mnco2p2std=mean(Licor_fco2_surface(inds(1):indt(end)))+(2*std(Licor_fco2_surface(inds(1):indt(end))));
%         set(gca, 'YLim', [mnco2m2std,mnco2p2std]);
        set(gca,'FontSize',12); 
        %plot bin start and end points and numerical labels
        plot(Licor_Datetime([indslag(1:end),indslag(1:end)]),[0,max(Licor_fco2_surface)],'Color','r','LineWidth',2);
        plot(Licor_Datetime([indtlag(1:end),indtlag(1:end)]),[0,max(Licor_fco2_surface)],'Color',colour_green,'LineWidth',2);
        if ptxtnum== 1
        text(Licor_Datetime(indslag), (ones(length(inds),1)*mnco2p2std), indtextlabel3,'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        text(Licor_Datetime(indtlag), (ones(length(inds),1)*mnco2p2std), indtextlabel3,'color',colour_green,'Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        end; hold off;
        %%
        subplot(5,1,4)%Salinity timseries
        %%
        %plot patch first underneath.
        for r=1:length(inds)
        x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
        y= [ max(CTDsalinity_interp) max(CTDsalinity_interp) 0 0];
        patch(Licor_Datetime(x),y,colour_greyshade); 
        end
        hold on
        plot(Licor_Datetime,CTDsalinity_interp,'LineWidth',2);
        dynamicDateTicks([], [], ' HH:MM');
        title('Salinity timeseires plot','FontSize',12);xlabel('Time','FontSize',12); ylabel(['Salinity' char(10) '(PSU)'],'FontSize',12); set(gca,'FontSize',12); hold on;
        %set axis limits as 2 stnd devs
        mnsalm2std=mean(CTDsalinity_interp(inds(1):indt(end)))-(2*std(CTDsalinity_interp(inds(1):indt(end))));
        mnsalp2std=mean(CTDsalinity_interp(inds(1):indt(end)))+(2*std(CTDsalinity_interp(inds(1):indt(end))));
        set(gca, 'YLim', [mnsalm2std,mnsalp2std])
        textbp( 'Salinity','color','b','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        %PLot bin divisers and labels
        plot(Licor_Datetime([indslag(1:end-1),indslag(1:end-1)]),[0,max(CTDsalinity_interp)],'Color','r','LineWidth',2);
        plot(Licor_Datetime([indtlag(2:end),indtlag(2:end)]),[0,max(CTDsalinity_interp)],'Color',colour_green,'LineWidth',2);
        if ptxtnum== 1
        text(Licor_Datetime(indslag), (ones(length(inds),1)*mnsalp2std), indtextlabel3,'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        text(Licor_Datetime(indtlag), (ones(length(inds),1)*mnsalp2std), indtextlabel3,'color',colour_green,'Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        end; hold off;
        %%
        subplot(5,1,5)%Flow rate timseries
        %%
        %plot patch first underneath.
        for r=1:length(inds)
        x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
        y= [ max(MeasuredFlow_interp) max(MeasuredFlow_interp) 0 0];
        patch(Licor_Datetime(x),y,colour_greyshade); 
        end
        hold on
        plot(Licor_Datetime,MeasuredFlow_interp,'LineWidth',2); hold on;
        plot(Licor_Datetime,smoothflow,'m','LineWidth',2);
        dynamicDateTicks([], [], ' HH:MM');
        title('Flowrate timseries plot','FontSize',12);xlabel('Time','FontSize',12); ylabel(['Flowrate' char(10) 'L/min'],'FontSize',12);set(gca,'FontSize',12);
        %set axis limits as 2 stnd devs
        mnflowm2std=mean(smoothflow(inds(1):indt(end)))-(2*std(smoothflow(inds(1):indt(end))));
        mnflowp2std=mean(smoothflow(inds(1):indt(end)))+(2*std(smoothflow(inds(1):indt(end))));
        set(gca, 'YLim', [mnflowm2std,mnflowp2std])
        textbp( 'Smoothed Flowrate','color','m','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        %PLot bin divisers and labels
        plot(Licor_Datetime([indslag(1:end-1),indslag(1:end-1)]),[0,max(CTDsalinity_interp)],'Color','r','LineWidth',2);
        plot(Licor_Datetime([indtlag(2:end),indtlag(2:end)]),[0,max(CTDsalinity_interp)],'Color',colour_green,'LineWidth',2);
        if ptxtnum== 1
        text(Licor_Datetime(indslag), (ones(length(inds),1)*mnflowp2std), indtextlabel3,'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        text(Licor_Datetime(indtlag), (ones(length(inds),1)*mnflowp2std), indtextlabel3,'color',colour_green,'Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        end
        
        
        
        
        figure(4) %timeseries data publication
        %%
        subplot(4,1,1)% depth
        %%
        %Plot depth profile and index points
        set(gcf, 'Color','w','Position', get(0,'Screensize'));
        plot(Licor_Datetime(ind(1):indt(end)),CTDdepth_interp(ind(1):indt(end)),'.','LineWidth',4);  hold on;
        %plot(Licor_Datetime(ind(1):indt(end)),smoothCTD(ind(1):indt(end)),'.r');
%         plot(Licor_Datetime(ind),smoothCTD(ind),'.g','MarkerSize',32)
        dynamicDateTicks([], [], ' HH:MM');
        set(gca,'FontSize',22);
        set(gca,'FontSize',22);
        set(gca,'YDir','reverse')
        ylim([0.15 5.3])
        ylabel('Depth (m)','Fontsize',22);
        text(0.02,1.02,'(a)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
        %%
        subplot(4,1,2) %temp
        %%
        %plot patch first underneath.
        for r=1:length(inds)
            x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
            y= [ max(Equilibrator_Temperature_interp) max(Equilibrator_Temperature_interp) 0 0];
            patch(Licor_Datetime(x),y,colour_greyshade);
        end
        %add CTD temperature to 2nd axis
        hold on;
        plot(Licor_Datetime(ind(1):indt(end)),CTDTemperature_interp(ind(1):indt(end)),'m','LineWidth',4);
        dynamicDateTicks([], [], ' HH:MM');
        ylabel({['Temperature (',num2str(degree_symbol),'C)']},'FontSize',22);
        %plot bin start and end points and numerical labels
        mmequm2std=mean(Equilibrator_Temperature_interp(inds(1):indt(end)))-(2*std(Equilibrator_Temperature_interp(inds(1):indt(end))));
        mnequp2std=mean(Equilibrator_Temperature_interp(inds(1):indt(end)))+(2*std(Equilibrator_Temperature_interp(inds(1):indt(end))));
        %Plot ctd temp
        mnCTDtempm2std=mean(CTDTemperature_interp(inds(1):indt(end)))-(3*std(CTDTemperature_interp(inds(1):indt(end))));
        mnCTDtempp2std=mean(CTDTemperature_interp(inds(1):indt(end)))+(3*std(CTDTemperature_interp(inds(1):indt(end))));
        set(gca, 'YLim', [mnCTDtempm2std,mnCTDtempp2std])
        set(gca,'FontSize',22); hold on;
        if ptxtnum== 1
            text(Licor_Datetime(indslag), (ones(length(inds),1)*mnCTDtempp2std), indtextlabel3,'color','r','FontSize',16,'Fontweight','bold','BackgroundColor','w');
            text(Licor_Datetime(indtlag), (ones(length(inds),1)*mnCTDtempp2std), indtextlabel3,'color',colour_green,'FontSize',16,'Fontweight','bold','BackgroundColor','w');
        end;
        text(0.02,1.02,'(b)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
        hold off;
        %%
        subplot(4,1,3)%salinity
        %%
        %plot patch first underneath.
        for r=1:length(inds)
            x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
            y= [ 2*max(CTDsalinity_interp) 2*max(CTDsalinity_interp) 0 0];
            patch(Licor_Datetime(x),y,colour_greyshade);
        end
        hold on
        plot(Licor_Datetime(ind(1):indt(end)),CTDsalinity_interp(ind(1):indt(end)),'LineWidth',4);
        dynamicDateTicks([], [], ' HH:MM');
        ylabel(['Salinity' char(10) '(PSU)'],'FontSize',22); set(gca,'FontSize',22); hold on;
        %set axis limits as 2 stnd devs
        mnsalm2std=mean(CTDsalinity_interp(inds(1):indt(end)))-(2*std(CTDsalinity_interp(inds(1):indt(end))));
        mnsalp2std=mean(CTDsalinity_interp(inds(1):indt(end)))+(2*std(CTDsalinity_interp(inds(1):indt(end))));
        % set(gca, 'YLim', [mnsalm2std,mnsalp2std])
        set(gca,'YLim',[ 35.4,35.5])
        text(0.02,1.02,'(c)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
        if ptxtnum== 1
            text(Licor_Datetime(indslag), (ones(length(inds),1)*mnsalp2std), indtextlabel3,'color','r','Fontsize',16,'Fontweight','bold','BackgroundColor','w');
            text(Licor_Datetime(indtlag), (ones(length(inds),1)*mnsalp2std), indtextlabel3,'color',colour_green,'Fontsize',16,'Fontweight','bold','BackgroundColor','w');
        end; hold off;
        %%
        subplot(4,1,4)% Licor timeseries
        %%
        %plot patch first underneath.
        for r=1:length(inds)
            x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
            y= [ max(Licor_fco2_surface_lag) max(Licor_fco2_surface_lag) 0 0];
            patch(Licor_Datetime(x),y,colour_greyshade);
        end
        hold on
        plot(Licor_Datetime(ind(1):indt(end)),Licor_fco2_surface_lag(ind(1):indt(end)),'k','LineWidth',4);
        dynamicDateTicks([], [], 'HH:MM');
        xlabel('Time (UTC)','FontSize',22)
        ylabel(['fCO_{2(sw)}' char(10) '(',num2str(micro_symbol),'atm)'],'FontSize',22);
        mnco2m2std=mean(Licor_fco2_surface(inds(1):indt(end)))-(3*std(Licor_fco2_surface(inds(1):indt(end))));
        mnco2p2std=mean(Licor_fco2_surface(inds(1):indt(end)))+(3*std(Licor_fco2_surface(inds(1):indt(end))));
        set(gca, 'YLim', [mnco2m2std,mnco2p2std]);
        set(gca,'FontSize',22);
        set(gca,'FontSize',22)
        text(0.02,1.02,'(d)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
        if ptxtnum== 1
            text(Licor_Datetime(indslag), (ones(length(inds),1)*mnco2p2std), indtextlabel3,'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w');
            text(Licor_Datetime(indtlag), (ones(length(inds),1)*mnco2p2std), indtextlabel3,'color',colour_green,'Fontsize',12,'Fontweight','bold','BackgroundColor','w');
        end; hold off;
        %%
        
        
        figure(5)
        set(gcf, 'Color','w','Position', get(0,'Screensize'));
        nsubplots=7;
%         subplot(nsubplots,1,1) %significant wave height
%         %%
%         %plot patch first underneath.
%         for r=1:length(inds)
%         x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
%         y= [ max(Sigwave_interp) max(Sigwave_interp) 0 0];
%         patch(Licor_Datetime(x),y,colour_greyshade); 
%         end
%         hold on
%         plot(Licor_Datetime(ind(1):indt(end)),Sigwave_interp(ind(1):indt(end)),'LineWidth',2);
%         title('Significant Wave height','FontSize',12);
%         xlabel('Time','FontSize',12); 
%         ylabel({['Wave' char(10) 'height(m)']},'FontSize',8);
%         set(gca,'FontSize',12)
%         dynamicDateTicks([], [], ' HH:MM');
%         %textbp( 'Waveheight' ,'Color','b','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
%         textbp('Bin Start','color','r','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
%         textbp('Bin End','color',colour_green,'Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
%         hold on
%         %plot bin start and end points and numerical labels
%         addaxisplot(Licor_Datetime([indslag(1:end),indslag(1:end)]),[0,max(Sigwave_interp)],1,'Color','r','LineWidth',2);
%         addaxisplot(Licor_Datetime([indtlag(1:end),indtlag(1:end)]),[0,max(Sigwave_interp)],1,'Color',colour_green,'LineWidth',2);
%         mmhsm2std=mean(Sigwave_interp(inds(1):indt(end)))-(2*std(Sigwave_interp(inds(1):indt(end))));
%         mnhsp2std=mean(Sigwave_interp(inds(1):indt(end)))+(2*std(Sigwave_interp(inds(1):indt(end))));
%         set(gca, 'YLim', [mmhsm2std,mnhsp2std]);
%         if ptxtnum== 1
%         text(Licor_Datetime(indslag), (ones(length(inds),1)*mnhsp2std), indtextlabel3,'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
%         text(Licor_Datetime(indtlag), (ones(length(inds),1)*mnhsp2std), indtextlabel3,'color',colour_green,'Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
%         end
        %%
        subplot(nsubplots,1,2)%Distance
        %%
        for r=1:length(inds)
        x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
        y= [ max(dist) max(dist) 0 0];
        patch(Licor_Datetime(x),y,colour_greyshade); 
        end
        hold on
        plot(Licor_Datetime(ind(1):indt(end)),dist(ind(1):indt(end)),'LineWidth',2);
        title('Distance from profile start point','FontSize',12);
        xlabel('Time','FontSize',12); 
        ylabel({['Distance' char(10) '(km)']},'FontSize',8);
        set(gca,'FontSize',12)
        dynamicDateTicks([], [], ' HH:MM');
        %textbp( 'Distance' ,'Color','b','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        hold on
        %plot bin start and end points and numerical labels
        addaxisplot(Licor_Datetime([indslag(1:end),indslag(1:end)]),[0,max(dist)],1,'Color','r','LineWidth',2);
        addaxisplot(Licor_Datetime([indtlag(1:end),indtlag(1:end)]),[0,max(dist)],1,'Color',colour_green,'LineWidth',2);
        mmdsm2std=mean(dist(inds(1):indt(end)))-(2*std(dist(inds(1):indt(end))));
        mndsp2std=mean(dist(inds(1):indt(end)))+(2*std(dist(inds(1):indt(end))));
        set(gca, 'YLim', [mmdsm2std,mndsp2std]);
        if ptxtnum== 1
        text(Licor_Datetime(indslag), (ones(length(inds),1)*mndsp2std), indtextlabel3,'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        text(Licor_Datetime(indtlag), (ones(length(inds),1)*mndsp2std), indtextlabel3,'color',colour_green,'Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        end
        %%
        subplot(nsubplots,1,3)%Chlorphyll a
        %%
        for r=1:length(inds)
        x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
        y= [ max(Chla_interp) max(Chla_interp) 0 0];
        patch(Licor_Datetime(x),y,colour_greyshade); 
        end
        hold on
        plot(Licor_Datetime(ind(1):indt(end)),Chla_interp(ind(1):indt(end)),'LineWidth',2);
        title('Chlorphyll-a','FontSize',12);
        xlabel('Time','FontSize',12); 
        ylabel({[' Chloprhyll-a' char(10) '(Mg M^{-3})']},'FontSize',8);
        set(gca,'FontSize',12)
        dynamicDateTicks([], [], ' HH:MM');
        %textbp( 'Flourescence' ,'Color','b','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        hold on
        %plot bin start and end points and numerical labels
        addaxisplot(Licor_Datetime([indslag(1:end),indslag(1:end)]),[0,max(Chla_interp)],1,'Color','r','LineWidth',2);
        addaxisplot(Licor_Datetime([indtlag(1:end),indtlag(1:end)]),[0,max(Chla_interp)],1,'Color',colour_green,'LineWidth',2);
        mmfluom2std=mean(Chla_interp(inds(1):indt(end)))-(2*std(Chla_interp(inds(1):indt(end))));
        mnfluop2std=mean(Chla_interp(inds(1):indt(end)))+(2*std(Chla_interp(inds(1):indt(end))));
         if ptxtnum== 1
        text(Licor_Datetime(indslag), (ones(length(inds),1)*mnfluop2std), indtextlabel3,'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        text(Licor_Datetime(indtlag), (ones(length(inds),1)*mnfluop2std), indtextlabel3,'color',colour_green,'Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        end
        %%
        subplot(nsubplots,1,4) %irradiance
        %%
        for r=1:length(inds)
        x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
        y= [ max(Parport_interp) max(Parport_interp) 0 0];
        patch(Licor_Datetime(x),y,colour_greyshade); 
        end
        hold on
        h1=plot(Licor_Datetime(ind(1):indt(end)),Parport_interp(ind(1):indt(end)),'LineWidth',2);

        title('PAR','FontSize',12);
        xlabel('Time','FontSize',12); 
        ylabel({[' Watts' char(10) 'M^{-3}']},'FontSize',8);
        set(gca,'FontSize',12)
        dynamicDateTicks([], [], ' HH:MM');
        hold on
        %plot bin start and end points and numerical labels
        addaxisplot(Licor_Datetime([indslag(1:end),indslag(1:end)]),[0,max(Parport_interp)],1,'Color','r','LineWidth',2);
        addaxisplot(Licor_Datetime([indtlag(1:end),indtlag(1:end)]),[0,max(Parport_interp)],1,'Color',colour_green,'LineWidth',2);
        mnirrp2std=mean(Parport_interp(inds(1):indt(end)))+(2*std(Parport_interp(inds(1):indt(end))));
        mnirrm2std=mean(Parport_interp(inds(1):indt(end)))-(2*std(Parport_interp(inds(1):indt(end))));
        set(gca, 'YLim', [mnirrm2std,mnirrp2std]);
        if ptxtnum== 1
        text(Licor_Datetime(indslag), (ones(length(inds),1)*mnirrp2std), indtextlabel3,'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        text(Licor_Datetime(indtlag), (ones(length(inds),1)*mnirrp2std), indtextlabel3,'color',colour_green,'Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
         end
        %%
        subplot(nsubplots,1,5) % Absolute wind speed
        %%
        for r=1:length(inds)
        x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
        y= [ max(Windspeedabs_interp) max(Windspeedabs_interp) 0 0];
        patch(Licor_Datetime(x),y,colour_greyshade); 
        end
        hold on
        plot(Licor_Datetime(ind(1):indt(end)),Windspeedabs_interp(ind(1):indt(end)),'LineWidth',2);
        title('Wind Speed','FontSize',12);
        xlabel('Time','FontSize',12); 
        ylabel({[' Wind' char(10) 'Speed {m/s}']},'FontSize',8);
        set(gca,'FontSize',12)
        dynamicDateTicks([], [], ' HH:MM');
        %textbp( 'Wind Speed' ,'Color','b','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        hold on
        %plot bin start and end points and numerical labels
        addaxisplot(Licor_Datetime([indslag(1:end),indslag(1:end)]),[0,max(Windspeedabs_interp)],1,'Color','r','LineWidth',2);
        addaxisplot(Licor_Datetime([indtlag(1:end),indtlag(1:end)]),[0,max(Windspeedabs_interp)],1,'Color',colour_green,'LineWidth',2);
        mnwindp2std=mean(Windspeedabs_interp(inds(1):indt(end)))+(2*std(Windspeedabs_interp(inds(1):indt(end))));
        mnwindm2std=mean(Windspeedabs_interp(inds(1):indt(end)))-(2*std(Windspeedabs_interp(inds(1):indt(end))));
        set(gca, 'YLim', [mnwindm2std,mnwindp2std]);
        if ptxtnum== 1
        text(Licor_Datetime(indslag), (ones(length(inds),1)*mnwindp2std), indtextlabel3,'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        text(Licor_Datetime(indtlag), (ones(length(inds),1)*mnwindp2std), indtextlabel3,'color',colour_green,'Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        end
        %%      
        subplot(nsubplots,1,6) %Air temperature
        %%
        for r=1:length(inds)
        x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
        y= [ max(AtmTemp_interp) max(AtmTemp_interp) 0 0];
        patch(Licor_Datetime(x),y,colour_greyshade); 
        end
        hold on
        plot(Licor_Datetime(ind(1):indt(end)),AtmTemp_interp(ind(1):indt(end)),'LineWidth',2);
        title('Air Temperature','FontSize',12);
        xlabel('Time','FontSize',12); 
        ylabel({['Air' char(10) 'Temperature (',num2str(degree_symbol),'C)']},'FontSize',8);         
        set(gca,'FontSize',12)
        dynamicDateTicks([], [], ' HH:MM');
        %textbp( 'Air Temperature' ,'Color','b','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        hold on
        %plot bin start and end points and numerical labels
        addaxisplot(Licor_Datetime([indslag(1:end),indslag(1:end)]),[0,max(AtmTemp_interp)],1,'Color','r','LineWidth',2);
        addaxisplot(Licor_Datetime([indtlag(1:end),indtlag(1:end)]),[0,max(AtmTemp_interp)],1,'Color',colour_green,'LineWidth',2);
        mnairtp2std=mean(AtmTemp_interp(inds(1):indt(end)))+(2*std(AtmTemp_interp(inds(1):indt(end))));
        mnairtm2std=mean(AtmTemp_interp(inds(1):indt(end)))-(2*std(AtmTemp_interp(inds(1):indt(end))));
        set(gca, 'YLim', [mnairtm2std,mnairtp2std]);
        if ptxtnum== 1
        text(Licor_Datetime(indslag), (ones(length(inds),1)*mnairtp2std), indtextlabel3,'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        text(Licor_Datetime(indtlag), (ones(length(inds),1)*mnairtp2std), indtextlabel3,'color',colour_green,'Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        end

        subplot(nsubplots,1,7) %PAR
        %%
       

        %%
        figure(6) %publication ready underway timeseries figure
        %%
        %find indexes for long underway timeseries
        val = Licor_Datetime(inds(1))- (48*3600/(24*60*60)); %value to find - change the hours to set start point of graph
        tmp = abs(Und_DT_interplongtime-val);
       [idx startindex] = min(tmp); %index of closest value
        val = Licor_Datetime(indt(end))+ (0*3600/(24*60*60));  %value to find - add time on do end of figur eto extend to after deployment
        tmp = abs(Und_DT_interplongtime-val);
        [idx endindex] = min(tmp); %index of closest value
        set(gcf, 'Color','w','Position', get(0,'Screensize'));
        nsubplots=3;
        plotthreshold=3600*6;
        %%
%         subplot(nsubplots,1,3) %significant wave height
%         %%
%         %plot patch first underneath. if statement to make the profile 1 patch if it is longer than amount set just below
%         plotthreshold=3600*6;
%         if (endindex-startindex)<(plotthreshold)
%             for r=1:length(inds)
%                 x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
%                 y= [ max(Sigwave_interplongtime) max(Sigwave_interplongtime) 0 0];
%                 patch(Licor_Datetime(x),y,colour_greyshade);
%             end
%         else
%             x=[indtlag(end) indslag(1) indslag(1) indtlag(end)];
%             y= [ max(Sigwave_interplongtime) max(Sigwave_interplongtime) 0 0];
%             patch(Licor_Datetime(x),y,colour_greyshade);
%         end
%         hold on
%         despiked_Sigwave_interplongtime=func_despike_phasespace3d( (Sigwave_interplongtime(startindex:endindex))); 
%         plot(Und_DT_interplongtime(startindex:endindex),despiked_Sigwave_interplongtime,'LineWidth',4);
%         xlabel('Time (UTC)','FontSize',22);          ylabel({['Wave' char(10) 'height(m)']},'FontSize',22);
%         set(gca,'FontSize',22)
%         dynamicDateTicks([], [], 'dd/mm');
%         mmhsm2std=mean(Sigwave_interplongtime(startindex:endindex))-(3*std(Sigwave_interplongtime(startindex:endindex)));
%         mnhsp2std=mean(Sigwave_interplongtime(startindex:endindex))+(3*std(Sigwave_interplongtime(startindex:endindex)));
%         set(gca, 'YLim', [mmhsm2std,mnhsp2std]);
%         if ptxtnum== 1
%         text(Licor_Datetime(indslag), (ones(length(inds),1)*mnhsp2std), indtextlabel3,'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
%         text(Licor_Datetime(indtlag), (ones(length(inds),1)*mnhsp2std), indtextlabel3,'color',colour_green,'Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
%         end
%         text(0.02,1.04,'(c)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
%         %%
%         subplot(nsubplots,1,1) %irradiance
%         %%
%         %if statement to make the profile 1 patch if it is longer thanamount set just below
%         if (endindex-startindex)<(plotthreshold)
%             for r=1:length(inds)
%                 x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
%                 y= [ max(Downirrstar_interplongtime) max(Downirrstar_interplongtime) 0 0];
%                 patch(Licor_Datetime(x),y,colour_greyshade);
%             end
%         else
%             x=[indtlag(end) indslag(1) indslag(1) indtlag(end)];
%             y= [ max(Downirrstar_interplongtime) max(Downirrstar_interplongtime) 0 0];
%             patch(Licor_Datetime(x),y,colour_greyshade);
%         end
%         hold on
%         h1=plot(Und_DT_interplongtime(startindex:endindex),Downirrstar_interplongtime(startindex:endindex)+32,'LineWidth',4);
%         ylabel({[' Irradiance' char(10) '(Wm^{-2})']},'FontSize',22);
%         set(gca,'FontSize',22)
%         dynamicDateTicks([], [], 'dd/mm');
%         hold on
% %         mnirrp2std=mean(Downirrstar_interplongtime(startindex:endindex))+(3*std(Downirrstar_interplongtime(startindex:endindex)));
% %         mnirrm2std=mean(Downirrstar_interplongtime(startindex:endindex))-(3*std(Downirrstar_interplongtime(startindex:endindex)));
% %         set(gca, 'YLim', [mnirrm2std,mnirrp2std]);
%       set(gca, 'YLim', [0,1250]);

        if ptxtnum== 1
        text(Licor_Datetime(indslag), (ones(length(inds),1)*mnirrp2std), indtextlabel3,'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        text(Licor_Datetime(indtlag), (ones(length(inds),1)*mnirrp2std), indtextlabel3,'color',colour_green,'Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        end
        text(0.02,1.02,'(a)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
        %%
        subplot(nsubplots,1,2) % Absolute wind speed
        %%
        %if statement to make the profile 1 patch if it is longer than amount set just below        
        if (endindex-startindex)<(plotthreshold)
            for r=1:length(inds)
                x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
                y= [ max(Windspeedabs_interplongtime) max(Windspeedabs_interplongtime) 0 0];
                patch(Licor_Datetime(x),y,colour_greyshade);
            end
        else
            x=[indtlag(end) indslag(1) indslag(1) indtlag(end)];
            y= [ max(Windspeedabs_interplongtime) max(Windspeedabs_interplongtime) 0 0];
            patch(Licor_Datetime(x),y,colour_greyshade);
        end
        hold on
        t=medfilt1 ((Windspeedabs_interplongtime(startindex:endindex)),60);
        despiked_Windspeedabs_interplongtime=func_despike_phasespace3d( (Windspeedabs_interplongtime(startindex:endindex)));
        plot(Und_DT_interplongtime(startindex:endindex),despiked_Windspeedabs_interplongtime,'LineWidth',4);
        ylabel({[' Wind' char(10) 'Speed (ms^{-1})']},'FontSize',22);
        set(gca,'FontSize',22)
        dynamicDateTicks([], [], 'dd/mm');
        hold on
        mnwindp2std=mean(Windspeedabs_interplongtime(startindex:endindex))+(3*std(Windspeedabs_interplongtime(startindex:endindex)));
        mnwindm2std=mean(Windspeedabs_interplongtime(startindex:endindex))-(3*std(Windspeedabs_interplongtime(startindex:endindex)));
        set(gca, 'YLim', [mnwindm2std,mnwindp2std]);
        if ptxtnum== 1
        text(Licor_Datetime(indslag), (ones(length(inds),1)*mnwindp2std), indtextlabel3,'color','r','Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        text(Licor_Datetime(indtlag), (ones(length(inds),1)*mnwindp2std), indtextlabel3,'color',colour_green,'Fontsize',12,'Fontweight','bold','BackgroundColor','w'); 
        end
        text(0.02,1.02,'(b)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
        %%      
        
        figure(7) % profiles of co2, temp and salinity
        set(gcf, 'Color','w','Position', get(0,'Screensize'));  
        subplot(2,2,1); %Temperature plot
        %%
        line(avgtemp,avgdepth,'LineWidth',0.5,'LineStyle','--');hold on;
        errorbar(avgtemp,avgdepth,(2*semdepth))
        herrorbar(avgtemp,avgdepth,(2*semtemp))
        if axT==1
        scatter(avgtemp,avgdepth,50,avgTstr,'filled');
        u=colorbar; cbdate(u,'HH:MM');ylabel(u,'Time since start(Hours Minutes)','FontSize',16);
        end
        if axT==0
        scatter(avgtemp,avgdepth,50,avgtime,'filled');
        u=colorbar; cbdate(u,'HH:MM');ylabel(u,'Time (Hours Minutes)','FontSize',16);
        end
        hold on;        
        xlabel(['Temperature (',num2str(degree_symbol),'C)'],'fontsize',22); ylabel('Depth (m)','fontsize',16);
        set(gca,'FontSize',16)
        set(gca,'FontSize',16)
        set(gca,'YDir','reverse'); hold off;
        %%
        subplot(2,2,2);%Salinity plot
        %%
        line(avgsal,avgdepth,'LineWidth',0.5,'LineStyle','--');hold on;
        errorbar(avgsal,avgdepth,(2*semdepth))
        herrorbar(avgsal,avgdepth,(2*semsal))
        if axT==1
        scatter(avgsal,avgdepth,50,avgTstr,'filled');hold on;
        z=colorbar; cbdate(z,'HH:MM');ylabel(z,'Time since start (Hours Minutes)','FontSize',16);
        end
        if axT==0
        scatter(avgsal,avgdepth,50,avgtime,'filled');hold on;
        z=colorbar; cbdate(z,'HH:MM');ylabel(z,'Time(Hours Minutes)','FontSize',16);
        end
        xlabel('Salinity (PSU)','fontsize',22); ylabel('Depth (m)','fontsize',16);
        set(gca,'FontSize',16)
        set(gca,'FontSize',16)
        set(gca,'YDir','reverse');hold off;
        %%
%         subplot(2,2,3);%FCO2 plot
%         %%
%         line(avgfco2,avgdepth,'LineWidth',0.5,'LineStyle','--');hold on;
%         errorbar(avgfco2,avgdepth,(2*semdepth))
%         herrorbar(avgfco2,avgdepth,(2*semfco2))
%         if axT==1
%         scatter(avgfco2,avgdepth,50,avgTstr,'filled');hold on;
%         y=colorbar; cbdate(y,'HH:MM');ylabel(y,'Time since start(Hours Minutes)','FontSize',16);
%         end
%         if axT==0
%         scatter(avgfco2,avgdepth,50,avgtime,'filled');hold on;
%         y=colorbar; cbdate(y,'HH:MM');ylabel(y,'Time(Hours Minutes)','FontSize',16);
%         end
%         xlabel(['fCO_2 (',num2str(micro_symbol),'atm)'],'fontsize',22); ylabel('Depth (m)','fontsize',16);
%         set(gca,'YDir','reverse');
%         set(gca,'FontSize',16);
%         set(gca,'FontSize',16);hold off;
%         %%
        subplot(2,2,4); %fco2 temp corrected
        %%
        line(avgfco2surf,avgdepth,'LineWidth',0.5,'LineStyle','--');hold on;
        errorbar(avgfco2surf,avgdepth,(2*semdepth))
        herrorbar(avgfco2surf,avgdepth,2*semavgfco2surf)
        if axT==1
        scatter(avgfco2surf,avgdepth,50,avgTstr,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');ylabel(v,'Time since start (Hours Minutes)','FontSize',16);
        end
        if axT==0
        scatter(avgfco2surf,avgdepth,50,avgtime,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');ylabel(v,'Time(Hours Minutes)','FontSize',16);
        end
        xlabel(['fCO_2 (',num2str(micro_symbol),'atm)'],'fontsize',22); ylabel('Depth (m)','fontsize',16);
        set(gca,'YDir','reverse'); 
        set(gca,'FontSize',16)
        set(gca,'FontSize',16)
        
        p=mtit(['Profiles for' ' '  fList(j,1:10) ],'fontsize',24,'xoff',-.025,'yoff',.025);

        if ~exist('Figures/Profiles', 'dir')
                 mkdir('Figures','Profiles');
        end
        saveas(gcf,[pwd '/Figures/Profiles/' fList(j,1:10) '.png']); hold off;
        %%
 
        %%
        figure(8) % near surface data
        %%
      plot(Licor_Datetime,T_03m_interp,'Color',[ 0+(0.3/7) , 0, 1-(0.3/7)]);hold on;
      plot(Licor_Datetime,T_06m_interp,'Color',[ 0+(0.6/7) ,0, 1-(0.6/7)]);
      plot(Licor_Datetime,T_15m_interp,'Color',[ 0+(1.5/7) ,0, 1-(1.5/7)]);
      plot(Licor_Datetime,T_35m_interp,'Color',[ 0+(3.5/7) ,0, 1-(3.5/7)]);
      plot(Und_DT_interp,SeaTemp_interp,'Color',[ 0+(5.5/7) ,0, 1-(5.5/7)]);
      colormatrix=[ 0+(avgdepth/7) ;zeros(1,length(avgdepth)); 1-(avgdepth/7)];
      scatter(avgtime',avgtemp',18,colormatrix','filled');
       dynamicDateTicks([], [], ' HH:MM');
      ylabel(['Temperature (',num2str(degree_symbol),'C)'],'FontSize',16);
       xlabel('Time','FontSize',16);
       set(gca,'FontSize',16);
       set(gca,'FontSize',16);
       %create custom colormap red to blue
       for z=1:1000;
       mymap(z,:)=[(0.001*z) (z*0) (1-(0.001*z))];
       end
       colormap(mymap)
       h = colorbar;
%        cblabel(h,'Depth(m)','FontSize',16)
       caxis([0,7])
       set(gca,'FontSize',16);
       set( h, 'YDir', 'reverse' );       % print distance from ship to CCS mooring on figure
       dist2mooring= pos2dist(49.4042,-8.5111,Latitude_interp(1),Longitude_interp(1),2);
%       textbp(['Distance from ship =' ,num2str(dist2mooring),'km'],'color','b','Fontsize',10,'Fontweight','bold','BackgroundColor','w');
%        jkl=legend('0.3m','0.6m','1.5m','3.5m','7.0m','5.5m Ship SST','NSOP SST')
%        set(jkl,'Location','NorthWest') 
        %%
        
       figure(9)% salinity and temp subplots with mooring data
       set(gcf, 'Color','w','Position', get(0,'Screensize'));
       %%
       subplot(3,4,[1 2 5 6]);%Salinity plot
        %%
        colormap default;
        line(avgsal,avgdepth,'LineWidth',0.5,'LineStyle','--');hold on;
        errorbar(avgsal,avgdepth,(2*semdepth))
        herrorbar(avgsal,avgdepth,(2*semsal))
       if axT==1
           p=scatter(avgsal,avgdepth,50,avgTstr,'filled');
           %cbfreeze;
           u=colorbar; cbdate(u,'HH:MM');
           ylabel(u,'Time since start(Hours Minutes)','FontSize',16);      
       end
       if axT==0
           p=scatter(avgsal,avgdepth,50,avgtime,'filled');
           u=colorbar; cbdate(u,'HH:MM'); 
       end
        xlabel('Salinity (PSU)','fontsize',16); ylabel('Depth (m)','fontsize',16);
        set(gca,'FontSize',16)
        set(gca,'FontSize',16)
        set(gca,'Xlim',[35.4 35.5]);
        set(gca,'YDir','reverse');
        text(0.02,0.95,'(a)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
        hold off;
%         cbfreeze;
%         freezeColors(p);
%         cblabel('Time','FontSize',16,);
       %%
       subplot(3,4,[3 4 7 8]); %Temperature plot
       %%
       line(avgtemp,avgdepth,'LineWidth',0.5,'LineStyle','--');hold on;
       errorbar(avgtemp,avgdepth,(2*semdepth))
       herrorbar(avgtemp,avgdepth,(2*semtemp))
       if axT==1
           p=scatter(avgtemp,avgdepth,50,avgTstr,'filled');
           u=colorbar; cbdate(u,'HH:MM');ylabel(u,'Time since start(Hours Minutes)','FontSize',16);
       end
       if axT==0
           p=scatter(avgtemp,avgdepth,50,avgtime,'filled');
           u=colorbar; cbdate(u,'HH:MM');
       end
       hold on;
%        freezeColors(p);
%        herrorbar(avgTdep03m,depth03m,(2*semTdep03m))
       q=scatter(avgT_03m,depth03m',80,avgtime,'v');
%        freezeColors(q);
%        herrorbar(avgTdep06m,depth06m,(2*semTdep06m))
       r=scatter(avgT_06m,depth06m',80,avgtime,'v');
%        freezeColors(r);
%        herrorbar(avgTdep15m,depth15m,(2*semTdep15m))
       s=scatter(avgT_15m,depth15m',80,avgtime,'v');
%        freezeColors(s);    
%        herrorbar(avgTdep35m,depth35m,(2*semTdep35m))
       t=scatter(avgT_35m,depth35m',80,avgtime,'v');
%        freezeColors(t);
%        herrorbar(avgTdep70m,depth70m,(2*semTdep70m))
%        u=scatter(avgTdep70m,depth70m,50,avgtime,'filled','v');freezeColors(u);
       xlabel(['Temperature (',num2str(degree_symbol),'C)'],'fontsize',16); ylabel('Depth (m)','fontsize',16);
       set(gca,'FontSize',16)
       set(gca,'FontSize',16)
       set(gca,'YDir','reverse'); 
       hold off;
%        cbfreeze;
%        cblabel('Time','FontSize',16);
       text(0.02,0.95,'(b)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
       %%
       subplot(3,4,9:12); %Temperature plot mooring
       %%
    %mooring data
    C=load('C:\Users\rps207\Documents\MATLAB\Temperature mooring output analysis\Data\NSTempmooringl4.mat');
    load('C:\Users\rps207\Documents\MATLAB\Temperature mooring output analysis\Data\NSTempmooringl4.mat','T_DT');
    names=({'T_03m','T_06m','T_15m','T_24m','T_35m','T_star_29m'})';
    %create new names for variables with _interp suffix
    
        %create new names for variables with _interp suffix
    for t=1:numel(names)
        rty=[names(t) '_interp'];x=horzcat(rty{:}) ;new_names(:,t)=cellstr(x);
    end
    % Interpolate variables in a loop using new variable names
    for i=1:numel(names)
        interpvars.(new_names{i})=interp1(T_DT,C.(names{i}),Licor_Datetime);
    end
       v2struct(interpvars)
      plot(Licor_Datetime,T_03m_interp,'Color',[1-(0.3/7)  , 0, 0+(0.3/7)]);hold on
      plot(Licor_Datetime,T_06m_interp,'Color',[1-(0.6/7)  ,0, 0+(0.6/7)])
      plot(Licor_Datetime,T_15m_interp,'Color',[1-(1.5/7)  ,0, 0+(1.5/7)])
      plot(Licor_Datetime,T_35m_interp,'Color',[1-(3.5/7)  ,0, 0+(3.5/7)])
      plot(Und_DT_interp,SeaTemp_interp,'--','Color',[ 1-(5.5/7) ,0,0+(5.5/7) ])
      colormatrix=[ 1-(avgdepth/7) ;zeros(1,length(avgdepth)); 0+(avgdepth/7)];
      scatter(avgtime',avgtemp',18,colormatrix','filled')
       dynamicDateTicks([], [], ' HH:MM');
%       jkl=legend('0.3m','0.6m','1.5m','3.5m','7.0m','5.5m Ship SST','NSOP SST')
%       set(jkl,'Location','NorthWest') 
      ylabel(['Temperature (',num2str(degree_symbol),'C)'],'FontSize',16)
       xlabel('Time','FontSize',16)
       set(gca,'FontSize',16);
       set(gca,'FontSize',16);
       %create custom colormap red to blue
       for z=1:1000;
       mymap(z,:)=[(1-(0.001*z)) (z*0) (0.001*z)];
       end
       colormap(mymap)
       h = colorbar;
%        cblabel(h,'Depth(m)','FontSize',16)
       caxis([0,7])
       set(gca,'FontSize',16);
       set( h, 'YDir', 'reverse' );       % print distance from ship to CCS mooring on figure
%         depthlab= ({'0.3m' '0.6m' '1.5m' '3.5m' '7m' 'Und'}); lablist=[Tdep03m_interp,Tdep06m_interp,Tdep15m_interp,Tdep35m_interp,Tdep70m_interp,SeaTemp_interp];
%         for tyu=1:length(depthlab);
%         textfit(Licor_Datetime(end), lablist(end,tyu), depthlab(tyu),'color','k','Fontsize',8,'Fontweight','bold','BackgroundColor','w'); 
%         end
        dist2mooring= pos2dist(49.40356,-8.606,Latitude_interp(inds(3)),Longitude_interp(inds(3)),2);
        %calculate distance form starting point of profile
        Lat1=ones(length(Latitude_interp),1).*49.40356;
        Long1=ones(length(Longitude_interp),1).*-8.606; dist2mooringall=[];
        for t=1:length(Longitude_interp)
            dist2mooringall(t) = pos2dist(Lat1(t),Long1(t),Latitude_interp(t),Longitude_interp(t),2);
        end
        text(0.02,0.90,'(c)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
        textbp(['Distance from ship =' ,num2str(dist2mooring),'km'],'color','k','Fontsize',16,'Fontweight','bold','BackgroundColor','w');
        %%
      

                       %%
       figure(10)
       %%
       set(gcf, 'Color','w','Position', get(0,'Screensize'));
       subplot(1,2,1);%density plot
        %%
        line(avgden,avgdepth,'LineWidth',0.5,'LineStyle','--');hold on;
        errorbar(avgden,avgdepth,(2*semdepth))
        herrorbar(avgden,avgdepth,(2*semden))
       if axT==1
           scatter(avgden,avgdepth,80,avgTstr,'filled');
           u=colorbar; cbdate(u,'HH:MM');ylabel(u,'Time since start(Hours Minutes)','FontSize',16);
       end
       if axT==0
           scatter(avgdensity,avgdepth,80,avgtime,'filled');
           u=colorbar; cbdate(u,'HH:MM');title(u,'          Time(UTC)','FontSize',22)
       end
        xlabel('Density (kg m^{-3})','fontsize',22); ylabel('Depth (m)','fontsize',22);
        NumTicks = 3;
        L = get(gca,'XLim');
        set(gca,'XTick',linspace(L(1),L(2),NumTicks))
        set(gca,'FontSize',22)
        set(gca,'FontSize',22)
        text(0.02,1.02,'(a)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized'); 
        set(gca,'YDir','reverse');hold off;
        %%
        subplot(1,2,2); %fco2 temp corrected
        %%
        line(avgfco2surf,avgdepth,'LineWidth',0.5,'LineStyle','--');hold on;
        errorbar(avgfco2surf,avgdepth,(2*semdepth))
        herrorbar(avgfco2surf,avgdepth,2*semavgfco2surf)
        if axT==1
        scatter(avgfco2surf,avgdepth,80,avgTstr,'filled');hold on;
        u=colorbar; cbdate(u,'HH:MM');title(u,'Time since start(Hours Minutes)','FontSize',22);
        end
        if axT==0
        scatter(avgfco2surf,avgdepth,80,avgtime,'filled');hold on;
        u=colorbar; cbdate(u,'HH:MM');title(u,'          Time(UTC)','FontSize',22)
        end
        xlabel(['fCO_2(sw) (',num2str(micro_symbol),'atm)'],'fontsize',22); ylabel('Depth (m)','fontsize',22);
        set(gca,'YDir','reverse'); 
        set(gca,'FontSize',22)
        set(gca,'FontSize',22)
        text(0.02,1.02,'(b)','color','k','Fontsize',22,'Fontweight','bold','BackgroundColor','w','units','normalized'); 

      %%
      
      figure(11)%location on map
      %%
      load coast.dat
      addpath('c:\Users\rps207\Documents\Matlab\Functions\m_map');
      m_proj('mercator','lon',[-10 -4],'lat',[48 52]);hold on
      clf
      m_grid('box','fancy','tickdir','out');
      m_line(coast(:,1),coast(:,2));
      set(gca,'fontsize',14)
      ylabel('Latitude');
      set(gca,'fontsize',14)
      xlabel('Longitude');
      set(gca,'fontsize',14)
      title('DY030 Cruise Track','fontsize',14);
      set(gca,'fontsize',14)
      hold on
      %Add cruise locations during profile to map
      m_scatter(Longitude_interp(inds(1):indt(end)),Latitude_interp(inds(1):indt(end)),[],Licor_Datetime(inds(1):indt(end)))
      hcb = colorbar('horiz');
      set(get(hcb,'Xlabel'),'String','DOY')
      %add Candyfloss to map
      [C,D]=m_ll2xy(-8.5111,49.4042);
      line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','red');
      text(C,D,'Candyfloss','HorizontalAlignment', 'left','VerticalAlignment', 'bottom','fontsize',14);
      %%
      
      figure(12)
      %%
      scatter(avgflow,avgfco2surf)
      xlabel('flow')
      ylabel('fco2')
      saveas(gca,[pwd '/Figures/Flowco2/' fList(j,1:10) '.png']);
     %% figure (13) TA/DIC
     %%
     if exist(fullpathTADIC)
         figure(13)
         %%
         line(avgfco2surf,avgdepth,'LineWidth',0.5,'LineStyle','--');hold on;
         errorbar(avgfco2surf,avgdepth,(2*semdepth))
         herrorbar(avgfco2surf,avgdepth,2*semavgfco2surf)
         if axT==1
             scatter(avgfco2surf,avgdepth,50,avgTstr,'filled');hold on;
             v=colorbar; cbdate(v,'HH:MM');ylabel(v,'Time since start (Hours Minutes)','FontSize',16);
         end
         if axT==0
             scatter(avgfco2surf,avgdepth,50,avgtime,'filled');hold on;
             v=colorbar; cbdate(v,'HH:MM');ylabel(v,'Time(Hours Minutes)','FontSize',16);
         end
         %plot ta/dic derived fco2
         l=scatter(samplefco2,sampledepth,50,meanfilltime,'filled','d');hold on;

         load('C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\underwaypCO2.mat','underway_DT','underwayfCO2_sw');
         [tyu , opl]=find(underway_DT>Licor_Datetime(inds(1))&underway_DT<Licor_Datetime(indt(end)));
         intakedepth=ones(length(tyu),1)*5.5;
         scatter(underwayfCO2_sw(tyu),intakedepth,50,underway_DT(tyu))
         
         xlabel(['FCO_2(sw) (',num2str(micro_symbol),'atm)'],'fontsize',22); ylabel('Depth (m)','fontsize',16);
         set(gca,'YDir','reverse');
         set(gca,'FontSize',16)
         set(gca,'FontSize',16)
         
%          TAprec=1.5017*(ones(6,1));
         figure(14)
%          scatter(TA,sampledepth,120,meanfilltime,'filled','d');hold on;
         scatter(TA,sampledepth,120,'k','filled','d');hold on;
         errorbar(TA,sampledepth,(2*semsampledepth),'.')
%          herrorbar(TA,sampledepth,TAprec,'.')
         
         xlabel(['Total alkalinity (',num2str(micro_symbol),' mol Kg^{-1})'],'fontsize',22)
         ylabel('Depth (m)','fontsize',22);
         set(gca,'YDir','reverse');
         set(gca,'FontSize',22)
         set(gca,'FontSize',22)
         xlim([2290 2310])
     else
         sprintf('NO TA/DIC file')
     end
     
     
     
     figure(15)
     plot(Licor_Datetime,CTDsalinity_interp)
     hold on
     plot(Licor_Datetime,Salinity_interp,'g')
     
     figure(16)
     plot(Licor_Datetime,CTDTemperature_interp)
     hold on
     plot(Licor_Datetime,SeaTemp_interp,'g')
     dynamicDateTicks([], [], ' HH:MM');
     Tempcorr=avgCTDTemperature_interp_lag - avgEquilibrator_Temperature_interp;
          plot(Licor_Datetime,underwayT_sea_interp,'r')
                plot(Licor_Datetime,Equilibrator_Temperature_interp,'y')
    
          
          

     figure(17)
     plot(Tempcorr,avgfco2,'b*')
     hold on
     [p , s] = polyfit(Tempcorr, avgfco2,1);
     x1 = linspace(-2,2);
     f1 = polyval(p,x1);
     plot(x1,f1,'b');
     tempxaxis=(-2:0.01:2);
     yaxis=(min(avgfco2)*(-0.0423*tempxaxis))+min(avgfco2);
     plot(tempxaxis,yaxis,'r');
     xlabel('Temperature');
     ylabel('fco2');
     legend('Profile data','Profile fit','Takahashi 0.0423');
     xlim([min(Tempcorr) max(Tempcorr)]);
     pflTCO2=p(1)/p(2);

     
     
     
     figure(18)
      m_proj('mercator','lon',[-4.35 -4],'lat',[50.2 50.4]);
%     clf
    m_grid('linestyle','none','tickdir','out','fontsize',22);hold on;
    m_patch(coast(:,1),coast(:,2),'b','EdgeColor','none'); hold on;
%     p=m_contour(u,i,bathymetrycelticsea',[-80:5:-15],'ShowText','on');

    m_scatter(Longitude_interp*-1,Latitude_interp,10,Licor_Datetime_doy)
    v=colorbar; cbdate(v,'HH:MM');ylabel(v,'Time since start (Hours Minutes)','FontSize',22);
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
    ring1=m_range_ring(-4.217,50.25,0.5);
    set(ring1,'color','r')
    ring2=m_range_ring(-4.217,50.25,1);
    set(ring2,'color','b')
    ring3=m_range_ring(-4.217,50.25,2);
    set(ring3,'color','g')
    
    figure(19)
      m_proj('mercator','lon',[-4.35 -4],'lat',[50.2 50.4]);
%     clf
    m_grid('linestyle','none','tickdir','out','fontsize',22);hold on;
    m_patch(coast(:,1),coast(:,2),'b','EdgeColor','none'); hold on;
%     p=m_contour(u,i,bathymetrycelticsea',[-80:5:-15],'ShowText','on');

    m_scatter(Longitude_interp(inds(1):indt(end))*-1,Latitude_interp(inds(1):indt(end)),10,Licor_Datetime_doy(inds(1):indt(end)))
    v=colorbar; cbdate(v,'HH:MM');ylabel(v,'Time since start (Hours Minutes)','FontSize',22);
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
    ring1=m_range_ring(-4.217,50.25,0.5);
    set(ring1,'color','r')
    ring2=m_range_ring(-4.217,50.25,1);
    set(ring2,'color','b')
    ring3=m_range_ring(-4.217,50.25,2);
    set(ring3,'color','g')
           
     figure(20)
        for r=1:length(inds)
        x=[indtlag(r) indslag(r) indslag(r) indtlag(r)];
        y= [ max(Licor_fco2_surface) max(Licor_fco2_surface) 0 0];
        patch(Licor_Datetime(x),y,colour_greyshade); 
        end
        hold on
        %
        plot(Licor_Datetime(ind(1):indt(end)),Licor_fco2_surface_lag(ind(1):indt(end)),'k','LineWidth',2); 
        plot(Licor_Datetime(ind(1):indt(end)),underwayfCO2_sw_interp(ind(1):indt(end)),'g','LineWidth',2); 
        plot(underway_DT,underwayfCO2_sw,'k*'); 
        set(gca,'FontSize',12)
        dynamicDateTicks([], [], 'HH:MM');
        %add text labels
        textbp( 'XCO_2','color','b','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        textbp( 'F CO_2','color','m','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        textbp( 'F CO_2 underway','color','m','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        textbp( ['F CO_2' char(10) 'insitu corrected'],'color','k','Fontsize',10,'Fontweight','bold','BackgroundColor','w'); 
        title('CO_2 Timeseires Plot','FontSize',12);xlabel('Time','FontSize',12); ylabel(['CO_2' char(10) '(',num2str(micro_symbol),'atm)'],'FontSize',12);
        set(gca,'FontSize',12); 
        xlim([Licor_Datetime(1) Licor_Datetime(end)])
    end
       
        
    % subplot figures of all deployments  
    
    load coast.dat %high resollution coastline for south west of uk
    
    load('Data\bathymetry.mat')
    u=longitudeceltic(longcelticsea);
    i=(latitudeceltic(latcelticsea))';
    
    load('C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\Plymouth_tides.mat')

    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%% co2 subplot %%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Subplotcolumns=4;
    figure(101)
    subplot(4,Subplotcolumns,j)
    line(avgfco2surf,avgdepth,'LineWidth',0.5,'LineStyle','--');hold on;
    errorbar(avgfco2surf,avgdepth,(2*semdepth))
    herrorbar(avgfco2surf,avgdepth,2*semavgfco2surf)
    if axT==1
        scatter(avgfco2surf,avgdepth,80,avgTstr,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');title(v,'Time(UTC)','FontSize',16);
    end
    if axT==0
        scatter(avgfco2surf,avgdepth,80,avgtime,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');title(v,'Time(UTC)','FontSize',16);
    end
    xlabel(['fCO_2(sw) (',num2str(micro_symbol),'atm) '],'fontsize',16); ylabel('Depth (m)','fontsize',16);
    set(gca,'YDir','reverse');
    set(gca,'FontSize',16)
    set(gca,'FontSize',16)
    title({[fList(j,9:10) '/' fList(j,6:7)]});
    avgfco2surfrange=max(avgfco2surf)-min(avgfco2surf);
    avgfco2surflowerlim=min(avgfco2surf)- (avgfco2surfrange*0.1);
    avgfco2surfupperlim=max(avgfco2surf)+ (avgfco2surfrange*0.1);
%     xlim([avgfco2surflowerlim avgfco2surfupperlim])
%     NumTicks=4;
%     L=get(gca,'XLim');
%     set(gca,'XTick',linspace(L(1),L(2),NumTicks))
%     set(gca,'XTickLabel',linspace(L(1),L(2),NumTicks))
    rotateXLabels(gca(),90)

    
    if j <9;
    tyu=j;
    figure(102)
    subplot(2,Subplotcolumns,tyu)
    line(avgfco2surf,avgdepth,'LineWidth',0.5,'LineStyle','--');hold on;
    errorbar(avgfco2surf,avgdepth,(2*semdepth))
    herrorbar(avgfco2surf,avgdepth,2*semavgfco2surf)
    if axT==1
        scatter(avgfco2surf,avgdepth,80,avgTstr,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');title(v,'Time(UTC)','FontSize',16);
    end
    if axT==0
        scatter(avgfco2surf,avgdepth,80,avgtime,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');title(v,'Time(UTC)','FontSize',16);
    end
    xlabel(['fCO_2_(_s_w_) (',num2str(micro_symbol),'atm) '],'fontsize',16); ylabel('Depth (m)','fontsize',16);
    set(gca,'YDir','reverse');
    set(gca,'FontSize',16)
    set(gca,'FontSize',16)
    title({[fList(j,9:10) '/' fList(j,6:7)]});
    avgfco2surfrange=max(avgfco2surf)-min(avgfco2surf);
    avgfco2surflowerlim=min(avgfco2surf)- (avgfco2surfrange*0.1);
    avgfco2surfupperlim=max(avgfco2surf)+ (avgfco2surfrange*0.1);
%     xlim([avgfco2surflowerlim avgfco2surfupperlim])
%     NumTicks=4;
%     L=get(gca,'XLim');
%     set(gca,'XTick',linspace(L(1),L(2),NumTicks))
%     set(gca,'XTickLabel',linspace(L(1),L(2),NumTicks))
    rotateXLabels(gca(),90)
    end
    
    if j > 8;
    tyu=j-8;
    figure(103)
    subplot(2,Subplotcolumns,tyu)
    line(avgfco2surf,avgdepth,'LineWidth',0.5,'LineStyle','--');hold on;
    errorbar(avgfco2surf,avgdepth,(2*semdepth))
    herrorbar(avgfco2surf,avgdepth,2*semavgfco2surf)
    if axT==1
        scatter(avgfco2surf,avgdepth,80,avgTstr,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');title(v,'Time(UTC)','FontSize',16);
    end
    if axT==0
        scatter(avgfco2surf,avgdepth,80,avgtime,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');title(v,'Time(UTC)','FontSize',16);
    end
    xlabel(['fCO_2_(_s_w_) (',num2str(micro_symbol),'atm) '],'fontsize',16); ylabel('Depth (m)','fontsize',16);
    set(gca,'YDir','reverse');
    set(gca,'FontSize',16)
    set(gca,'FontSize',16)
    title({[fList(j,9:10) '/' fList(j,6:7)]});
    avgfco2surfrange=max(avgfco2surf)-min(avgfco2surf);
    avgfco2surflowerlim=min(avgfco2surf)- (avgfco2surfrange*0.1);
    avgfco2surfupperlim=max(avgfco2surf)+ (avgfco2surfrange*0.1);
%     xlim([avgfco2surflowerlim avgfco2surfupperlim])
%     NumTicks=4;
%     L=get(gca,'XLim');
%     set(gca,'XTick',linspace(L(1),L(2),NumTicks))
%     set(gca,'XTickLabel',linspace(L(1),L(2),NumTicks))
    rotateXLabels(gca(),90)
    end

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%% temp subplot %%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    figure(104)
    subplot(4,Subplotcolumns,j)
    line(avgtemp,avgdepth,'LineWidth',0.5,'LineStyle','--');hold on;
    errorbar(avgtemp,avgdepth,(2*semdepth))
    herrorbar(avgtemp,avgdepth,(2*semtemp))
    if axT==1
        scatter(avgtemp,avgdepth,80,avgTstr,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');title(v,'Time(UTC)','FontSize',16);
    end
    if axT==0
        scatter(avgtemp,avgdepth,80,avgtime,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');title(v,'Time(UTC)','FontSize',16);
    end
    xlabel({['Temperature (',num2str(degree_symbol),'C)']},'fontsize',16); ylabel('Depth (m)','fontsize',16);
    set(gca,'YDir','reverse');
    set(gca,'FontSize',16)
    set(gca,'FontSize',16)
    title([fList(j,1:10)]);
    avgtemprange=max(avgtemp)-min(avgtemp);
    avgtemplowerlim=min(avgtemp)- (avgtemprange*0.1);
    avgtempupperlim=max(avgtemp)+ (avgtemprange*0.1);
%     xlim([avgtemplowerlim avgtempupperlim])
%     NumTicks=4;
%     L=get(gca,'XLim');
%     set(gca,'XTick',linspace(L(1),L(2),NumTicks))
%     set(gca,'XTickLabel',linspace(L(1),L(2),NumTicks))
    rotateXLabels(gca(),90)
    
    if j < 9;
    tyu=j;
    figure(105)
    subplot(2,Subplotcolumns,tyu)
    line(avgtemp,avgdepth,'LineWidth',0.5,'LineStyle','--');hold on;
    errorbar(avgtemp,avgdepth,(2*semdepth))
    herrorbar(avgtemp,avgdepth,(2*semtemp))
    if axT==1
        scatter(avgtemp,avgdepth,80,avgTstr,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');title(v,'Time(UTC)','FontSize',16);
    end
    if axT==0
        scatter(avgtemp,avgdepth,80,avgtime,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');title(v,'Time(UTC)','FontSize',16);
    end
    xlabel({['Temperature (',num2str(degree_symbol),'C)']},'fontsize',16); ylabel('Depth (m)','fontsize',16);
    set(gca,'YDir','reverse');
    set(gca,'FontSize',16)
    set(gca,'FontSize',16)
    title({[fList(j,9:10) '/' fList(j,6:7)]});
    avgtemprange=max(avgtemp)-min(avgtemp);
    avgtemplowerlim=min(avgtemp)- (avgtemprange*0.1);
    avgtempupperlim=max(avgtemp)+ (avgtemprange*0.1);
%     xlim([avgtemplowerlim avgtempupperlim])
%     NumTicks=4;
%     L=get(gca,'XLim');
%     set(gca,'XTick',linspace(L(1),L(2),NumTicks))
%     set(gca,'XTickLabel',linspace(L(1),L(2),NumTicks))
    rotateXLabels(gca(),90)
    end
    
    if j > 8;
    tyu=j-8;
    figure(106)
    subplot(2,Subplotcolumns,tyu)
    line(avgtemp,avgdepth,'LineWidth',0.5,'LineStyle','--');hold on;
    errorbar(avgtemp,avgdepth,(2*semdepth))
    herrorbar(avgtemp,avgdepth,(2*semtemp))
    if axT==1
        scatter(avgtemp,avgdepth,80,avgTstr,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');title(v,'Time(UTC)','FontSize',16);
    end
    if axT==0
        scatter(avgtemp,avgdepth,80,avgtime,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');title(v,'Time(UTC)','FontSize',16);
    end
    xlabel({['Temperature (',num2str(degree_symbol),'C)']},'fontsize',16); ylabel('Depth (m)','fontsize',16);
    set(gca,'YDir','reverse');
    set(gca,'FontSize',16)
    set(gca,'FontSize',16)
    title({[fList(j,9:10) '/' fList(j,6:7)]});
    avgtemprange=max(avgtemp)-min(avgtemp);
    avgtemplowerlim=min(avgtemp)- (avgtemprange*0.1);
    avgtempupperlim=max(avgtemp)+ (avgtemprange*0.1);
%     xlim([avgtemplowerlim avgtempupperlim])
%     NumTicks=4;
%     L=get(gca,'XLim');
%     set(gca,'XTick',linspace(L(1),L(2),NumTicks))
%     set(gca,'XTickLabel',linspace(L(1),L(2),NumTicks))
    rotateXLabels(gca(),90)
    end
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%% salinity subplot %%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(107)
    subplot(4,Subplotcolumns,j);
    line(avgsal,avgdepth,'LineWidth',0.5,'LineStyle','--');hold on;
    errorbar(avgsal,avgdepth,(2*semdepth));
    herrorbar(avgsal,avgdepth,(2*semsal));
    if axT==1
        scatter(avgsal,avgdepth,80,avgTstr,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');title(v,'Time(UTC)','FontSize',16);
    end
    if axT==0;
        scatter(avgsal,avgdepth,80,avgtime,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');title(v,'Time(UTC)','FontSize',16);
    end
    xlabel(['Salinity(PSU)'],'fontsize',16); ylabel('Depth (m)','fontsize',16);
    set(gca,'YDir','reverse');
    set(gca,'FontSize',16);
    set(gca,'FontSize',16);
    title([fList(j,1:10)]);
    avgsalrange=max(avgsal)-min(avgsal);
    avgsallowerlim=min(avgsal)- (avgsalrange*0.1);
    avgsalupperlim=max(avgsal)+ (avgsalrange*0.1);
%     xlim([avgsallowerlim avgsalupperlim]);
%     NumTicks=4;
%     L=get(gca,'XLim');
%     set(gca,'XTick',linspace(L(1),L(2),NumTicks));
%     set(gca,'XTickLabel',linspace(L(1),L(2),NumTicks));
    rotateXLabels(gca(),90);
     
    if j < 9;
    tyu=j;
    figure(108)
    subplot(2,Subplotcolumns,tyu)
    line(avgsal,avgdepth,'LineWidth',0.5,'LineStyle','--');hold on;
    errorbar(avgsal,avgdepth,(2*semdepth));
    herrorbar(avgsal,avgdepth,(2*semsal));
    if axT==1
        scatter(avgsal,avgdepth,80,avgTstr,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');title(v,'Time(UTC)','FontSize',16);
    end
    if axT==0;
        scatter(avgsal,avgdepth,80,avgtime,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');title(v,'Time(UTC)','FontSize',16);
    end
    xlabel(['Salinity(PSU)'],'fontsize',16); ylabel('Depth (m)','fontsize',16);
    set(gca,'YDir','reverse');
    set(gca,'FontSize',16);
    set(gca,'FontSize',16);
    title({[fList(j,9:10) '/' fList(j,6:7)]});
    avgsalrange=max(avgsal)-min(avgsal);
    avgsallowerlim=min(avgsal)- (avgsalrange*0.1);
    avgsalupperlim=max(avgsal)+ (avgsalrange*0.1);
%     xlim([avgsallowerlim avgsalupperlim]);
%     NumTicks=4;
%     L=get(gca,'XLim');
%     set(gca,'XTick',linspace(L(1),L(2),NumTicks));
%     set(gca,'XTickLabel',linspace(L(1),L(2),NumTicks));
    rotateXLabels(gca(),90);
    end
    
    if j > 8;
    tyu=j-8;
    figure(109)
    subplot(2,Subplotcolumns,tyu)
    subplot(2,Subplotcolumns,tyu)
    line(avgsal,avgdepth,'LineWidth',0.5,'LineStyle','--');hold on;
    errorbar(avgsal,avgdepth,(2*semdepth));
    herrorbar(avgsal,avgdepth,(2*semsal));
    if axT==1
        scatter(avgsal,avgdepth,80,avgTstr,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');title(v,'Time(UTC)','FontSize',16);
    end
    if axT==0;
        scatter(avgsal,avgdepth,80,avgtime,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');title(v,'Time(UTC)','FontSize',16);
    end
    xlabel(['Salinity(PSU)'],'fontsize',16); ylabel('Depth (m)','fontsize',16);
    set(gca,'YDir','reverse');
    set(gca,'FontSize',16);
    set(gca,'FontSize',16);
    title({[fList(j,9:10) '/' fList(j,6:7)]});
    avgsalrange=max(avgsal)-min(avgsal);
    avgsallowerlim=min(avgsal)- (avgsalrange*0.1);
    avgsalupperlim=max(avgsal)+ (avgsalrange*0.1);
%     xlim([avgsallowerlim avgsalupperlim]);
%     NumTicks=4;
%     L=get(gca,'XLim');
%     set(gca,'XTick',linspace(L(1),L(2),NumTicks));
%     set(gca,'XTickLabel',linspace(L(1),L(2),NumTicks));
    rotateXLabels(gca(),90);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%   flou subplots %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if j < 9 & j>1;
    figure(110)
    subplot(2,Subplotcolumns,j);
    line(avgchl,avgdepth,'LineWidth',0.5,'LineStyle','--');hold on;
    errorbar(avgchl,avgdepth,(2*semdepth));
    herrorbar(avgchl,avgdepth,(2*semchl));
    if axT==1;
        scatter(avgchl,avgdepth,80,avgTstr,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');ylabel(v,'Time since start (Hours Minutes)','FontSize',16);
    end
    if axT==0;
        scatter(avgchl,avgdepth,80,avgtime,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');ylabel(v,'Time(Hours Minutes)','FontSize',16);
    end
    xlabel(['Chl'],'fontsize',16); ylabel('Depth (m)','fontsize',16);
    set(gca,'YDir','reverse');
    set(gca,'FontSize',16)
    set(gca,'FontSize',16)
    title([fList(j,1:10)]);
    
    avgchlrange=max(avgchl)-min(avgchl);
    avgchllowerlim=min(avgchl)- (avgchlrange*0.1);
    avgchlupperlim=max(avgchl)+ (avgchlrange*0.1);
    xlim([avgchllowerlim avgchlupperlim])
    ylim([0 5])
    NumTicks=4;
    L=get(gca,'XLim');
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    set(gca,'XTickLabel',linspace(L(1),L(2),NumTicks))
    rotateXLabels(gca(),90)
    end
    
     Tempcorr=avgCTDTemperature_interp_lag - avgEquilibrator_Temperature_interp;

     
    if j > 8;
    tyu=j-8;
    figure(111)
    subplot(2,Subplotcolumns,tyu)
    line(avgchl,avgdepth,'LineWidth',0.5,'LineStyle','--');hold on;
    errorbar(avgchl,avgdepth,(2*semdepth))
    herrorbar(avgchl,avgdepth,(2*semchl))
    if axT==1
        scatter(avgchl,avgdepth,80,avgTstr,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');ylabel(v,'Time since start (Hours Minutes)','FontSize',16);
    end
    if axT==0
        scatter(avgchl,avgdepth,80,avgtime,'filled');hold on;
        v=colorbar; cbdate(v,'HH:MM');ylabel(v,'Time(Hours Minutes)','FontSize',16);
    end
    xlabel(['chl'],'fontsize',16); ylabel('Depth (m)','fontsize',16);
    set(gca,'YDir','reverse');
    set(gca,'FontSize',16)
    set(gca,'FontSize',16)
    title([fList(j,1:10)]);
    
    avgchlrange=max(avgchl)-min(avgchl);
    avgchllowerlim=min(avgchl)- (avgchlrange*0.1);
    avgchlupperlim=max(avgchl)+ (avgchlrange*0.1);
    xlim([avgchllowerlim avgchlupperlim])
    ylim([0 5])
    NumTicks=4;
    L=get(gca,'XLim');
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    set(gca,'XTickLabel',linspace(L(1),L(2),NumTicks))
    rotateXLabels(gca(),90)
    end

    %these figures are fine but don't show anything useful relating to
    %profiles except that we were on station. Just extra computation for no
    %real reason tbh. 

    
% %     if j < 9;
% %         figure(112)
% %         subplot(2,Subplotcolumns,j)
% %         m_proj('mercator','lon',[-4.35 -4],'lat',[50.2 50.4]);
% %         %     clf
% %         m_grid('linestyle','none','tickdir','out','fontsize',22)
% %         m_patch(coast(:,1),coast(:,2),'b','EdgeColor','none'); hold on;
% %         %     p=m_contour(u,i,bathymetrycelticsea',[-80:5:-15],'ShowText','on');
% %         m_scatter(Longitude_interp(inds(1):indt(end))*-1,Latitude_interp(inds(1):indt(end)),10,Licor_Datetime_doy(inds(1):indt(end)))
% %         set(gca,'fontsize',22)
% %         y=ylabel('Latitude');
% %         set(y,'Units','Normalized','Position',[-0.1,0.5,0]);
% %         set(gca,'fontsize',22)
% %         title([fList(j,1:10)]);
% %         xlabel('Longitude');
% %         set(gca,'fontsize',22)
% %         set(gca,'fontsize',22)
% %         hold on
% %         %add L4
% %         [C,D]=m_ll2xy(-4.217,50.25);
% %         line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% %         [CC,D]=m_ll2xy(-4.217,50.25);
% %         text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',24);
% %         ring1=m_range_ring(-4.217,50.25,0.5);
% %         set(ring1,'color','r')
% %         ring2=m_range_ring(-4.217,50.25,1);
% %         set(ring2,'color','b')
% %         ring3=m_range_ring(-4.217,50.25,2);
% %         set(ring3,'color','g')
% %     end
% %     
% %     
% %     if j > 8;
% %         tyu=j-8;
% %         figure(113)
% %         subplot(2,Subplotcolumns,tyu)
% %         m_proj('mercator','lon',[-4.35 -4],'lat',[50.2 50.4]);
% %         %      clf
% %         m_grid('linestyle','none','tickdir','out','fontsize',22)
% %         m_patch(coast(:,1),coast(:,2),'b','EdgeColor','none'); hold on;
% %         %     p=m_contour(u,i,bathymetrycelticsea',[-80:5:-15],'ShowText','on');
% %         m_scatter(Longitude_interp(inds(1):indt(end))*-1,Latitude_interp(inds(1):indt(end)),10,Licor_Datetime_doy(inds(1):indt(end)))
% %         set(gca,'fontsize',22)
% %         y=ylabel('Latitude');
% %         set(y,'Units','Normalized','Position',[-0.1,0.5,0]);
% %         set(gca,'fontsize',22)
% %         title([fList(j,1:10)]);
% %         xlabel('Longitude');
% %         set(gca,'fontsize',22)
% %         set(gca,'fontsize',22)
% %         hold on
% %         %add L4
% %         [C,D]=m_ll2xy(-4.217,50.25);
% %         line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% %         [CC,D]=m_ll2xy(-4.217,50.25);
% %         text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',24);
% %         ring1=m_range_ring(-4.217,50.25,0.5);
% %         set(ring1,'color','r')
% %         ring2=m_range_ring(-4.217,50.25,1);
% %         set(ring2,'color','b')
% %         ring3=m_range_ring(-4.217,50.25,2);
% %         set(ring3,'color','g')
% %     end
    
    figure(114)
    if j==1 | j==2 | j==3 | j==11
        %dont plot when there is no underway data to compare against
    else
        counter=counter+1;
        Subplotcolumns=4;
         subplot(3,Subplotcolumns,counter)%plot patch first underneath.
        x=[indt(end) inds(1) inds(1) indt(end)];
        y= [ max(Licor_fco2_surface_lag) max(Licor_fco2_surface_lag) 0 0];
        patch(Licor_Datetime(x),y,colour_greyshade);
        hold on
        
        %return journey
        if j==6;
        elseif j==14;
        else
            %Underway return
            x=[undind4 undind3 undind3 undind4];
            y= [ max(Licor_fco2_surface_lag) max(Licor_fco2_surface_lag) 0 0];
            patch(Licor_Datetime(x),y,colour_lightblue);
        end
        
        %outward journey
        if  j==4; %no data
        elseif j==1 %no data
        elseif j==2 %no data
        elseif j==8 %not enough data
        elseif j==5 %mfc not working on way out
        elseif j==9 %water bad agreement
        elseif j==11 %water bad agreement
        elseif j==12 %water bad agreement
        elseif j==16 %no data
        else
            %underway out
            x=[undind2 undind1 undind1 undind2];
            y= [ max(Licor_fco2_surface_lag) max(Licor_fco2_surface_lag) 0 0];
            patch(Licor_Datetime(x),y,colour_lightblue);
        end
        %
        %         plot(Licor_Datetime,Licor_fco2_surface_lag,'k','LineWidth',2);

        
        %use this to plot all data
        %     plot(Licor_Datetime,Licor_fco2_combined,'k','LineWidth',2);
        %use these to plot only actual data
        plot(Licor_Datetime(undind1:undind2),Licor_fco2_combined(undind1:undind2),'k','LineWidth',2);
        plot(Licor_Datetime(undind3:undind4),Licor_fco2_combined(undind3:undind4),'k','LineWidth',2);
        plot(Licor_Datetime(inds(1):indt(end)),Licor_fco2_combined(inds(1):indt(end)),'k','LineWidth',2);
        if  j==2; %no data
        elseif j==1 %no data dont trust kept turning off
        elseif j==3 %no data
        elseif j==11 %no data
        else
            plot(Licor_Datetime,underwayfCO2_sw_interp,'g','LineWidth',2);
        end
        set(gca,'FontSize',12)
        dynamicDateTicks([], [], 'HH:MM');
        %add text labels
    title({[fList(j,9:10) '/' fList(j,6:7)]});
        ylabel(['CO_2' char(10) '(',num2str(micro_symbol),'atm)'],'FontSize',12);
        set(gca,'FontSize',12);
        plot(underway_DT,underwayfCO2_sw,'k*');
        xlim([Licor_Datetime(1) Licor_Datetime(end)])
        %set y limits for each subplot
        if j==4;
            ylim([330 370]);
        elseif j==5 ;
            ylim([320 370]);
        elseif j==6 ;
            ylim([335 390]);
        elseif j==7 ;
            ylim([340 380]);
        elseif j==8 ;
            ylim([340 420]);
        elseif j==9 ;
            ylim([350 410]);
        elseif j==10;
            ylim([335 445]);
        elseif j==11;
            ylim([330 415]);
        elseif j==12;
            ylim([320 360]);
        elseif j==13;
            ylim([370 420]);
        elseif j==14 ;
            ylim([380 420]);
        elseif j==15 ;
            ylim([370 440]);
        elseif j==16 ;
            ylim([370 460])
        else
        end
    end

    if j < 9;
    tyu=j;
    figure(115)
    subplot(Subplotcolumns,2,tyu)
    if j==1 | j==3 | j==5 | j==7
    ylabel({['Tidal' char(10) 'height(m)']},'FontSize',16);
    end
    if j ==7 | j==8
    xlabel('Time (UTC)','FontSize',16)    
    end
    x=[indt(end) inds(1) inds(1) indt(end)];
    y= [ max(Plymouth_tidal_height) max(Plymouth_tidal_height) 0 0];
    patch(Licor_Datetime(x),y,colour_greyshade);
    hold on
    [f ~]=find(Plymouth_tidal_height_dt>=Licor_Datetime(1)-0.25);
    tideind1=min(f);
    [h ~]=find(Plymouth_tidal_height_dt<Licor_Datetime(end)+0.25);
    tideind2=max(h);
    plot(Plymouth_tidal_height_dt(tideind1:tideind2),Plymouth_tidal_height(tideind1:tideind2),'-k','Linewidth',8);
    %add text labels
    ylim([0 7]);
    dynamicDateTicks([], [], 'HH');
    setDateAxes(gca,'XLim',[Plymouth_tidal_height_dt(tideind1) Plymouth_tidal_height_dt(tideind2)]);
    set(gca,'FontSize',16);
    set(gca,'FontSize',16);
    title({[fList(j,9:10) '/' fList(j,6:7)]});
    end

    if j > 8;
    tyu=j-8;
    figure(116);
    subplot(Subplotcolumns,2,tyu);
    if j >14
    xlabel('Time (UTC)','FontSize',16)  ;  
    end
    if j==9 | j==11 | j==13 | j==15
    ylabel({['Tidal' char(10) 'height(m)']},'FontSize',16);
    end
    x=[indt(end) inds(1) inds(1) indt(end)];
    y= [ max(Plymouth_tidal_height) max(Plymouth_tidal_height) 0 0];
    patch(Licor_Datetime(x),y,colour_greyshade);
    hold on
    [f ~]=find(Plymouth_tidal_height_dt>=Licor_Datetime(1)-0.25);
    tideind1=min(f);
    [h ~]=find(Plymouth_tidal_height_dt<Licor_Datetime(end)+0.25);
    tideind2=max(h);
    plot(Plymouth_tidal_height_dt(tideind1:tideind2),Plymouth_tidal_height(tideind1:tideind2),'-k','Linewidth',8);
    %add text labels
    ylim([0 7]);
    dynamicDateTicks([], [], 'HH');
    setDateAxes(gca,'XLim',[Plymouth_tidal_height_dt(tideind1) Plymouth_tidal_height_dt(tideind2)]);
    set(gca,'FontSize',16);
    set(gca,'FontSize',16);
    title({[fList(j,9:10) '/' fList(j,6:7)]});
    end
    
           [ ind_temp_lowtide_und,~  ]=find(Plymouth_tidal_height_dt(ind_low_tide)<undtimeavg);
       nearest_lowtide_before_und=Plymouth_tidal_height_dt(ind_low_tide(max(ind_temp_lowtide_und)));
       dt_since_lowtide_und_ret=undtimeavg-nearest_lowtide_before_und;
       

       %plotted in get_plot_seasonalstudy.m
 
       
% %        
% %     %plot of all co2 by lat/long
% %     figure(117)
% %     subsample=30;
% %     m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
% %     m_grid('linestyle','none','tickdir','out','fontsize',22)
% %     m_patch(coast(:,1),coast(:,2),'b','EdgeColor','none'); hold on;
% %      if j==1 | j==2 | j==3 | j==11 ;
% %           else
% %          if j==6;
% %          elseif j==14;
% %          else
% %         undtimeavg=mean(Licor_Datetime(undind3:undind4))
% %              
% %         
% % 
% %        
% %        if (dt_since_lowtide_und_ret>0 & dt_since_lowtide_und_ret <1/8);
% %            markersty='h'
% %        elseif   (dt_since_lowtide_und_ret>1/8 & dt_since_lowtide_und_ret <2/8);
% %            markersty='s'
% %        elseif   (dt_since_lowtide_und_ret>2/8 & dt_since_lowtide_und_ret <3/8);
% %            markersty='d'  
% %        elseif   (dt_since_lowtide_und_ret>3/8 & dt_since_lowtide_und_ret <4/8);
% %            markersty='p' 
% %        else
% %            markersty='p' ;
% %        end
% %  
% %         m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,Licor_fco2_combined(undind3:subsample:undind4),'marker',markersty)
% %         text(Longitude_interp(undind4)*-1,Latitude_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% %          end
% %          hold on
% %          if  j==4;
% %          elseif j==5 %mfc not working on way out
% %          elseif j==8 %not enough data
% %          elseif j==9 %water bad agreement
% %          elseif j==12 %water bad agreement
% %          elseif j==16 %no data
% %          else
% %              %underway out
% %         m_scatter(Longitude_interp(undind1:undind2)*-1,Latitude_interp((undind1:undind2)),0.5,Licor_fco2_combined(undind1:undind2))
% %          text(Longitude_interp(undind2)*-1,Latitude_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% % 
% %          end
% %     hold on
% %     set(gca,'fontsize',22)
% %     y=ylabel('Latitude');
% %     set(y,'Units','Normalized','Position',[-0.1,0.5,0]);
% %     set(gca,'fontsize',22)
% %     xlabel('Longitude');
% %     set(gca,'fontsize',22)
% %     set(gca,'fontsize',22)
% %     hold on
% %     %add L4
% %     [C,D]=m_ll2xy(-4.217,50.25);
% %     line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% %     [CC,D]=m_ll2xy(-4.217,50.25);
% %     text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',24);
% %     
% %     %add penlee
% %     [C,D]=m_ll2xy(-4.188801,50.318280);
% %     line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% %     [CC,D]=m_ll2xy(-4.188801,50.318280);
% %     text(CC,D,'Penlee','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',24);
% %     %add breakwater
% %     [C,D]=m_ll2xy(-4.158932,50.334499);
% %     line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% %     [CC,D]=m_ll2xy(-4.158932,50.334499);
% %     text(CC,D,'Breakwater','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',24);
% %     
% %     ring1=m_range_ring(-4.217,50.25,0.5);
% %     set(ring1,'color','r')
% %     ring2=m_range_ring(-4.217,50.25,1);
% %     set(ring2,'color','b')
% %     ring3=m_range_ring(-4.217,50.25,2);
% %     set(ring3,'color','g')
% %     hj=colorbar;
% %     ylabel(hj,['CO_2''(',num2str(micro_symbol),'atm)'],'FontSize',24);
% %      end

    
         %plotted in get_plot_seasonalstudy.m
% % 
% %      %figure of co2 rleative to L4
% %      figure(119)
% %      if j==1 | j==2 | j==3 | j==11 %no data
% %      else
% %          %return voyages
% %          if j==6;
% %          elseif j==14;
% %          else
% %              scatter(dist(undind3:undind4), Licor_fco2_combined(undind3:undind4)-Licor_fco2_combined(undind3),25,Licor_Datetime((undind3:undind4)))
% %              %add date label
% %              text(dist(undind4),Licor_fco2_combined(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% %          end
% %          hold on
% %          %Outward voyages
% %          if  j==4;
% %          elseif j==5 %mfc not working on way out
% %          elseif j==8 %not enough data
% %          elseif j==9 %water bad agreement
% %          elseif j==12 %water bad agreement
% %          elseif j==16 %no data
% %          elseif j==10 %short    
% %          elseif j==14 %short  
% %          elseif j==7 %no data
% %          elseif j==10 %no data
% %          elseif j==14 %no data     
% %          else
% %              %underway out
% %              scatter(dist(undind1:undind2), Licor_fco2_combined(undind1:undind2)-Licor_fco2_combined(undind2),25,Licor_Datetime((undind1:undind2)))
% %              %add date label
% %               text(dist(undind2),Licor_fco2_combined(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% %         end
% %          xlabel('Distance from L4 (km)');
% %          set(gca,'fontsize',22)
% %          set(gca,'fontsize',22)
% %          ylabel(['\Delta fCO_2_(_s_w_)' '(',num2str(micro_symbol),'atm)'],'FontSize',24);
% %          set(gca,'Xlim',[-1 11]);
% %           u=colorbar; cbdate(u,'dd/mm');
% %          ylabel(u,'Date(dd/mm)','FontSize',22);    hold on
% %          plot(zeros(length(-60:1:120)),-60:1:120,'k','LineWidth',4)
% %          plot(10*ones(length(-60:1:120)),-60:1:120,'k','LineWidth',4)
% %      end  
% %         
% %      
% %      load coast.dat %high resollution coastline for south west of uk
% %      load('Data\bathymetry.mat')
% %      u=longitudeceltic(longcelticsea)
% %      i=(latitudeceltic(latcelticsea))'
% %      
% %     %plot of all co2 by lat/long
% %     figure(120)
% %     m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
% %     m_grid('linestyle','none','tickdir','out','fontsize',22)
% %     m_contourf(u,i,bathymetrycelticsea',[0 0]); hold on;
% % p=m_contour(u,i,bathymetrycelticsea',[-50:5:20],'ShowText','on');
% % 
% %      if j==1 | j==2 | j==3 | j==11
% %           else
% %          if j==6;
% %          elseif j==14;
% %          else
% %         m_scatter(Longitude_interp((undind3:undind4))*-1,Latitude_interp((undind3:undind4)),0.5,SeaTemp_interp(undind3:undind4)-SeaTemp_interp(undind3)); hold on;
% % %         m_text(Longitude_interp(undind4)*-1,Latitude_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);         
% %          end
% %          hold on
% %          if  j==4;
% %          elseif j==5 %mfc not working on way out
% %          elseif j==8 %not enough data
% %          elseif j==9 %water bad agreement
% %          elseif j==12 %water bad agreement
% %          elseif j==16 %no data
% %          elseif j==10 %short    
% %          elseif j==14 %short       
% %          else
% %              %underway out
% %         m_scatter(Longitude_interp(undind1:undind2)*-1,Latitude_interp((undind1:undind2)),0.5,SeaTemp_interp(undind1:undind2)-SeaTemp_interp(undind2))
% % %         m_text(Longitude_interp(undind2)*-1,Latitude_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% %          end
% %     hold on
% %     set(gca,'fontsize',22)
% %     y=ylabel('Latitude');
% %     set(y,'Units','Normalized','Position',[-0.15,0.5,0]);
% %     x=xlabel('Longitude');
% %     set(x,'Units','Normalized','Position',[0.5,-0.1,0]);
% %     set(gca,'fontsize',22)
% %     set(gca,'fontsize',22)
% %     hold on
% %     
% %     
% %     
% %     
% %     %add L4
% %     [C,D]=m_ll2xy(-4.217,50.25);
% %     line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% %     [CC,D]=m_ll2xy(-4.217,50.25);
% %     text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',24);
% %        %add penlee
% %     [C,D]=m_ll2xy(-4.188801,50.318280);
% %     line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% %     [CC,D]=m_ll2xy(-4.188801,50.318280);
% %     text(CC,D,'Penlee','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',24);
% %     %add breakwater
% %     [C,D]=m_ll2xy(-4.158932,50.334499);
% %     line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% %     [CC,D]=m_ll2xy(-4.158932,50.334499);
% %     text(CC,D,'Breakwater','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',24);
% %     ring1=m_range_ring(-4.217,50.25,0.5);
% %     set(ring1,'color','r')
% %     ring2=m_range_ring(-4.217,50.25,1);
% %     set(ring2,'color','b')
% %     hj=colorbar;
% %     ylabel(hj,['Temperature (',num2str(degree_symbol),'C)'],'FontSize',24);
% %     
% %      end
% % 
% %      


    %plotted in get_plot_seasonalstudy.m
% % 
% %     figure(122)
% %      if j==1 | j==2 | j==3 | j==11 %no data
% %      else
% %          if j==6;
% %          elseif j==14;
% %          else
% %              scatter(dist(undind3:undind4), SeaTemp_interp(undind3:undind4)-SeaTemp_interp(undind3),25,Licor_Datetime((undind3:undind4)))
% %              %add date label
% %              text(dist(undind4),SeaTemp_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% %          end
% %          hold on
% %          if  j==4;
% %          elseif j==5 %mfc not working on way out
% %          elseif j==8 %not enough data
% %          elseif j==9 %water bad agreement
% %          elseif j==12 %water bad agreement
% %          elseif j==16 %no data
% %          elseif j==7 %no data
% %          elseif j==10 %no data
% %          elseif j==14 %no data     
% %          else
% %              %underway out
% %              scatter(dist(undind1:undind2), SeaTemp_interp(undind1:undind2)-SeaTemp_interp(undind2),25,Licor_Datetime((undind1:undind2)))
% %              %add date label
% %               text(dist(undind2),SeaTemp_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% %         end
% %          xlabel('Distance from L4 (km)');
% %          set(gca,'fontsize',22)
% %          set(gca,'fontsize',22)
% %          ylabel(['\Delta Temperature (',num2str(degree_symbol),'C)'],'FontSize',24);
% %           set(gca,'Xlim',[-1 11]);
% %          u=colorbar; cbdate(u,'dd/mm');
% %          ylabel(u,'Date(dd/mm)','FontSize',22);    hold on
% %          plot(zeros(length(-1.5:1:1.5)),-1.5:1:1.5,'k','LineWidth',4)
% %          plot(10*ones(length(-1.5:1:1.5)),-1.5:1:1.5,'k','LineWidth',4)
% %          
% %      end
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%% Save ouput as loadable mat files
% %     if ~exist('Data', 'dir')
% %         mkdir('Data');
% %     end
% %     save([pwd '/Data/' fList(j,1:10) '.mat']);

 
    %plotted in get_plot_seasonalstudy.m

% 
%     %plot of all salinity by lat/long
%     figure(123)
%     %%
%     m_proj('Lambert','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
%     m_grid('linestyle','none','tickdir','out','fontsize',16)
%     m_patch(coast(:,1),coast(:,2),'b','EdgeColor','none'); hold on;
%      if j==1 | j==2 | j==3 | j==11
%           else
%          if j==6;
%          elseif j==14;
%          else
%         m_scatter(Longitude_interp((undind3:undind4))*-1,Latitude_interp((undind3:undind4)),0.5,Salinity_interp(undind3:undind4)); hold on;
% %         m_text(Longitude_interp(undind4)*-1,Latitude_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
%          end
%          hold on
%          if  j==4;
%          elseif j==5 %mfc not working on way out
%          elseif j==8 %not enough data
%          elseif j==9 %water bad agreement
%          elseif j==12 %water bad agreement
%          elseif j==16 %no data
%          elseif j==7 %no data
%          elseif j==10 %no data
%          elseif j==14 %no data     
%          else
%              %underway out
%         m_scatter(Longitude_interp(undind1:undind2)*-1,Latitude_interp((undind1:undind2)),0.5,Salinity_interp(undind1:undind2))
% %         m_text(Longitude_interp(undind2)*-1,Latitude_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% 
%          end
%     hold on
%     y=ylabel('Latitude');
%     x=xlabel('Longitude');
%     set(gca,'fontsize',24)
%     set(gca,'fontsize',24)
%     hold on
%     %add L4
%     [C,D]=m_ll2xy(-4.217,50.25);
%     line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
%     [CC,D]=m_ll2xy(-4.217,50.25);
%     text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',24);
%     ring1=m_range_ring(-4.217,50.25,0.5);
%     set(ring1,'color','r')
%     ring2=m_range_ring(-4.217,50.25,1);
%     set(ring2,'color','b')
%     hj=colorbar;
%     ylabel(hj,['Salinity (PSU)'],'FontSize',24);
%     %add penlee
%     [C,D]=m_ll2xy(-4.188801,50.318280);
%     line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
%     [CC,D]=m_ll2xy(-4.188801,50.318280);
%     text(CC,D,'Penlee','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',24);
%     %add breakwater
%     [C,D]=m_ll2xy(-4.158932,50.334499);
%     line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
%     [CC,D]=m_ll2xy(-4.158932,50.334499);
%     text(CC,D,'Breakwater','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',24);
%      end
% %     subplot(1,2,2)
% %     m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
% %     m_grid('linestyle','none','tickdir','out','fontsize',22)
% %     m_patch(coast(:,1),coast(:,2),'b','EdgeColor','none'); hold on;
% %      if j==1 | j==2 | j==3 | j==11
% %           else
% %          if j==6;
% %          elseif j==14;
% %          else
% %         m_scatter(Longitude_interp((undind3:undind4))*-1,Latitude_interp((undind3:undind4)),0.5,Salinity_interp(undind3:undind4)-Salinity_interp(undind3)); hold on;
% %         m_text(Longitude_interp(undind4)*-1,Latitude_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% %          end
% %          hold on
% %          if  j==4;
% %          elseif j==5 %mfc not working on way out
% %          elseif j==8 %not enough data
% %          elseif j==9 %water bad agreement
% %          elseif j==12 %water bad agreement
% %          elseif j==16 %no data
% %          else
% %              %underway out
% %         m_scatter(Longitude_interp(undind1:undind2)*-1,Latitude_interp((undind1:undind2)),0.5,Salinity_interp(undind1:undind2)-Salinity_interp(undind2))
% %         m_text(Longitude_interp(undind2)*-1,Latitude_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% % 
% %          end
% %     hold on
% %     set(gca,'fontsize',22)
% %     y=ylabel('Latitude');
% %     set(y,'Units','Normalized','Position',[-0.1,0.5,0]);
% %     set(gca,'fontsize',22)
% %     xlabel('Longitude');
% %     set(gca,'fontsize',22)
% %     set(gca,'fontsize',22)
% %     hold on
% %     %add L4
% %     [C,D]=m_ll2xy(-4.217,50.25);
% %     line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% %     [CC,D]=m_ll2xy(-4.217,50.25);
% %     text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',24);
% %     ring1=m_range_ring(-4.217,50.25,0.5);
% %     set(ring1,'color','r')
% %     ring2=m_range_ring(-4.217,50.25,1);
% %     set(ring2,'color','b')
% %     ring3=m_range_ring(-4.217,50.25,2);
% %     set(ring3,'color','g')
% %     hj=colorbar;
% %     ylabel(hj,['Salinity (PSU)'],'FontSize',24);
% %      end

    %plotted in get_plot_seasonalstudy.m


% %      
% %     figure(125)
% %      if j==1 | j==2 | j==3 | j==11 %no data
% %      else
% %          if j==6;
% %          elseif j==14;
% %          else
% %              plot(dist(undind3:undind4), Salinity_interp(undind3:undind4)-Salinity_interp(undind3))
% %              %add date label
% %              text(dist(undind4),Salinity_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% %          end
% %          hold on
% %          if  j==4;
% %          elseif j==5 %mfc not working on way out
% %          elseif j==8 %not enough data
% %          elseif j==9 %water bad agreement
% %          elseif j==12 %water bad agreement
% %          elseif j==16 %no data
% %          else
% %              %underway out
% %              plot(dist(undind1:undind2), Salinity_interp(undind1:undind2)-Salinity_interp(undind2))
% %              %add date label
% %              text(dist(undind2),Salinity_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% %         end
% %          xlabel('Distance from L4 (km)');
% %          set(gca,'fontsize',22)
% %          set(gca,'fontsize',22)
% %          ylabel(['Salinity (PSU)'],'FontSize',24);
% %          plot(zeros(length(-1:1:0.2)),-1:1:0.2,'k','LineWidth',4)
% %          plot(10*ones(length(-1:1:0.2)),-1:1:0.2,'k','LineWidth',4)
% %      end
     
       
    %plotted in get_plot_seasonalstudy.m


% %      figure(126)
% %      subsample=30;
% %     m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
% %     m_grid('linestyle','none','xtick',[],'ytick',[])
% %     m_patch(coast(:,1),coast(:,2),'b','EdgeColor','none'); hold on;
% %      if j==1 | j==2 | j==3 | j==11 | j==14
% %      else
% %         counterF=counterF+1;
% %         Subplotcolumns=5;
% %         subplot(3,Subplotcolumns,counterF)%plot patch first underneath.  
% %         m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
% %         m_grid('linestyle','none','xtick',[],'ytick',[])
% %         m_patch(coast(:,1),coast(:,2),'b','EdgeColor','none'); hold on;
% %          if j==6;
% %          elseif j==14;
% %          else
% %         undtimeavg=mean(Licor_Datetime(undind3:undind4))
% %        [ ind_temp_lowtide_und,~  ]=find(Plymouth_tidal_height_dt(ind_low_tide)<undtimeavg);
% %        nearest_lowtide_before_und=Plymouth_tidal_height_dt(ind_low_tide(max(ind_temp_lowtide_und)));
% %        dt_since_lowtide_und_ret=undtimeavg-nearest_lowtide_before_und;
% %        dt_since_lowtide_und_return=[dt_since_lowtide_und_return,dt_since_lowtide_und_ret];
% %        if dt_since_lowtide_und_ret>10.5/24 || dt_since_lowtide_und_ret <1.5/24; %9.5 - 1.5 hrs
% %            markersty='h'
% %        elseif   (dt_since_lowtide_und_ret>1.5/24 & dt_since_lowtide_und_ret <4.5/24); %1.5 - 4.5
% %            markersty='s'
% %        elseif   (dt_since_lowtide_und_ret>4.5/24 & dt_since_lowtide_und_ret <7.5/24);% 4.5 - 7.5
% %            markersty='d'  
% %        elseif   (dt_since_lowtide_und_ret>7.5/24 & dt_since_lowtide_und_ret <10.5/24);%7.5 - 10.5
% %            markersty='p' 
% %        else
% %            markersty='p' ;
% %        end
% %        dt_since_lowtide_und_retHHMM=datestr(dt_since_lowtide_und_ret,'HH:MM')
% %         m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,Licor_fco2_combined(undind3:subsample:undind4),'marker',markersty)
% %         text(Longitude_interp(undind4)*-1,Latitude_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% %           title({[fList(j,9:10) '/' fList(j,6:7) 'LW+' dt_since_lowtide_und_retHHMM(1:5) 'hrs']},'FontSize',16);hold on
% %          set(gca,'fontsize',16);
% %        %add L4
% %        [C,D]=m_ll2xy(-4.217,50.25);line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% %        [CC,D]=m_ll2xy(-4.217,50.25);text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',16);
% %        ring1=m_range_ring(-4.217,50.25,0.5);set(ring1,'color','r');ring2=m_range_ring(-4.217,50.25,1);set(ring2,'color','b')
% %              hj=colorbar;
% %                ylabel(hj,['fCO_2_(_',num2str(micro_symbol),'_a_t_m_)'],'FontSize',12);
% %          end
% %          hold on
% %          
% %          if  j==4;
% %          elseif j==5 %mfc not working on way out
% %          elseif j==8 %not enough data
% % %          elseif j==9 %water bad agreement
% % %          elseif j==12 %water bad agreement
% %          elseif j==16 %no data
% %          elseif j==14 %no data    
% %          else
% %              %underway out
% %              if j==9 | j==13 | j==15
% %              counterF=counterF+1;
% %              end
% %        undtimeavg=mean(Licor_Datetime(undind1:undind2));
% %        [ ind_temp_lowtide_und,~  ]=find(Plymouth_tidal_height_dt(ind_low_tide)<undtimeavg);
% %        nearest_lowtide_before_und=Plymouth_tidal_height_dt(ind_low_tide(max(ind_temp_lowtide_und)));
% %        dt_since_lowtide_und_out=undtimeavg-nearest_lowtide_before_und;   
% %        dt_since_lowtide_und_outbound=[dt_since_lowtide_und_outbound,dt_since_lowtide_und_out];
% %        subplot(3,Subplotcolumns,counterF)%plot patch first underneath.  
% %         if dt_since_lowtide_und_out>10.5/24 || dt_since_lowtide_und_out <1.5/24; %9.5 - 1.5 hrs
% %            markersty='h'
% %        elseif   (dt_since_lowtide_und_out>1.5/24 & dt_since_lowtide_und_out <4.5/24); %1.5 - 4.5
% %            markersty='s'
% %        elseif   (dt_since_lowtide_und_out>4.5/24 & dt_since_lowtide_und_out <7.5/24);% 4.5 - 7.5
% %            markersty='d'  
% %        elseif   (dt_since_lowtide_und_out>7.5/24 & dt_since_lowtide_und_out <10.5/24);%7.5 - 10.5
% %            markersty='p' 
% %        else
% %            markersty='p' ;
% %         end  
% %         dt_since_lowtide_und_outHHMM=datestr(dt_since_lowtide_und_out,'HH:MM')
% %         m_scatter(Longitude_interp(undind1:subsample:undind2)*-1,Latitude_interp((undind1:subsample:undind2)),25,Licor_fco2_combined(undind1:subsample:undind2))
% %         text(Longitude_interp(undind2)*-1,Latitude_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% %          title({[fList(j,9:10) '/' fList(j,6:7)] 'LW+' dt_since_lowtide_und_outHHMM(1:5) 'hrs'},'FontSize',20);hold on
% %          end
% %     hold on
% % %     y=ylabel('Latitude','fontsize',16);
% % %     xlabel('Longitude','fontsize',16);
% %     set(gca,'fontsize',16)
% %     set(gca,'fontsize',16)
% %     %add L4
% %     [C,D]=m_ll2xy(-4.217,50.25);
% %     line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% %     [CC,D]=m_ll2xy(-4.217,50.25);
% %     text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',16);
% %     ring1=m_range_ring(-4.217,50.25,0.5);
% %     set(ring1,'color','r')
% %     ring2=m_range_ring(-4.217,50.25,1);
% %     set(ring2,'color','b')
% %     hj=colorbar;
% %      ylabel(hj,['fCO_2_(_',num2str(micro_symbol),'_a_t_m_)'],'FontSize',16);
% %      end 
        
     %plotted in get_plot_seasonalstudy.m
       
% %     
% %      figure(127)
% %      subsample=30;
% %     m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
% %     m_grid('linestyle','none','xtick',[],'ytick',[])
% %     m_patch(coast(:,1),coast(:,2),'b','EdgeColor','none'); hold on;
% %      if j==1 | j==2 | j==3 | j==11 | j==14
% %      else
% %         counterG=counterG+1;
% %         Subplotcolumns=5;
% %         subplot(3,Subplotcolumns,counterG)%plot patch first underneath.  
% %         m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
% %         m_grid('linestyle','none','xtick',[],'ytick',[])
% %         m_patch(coast(:,1),coast(:,2),'b','EdgeColor','none'); hold on;
% %          if j==6;
% %          elseif j==14;
% %          else
% %         undtimeavg=mean(Licor_Datetime(undind3:undind4))
% %        [ ind_temp_lowtide_und,~  ]=find(Plymouth_tidal_height_dt(ind_low_tide)<undtimeavg);
% %        nearest_lowtide_before_und=Plymouth_tidal_height_dt(ind_low_tide(max(ind_temp_lowtide_und)));
% %        dt_since_lowtide_und_ret=undtimeavg-nearest_lowtide_before_und;
% %        if dt_since_lowtide_und_ret>10.5/24 || dt_since_lowtide_und_ret <1.5/24; %9.5 - 1.5 hrs
% %            markersty='h'
% %        elseif   (dt_since_lowtide_und_ret>1.5/24 & dt_since_lowtide_und_ret <4.5/24); %1.5 - 4.5
% %            markersty='s'
% %        elseif   (dt_since_lowtide_und_ret>4.5/24 & dt_since_lowtide_und_ret <7.5/24);% 4.5 - 7.5
% %            markersty='d'  
% %        elseif   (dt_since_lowtide_und_ret>7.5/24 & dt_since_lowtide_und_ret <10.5/24);%7.5 - 10.5
% %            markersty='p' 
% %        else
% %            markersty='p' ;
% %        end
% %        dt_since_lowtide_und_retHHMM=datestr(dt_since_lowtide_und_ret,'HH:MM')
% %        m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,Salinity_interp(undind3:subsample:undind4),'marker',markersty)
% %        text(Longitude_interp(undind4)*-1,Latitude_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% %        title({[fList(j,9:10) '/' fList(j,6:7) '  LW+' dt_since_lowtide_und_retHHMM(1:5) 'hrs']},'FontSize',18);hold on
% %        hj=colorbar; ylabel(hj,['Salinity (PSU)'],'FontSize',12);    set(gca,'fontsize',16);
% %        %add L4
% %        [C,D]=m_ll2xy(-4.217,50.25);line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% %        [CC,D]=m_ll2xy(-4.217,50.25);text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',16);
% %        ring1=m_range_ring(-4.217,50.25,0.5);set(ring1,'color','r');ring2=m_range_ring(-4.217,50.25,1);set(ring2,'color','b')
% %          end
% %          hold on
% %          if  j==4;
% %          elseif j==5 %mfc not working on way out
% %          elseif j==8 %not enough data
% % %          elseif j==9 %water bad agreement elseif j==12 %water bad
% % %          agreement
% %          elseif j==14 %no data    
% %          elseif j==16 %no data
% %          else
% %              %underway out underway out
% %              if j==9 | j==13 | j==15
% %              counterG=counterG+1;      
% %              end
% %        subplot(3,Subplotcolumns,counterG)%plot patch first underneath.  
% %         if dt_since_lowtide_und_out>10.5/24 || dt_since_lowtide_und_out <1.5/24; %9.5 - 1.5 hrs
% %            markersty='h'
% %        elseif   (dt_since_lowtide_und_out>1.5/24 & dt_since_lowtide_und_out <4.5/24); %1.5 - 4.5
% %            markersty='s'
% %        elseif   (dt_since_lowtide_und_out>4.5/24 & dt_since_lowtide_und_out <7.5/24);% 4.5 - 7.5
% %            markersty='d'  
% %        elseif   (dt_since_lowtide_und_out>7.5/24 & dt_since_lowtide_und_out <10.5/24);%7.5 - 10.5
% %            markersty='p' 
% %        else
% %            markersty='p' ;
% %         end   
% %         dt_since_lowtide_und_outHHMM=datestr(dt_since_lowtide_und_out,'HH:MM')
% %         m_scatter(Longitude_interp(undind1:subsample:undind2)*-1,Latitude_interp((undind1:subsample:undind2)),25,Salinity_interp(undind1:subsample:undind2))
% %         text(Longitude_interp(undind2)*-1,Latitude_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% %         title({[fList(j,9:10) '/' fList(j,6:7) '  LW+' dt_since_lowtide_und_outHHMM(1:5) 'hrs']},'FontSize',20);hold on
% %          end
% %     hold on
% % %     y=ylabel('Latitude','fontsize',16);
% % %     xlabel('Longitude','fontsize',16);
% %     set(gca,'fontsize',16)
% %     set(gca,'fontsize',16)
% %     %add L4
% %     [C,D]=m_ll2xy(-4.217,50.25);
% %     line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% %     [CC,D]=m_ll2xy(-4.217,50.25);
% %     text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',16);
% %     ring1=m_range_ring(-4.217,50.25,0.5);
% %     set(ring1,'color','r')
% %     ring2=m_range_ring(-4.217,50.25,1);
% %     set(ring2,'color','b')
% %     hj=colorbar;
% %     ylabel(hj,['Salinity (PSU)'],'FontSize',16);
% %      end    
     
         %plotted in get_plot_seasonalstudy.m
% % 
% %     figure(128)
% %      subsample=30;
% %     m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
% %     m_grid('linestyle','none','xtick',[],'ytick',[])
% %     m_patch(coast(:,1),coast(:,2),'b','EdgeColor','none'); hold on;
% %      if j==1 | j==2 | j==3 | j==11 | j==14
% %      else
% %         counterH=counterH+1;
% %         Subplotcolumns=5;
% %         subplot(3,Subplotcolumns,counterH)%plot patch first underneath.  
% %         m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
% %         m_grid('linestyle','none','xtick',[],'ytick',[])
% %         m_patch(coast(:,1),coast(:,2),'b','EdgeColor','none'); hold on;
% %          if j==6;
% %          elseif j==14;
% %          else
% %         undtimeavg=mean(Licor_Datetime(undind3:undind4))
% %        [ ind_temp_lowtide_und,~  ]=find(Plymouth_tidal_height_dt(ind_low_tide)<undtimeavg);
% %        nearest_lowtide_before_und=Plymouth_tidal_height_dt(ind_low_tide(max(ind_temp_lowtide_und)));
% %        dt_since_lowtide_und_ret=undtimeavg-nearest_lowtide_before_und;
% %        if dt_since_lowtide_und_ret>10.5/24 || dt_since_lowtide_und_ret <1.5/24; %9.5 - 1.5 hrs
% %            markersty='h'
% %        elseif   (dt_since_lowtide_und_ret>1.5/24 & dt_since_lowtide_und_ret <4.5/24); %1.5 - 4.5
% %            markersty='s'
% %        elseif   (dt_since_lowtide_und_ret>4.5/24 & dt_since_lowtide_und_ret <7.5/24);% 4.5 - 7.5
% %            markersty='d'  
% %        elseif   (dt_since_lowtide_und_ret>7.5/24 & dt_since_lowtide_und_ret <10.5/24);%7.5 - 10.5
% %            markersty='p' 
% %        else
% %            markersty='p' ;
% %        end
% %        dt_since_lowtide_und_retHHMM=datestr(dt_since_lowtide_und_ret,'HH:MM')
% %        m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,Licor_fco2_combined(undind3:subsample:undind4),'marker',markersty)
% %        text(Longitude_interp(undind4)*-1,Latitude_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% %        title({[fList(j,9:10) '/' fList(j,6:7) '  LW+' dt_since_lowtide_und_retHHMM(1:5) 'hrs']},'FontSize',18);hold on
% %        hj=colorbar;  ylabel(hj,['fCO_2_(_',num2str(micro_symbol),'_a_t_m_)'],'FontSize',12); set(gca,'fontsize',16);
% %        %add L4
% %        [C,D]=m_ll2xy(-4.217,50.25);line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% %        [CC,D]=m_ll2xy(-4.217,50.25);text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',16);
% %        ring1=m_range_ring(-4.217,50.25,0.5);set(ring1,'color','r');ring2=m_range_ring(-4.217,50.25,1);set(ring2,'color','b')
% %          end
% %          hold on
% %          if  j==4;
% %          elseif j==5 %mfc not working on way out
% %          elseif j==8 %not enough data
% % %          elseif j==9 %water bad agreement
% % %          elseif j==12 %water bad agreement
% %          elseif j==14 %no data    
% %          elseif j==16 %no data
% %          else
% %              %underway out
% %              %underway out
% %              if j==9 | j==13 | j==15
% %              counterH=counterH+1;      
% %              end
% %        subplot(3,Subplotcolumns,counterH)%plot patch first underneath.  
% %         if dt_since_lowtide_und_out>10.5/24 || dt_since_lowtide_und_out <1.5/24; %9.5 - 1.5 hrs
% %            markersty='h'
% %        elseif   (dt_since_lowtide_und_out>1.5/24 & dt_since_lowtide_und_out <4.5/24); %1.5 - 4.5
% %            markersty='s'
% %        elseif   (dt_since_lowtide_und_out>4.5/24 & dt_since_lowtide_und_out <7.5/24);% 4.5 - 7.5
% %            markersty='d'  
% %        elseif   (dt_since_lowtide_und_out>7.5/24 & dt_since_lowtide_und_out <10.5/24);%7.5 - 10.5
% %            markersty='p' 
% %        else
% %            markersty='p' ;
% %         end   
% %         dt_since_lowtide_und_outHHMM=datestr(dt_since_lowtide_und_out,'HH:MM')
% %         m_scatter(Longitude_interp(undind1:subsample:undind2)*-1,Latitude_interp((undind1:subsample:undind2)),25,Licor_fco2_combined(undind1:subsample:undind2))
% %         text(Longitude_interp(undind2)*-1,Latitude_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% %         title({[fList(j,9:10) '/' fList(j,6:7) '  LW+' dt_since_lowtide_und_outHHMM(1:5) 'hrs']},'FontSize',20);hold on
% %          end
% %     hold on
% % %     y=ylabel('Latitude','fontsize',16);
% % %     xlabel('Longitude','fontsize',16);
% %     set(gca,'fontsize',16)
% %     set(gca,'fontsize',16)
% %     %add L4
% %     [C,D]=m_ll2xy(-4.217,50.25);
% %     line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% %     [CC,D]=m_ll2xy(-4.217,50.25);
% %     text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',16);
% %     ring1=m_range_ring(-4.217,50.25,0.5);
% %     set(ring1,'color','r')
% %     ring2=m_range_ring(-4.217,50.25,1);
% %     set(ring2,'color','b')
% %     hj=colorbar;
% %      ylabel(hj,['fCO_2_(_',num2str(micro_symbol),'_a_t_m_)'],'FontSize',16);
% %      end     
           
        
        
        counterx=[];
        figure(130)
        if j==1 | j==2 | j==3 | j==11
            %dont plot when there is no underway data to compare against
        else
            counterx=counterx+1;
            Subplotcolumns=4;
            subplot(3,Subplotcolumns,counter)%plot patch first underneath.
            if j==6;
            elseif j==14;
            else
                %Underway return
                scatter(Licor_fco2_combined(undind3:undind4)-Licor_fco2_combined(undind3),Salinity_interp(undind3:undind4)-Salinity_interp(undind3),8,dist(undind3:undind4))
                xlabel(['fCO_2(',num2str(micro_symbol),'atm)'],'fontsize',22)
                set(gca,'fontsize',22)
                set(gca,'fontsize',22)
                ylabel('Salinity','fontsize',22)
                hj=colorbar;
                title(hj,'Distance from L4 (km)','FontSize',24);
            end
            if  j==4; %no data
            elseif j==1 %no data
            elseif j==2 %no data
            elseif j==8 %not enough data
            elseif j==5 %mfc not working on way out
            elseif j==9 %water bad agreement
            elseif j==11 %water bad agreement
            elseif j==12 %water bad agreement
            elseif j==16 %no data
            else
                %underway out
                scatter(Licor_fco2_combined(undind1:undind2)-Licor_fco2_combined(undind2),Salinity_interp(undind1:undind2)-Salinity_interp(undind2),8,dist(undind1:undind2))
                xlabel(['fCO_2(',num2str(micro_symbol),'atm)'],'fontsize',22)
                set(gca,'fontsize',22)
                set(gca,'fontsize',22)
                ylabel('Salinity','fontsize',22)
                hj=colorbar;
                title(hj,'Distance from L4 (km)','FontSize',24);
            end
        end
     
        countery=[];
        figure(131)
        if j==1 | j==2 | j==3 | j==11
            %dont plot when there is no underway data to compare against
        else
            countery=countery+1;
            Subplotcolumns=4;
            subplot(3,Subplotcolumns,counter)%plot patch first underneath.
            if j==6;
            elseif j==14;
            else
                %Underway return
                scatter(SeaTemp_interp(undind3:undind4)-SeaTemp_interp(undind3),Salinity_interp(undind3:undind4)-Salinity_interp(undind3),8,Licor_fco2_combined(undind3:undind4)-Licor_fco2_combined(undind3))
                xlabel(['Temperature (',num2str(degree_symbol),'C)'],'fontsize',22)
                set(gca,'fontsize',22)
                set(gca,'fontsize',22)
                ylabel('Salinity','fontsize',22)
                hj=colorbar;
                title(hj,['fCO_2(',num2str(micro_symbol),'atm)'],'FontSize',24);
            end
            if  j==4; %no data
            elseif j==1 %no data
            elseif j==2 %no data
            elseif j==8 %not enough data
            elseif j==5 %mfc not working on way out
            elseif j==9 %water bad agreement
            elseif j==11 %water bad agreement
            elseif j==12 %water bad agreement
            elseif j==16 %no data
            else
                %underway out
                scatter(SeaTemp_interp(undind1:undind2)-SeaTemp_interp(undind2),Salinity_interp(undind1:undind2)-Salinity_interp(undind2),8,Licor_fco2_combined(undind1:undind2)-Licor_fco2_combined(undind2))
                xlabel(['Temperature (',num2str(degree_symbol),'C)'],'fontsize',22)
                set(gca,'fontsize',22)
                set(gca,'fontsize',22)
                ylabel('Salinity','fontsize',22)
                hj=colorbar;
                title(hj,['fCO_2(',num2str(micro_symbol),'atm)'],'FontSize',24);
            end
        end
     
%      figure(131)
%      subplot(Subplotcolumns,2,tyu)
%      scatter(Licor_fco2_combined(undind3:undind4)-Licor_fco2_combined(undind3),SeaTemp_interp(undind3:undind4)-SeaTemp_interp(undind3),8,dist(undind3:undind4))
%      xlabel(['fCO_2(',num2str(micro_symbol),'atm)'],'fontsize',22)
%      set(gca,'fontsize',22)
%      set(gca,'fontsize',22)
%      ylabel(['Temperature (',num2str(degree_symbol),'C)'],'FontSize',24);
%      hj=colorbar;
%      title(hj,'Distance from L4 (km)','FontSize',24);
%      
%      figure(132)
%      subplot(Subplotcolumns,2,tyu)
%      scatter(Salinity_interp(undind3:undind4)-Salinity_interp(undind3),SeaTemp_interp(undind3:undind4)-SeaTemp_interp(undind3),8,Licor_fco2_combined(undind3:undind4)-Licor_fco2_combined(undind3))
%      xlabel('Salinity','fontsize',22)
%      set(gca,'fontsize',22)
%      set(gca,'fontsize',22)
%      ylabel(['Temperature (',num2str(degree_symbol),'C)'],'FontSize',24);
%      hj=colorbar;
%      title(hj,['fCO_2(',num2str(micro_symbol),'atm)'],'FontSize',24);
     

% jumbo figure for paper combining all the transects as distance plots




% subplot = @(m,n,p) subtightplot(m,n,p,[0.035 0.05], [0.07 0.03], [0.04 0.02]); 

% 
% figure(132)
% seasonalstudy_dt_split1= datenum('2016-04-28 00:00:00','yyyy-mm-dd HH:MM:SS');
% seasonalstudy_dt_split2= datenum('2016-07-01 00:00:00','yyyy-mm-dd HH:MM:SS');
% seasonalstudy_dt_split3= datenum('2016-07-01 00:00:00','yyyy-mm-dd HH:MM:SS');
% seasonalstudy_dt_split4= datenum('2016-10-22 00:00:00','yyyy-mm-dd HH:MM:SS');
% underway_co2_pt1ss=find (underwayLat>50.20 & underway_DT<seasonalstudy_dt_split2 & underway_DT>seasonalstudy_dt_split1)
% underway_co2_pt2ss=find (underwayLat>50.20 & underway_DT>seasonalstudy_dt_split3 & underway_DT<seasonalstudy_dt_split4)
% 
% subplot(4, 2, 1,'Position',[0.03 0.77  0.45 0.23])
% hold on; 
% if j==1 | j==2 | j==3 | j==11 %no data
% else
% if j==6;
% elseif j==14;
% elseif j==9 | j==10 | j==11 | j==12 | j==13 | j==14 | j==15 | j==16 % second half of data block
% else
% scatter(dist(undind3:undind4), Salinity_interp(undind3:undind4),25,Licor_Datetime((undind3:undind4)))
% hold on; 
% %add date label
% text(dist(undind4),Salinity_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% end
% hold on
% if  j==4;
% elseif j==5 %mfc not working on way out
% elseif j==8 %not enough data
% elseif j==9 %water bad agreement
% elseif j==12 %water bad agreement
% elseif j==16 %no data
% elseif j==7 %no data
% elseif j==10 %no data
% elseif j==14 %no data
% elseif j==9 | j==10 | j==11 | j==12 | j==13 | j==14 | j==15 | j==16 % second half of data block
% else
% %underway out
% scatter(dist(undind1:undind2), Salinity_interp(undind1:undind2),25,Licor_Datetime((undind1:undind2)))
% hold on; 
% %add date label
% text(dist(undind2),Salinity_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% end
% set(gca,'fontsize',12)
% set(gca,'fontsize',12)
% ylabel(['Salinity'],'FontSize',12);
% % set(gca,'Xlim',[-1 11]);
% % set(gca,'Ylim',[34 35.2]);
% u=colorbar; cbdate(u,'dd/mm');
% title(u,'Date(dd/mm)','FontSize',12);    hold on
% plot(zeros(length(34:0.1:35.4)),34:0.1:35.4,'k','LineWidth',2)
% plot(10*ones(length(34:0.1:35.4)),34:0.1:35.4,'k','LineWidth',2)
% set(gca,'Xticklabel',[]) %remove tick labels but keep marks
% end
% subplot(4,2,2)
% if j==1 | j==2 | j==3 | j==11 %no data
% else
% if j==6;
% elseif j==14;
% elseif j==1 | j==2 | j==3 | j==4 | j==5 | j==6 | j==7 | j==8; % second half of data block
% else
% scatter(dist(undind3:undind4), Salinity_interp(undind3:undind4),25,Licor_Datetime((undind3:undind4)))
% %add date label
% text(dist(undind4),Salinity_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% end
% hold on
% if  j==4;
% elseif j==5 %mfc not working on way out
% elseif j==8 %not enough data
% elseif j==9 %water bad agreement
% elseif j==12 %water bad agreement
% elseif j==16 %no data
% elseif j==7 %no data
% elseif j==10 %no data
% elseif j==14 %no data
% elseif j==1 | j==2 | j==3 | j==4 | j==5 | j==6 | j==7 | j==8; % second half of data block
% else
% %underway out
% scatter(dist(undind1:undind2), Salinity_interp(undind1:undind2),25,Licor_Datetime((undind1:undind2)))
% %add date label
% text(dist(undind2),Salinity_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% end
% set(gca,'fontsize',12)
% set(gca,'fontsize',12)
% %          ylabel(['Salinity'],'FontSize',12);
% set(gca,'Xlim',[-1 11]);
% set(gca,'Ylim',[34 35.2]);
% u=colorbar; cbdate(u,'dd/mm');
% title(u,'Date(dd/mm)','FontSize',12);    hold on
% plot(zeros(length(34:0.1:35.4)),34:0.1:35.4,'k','LineWidth',2)
% plot(10*ones(length(34:0.1:35.4)),34:0.1:35.4,'k','LineWidth',2)
% set(gca,'Xticklabel',[]) %remove tick labels but keep marks
% end
% subplot(4,2,3)
% if j==1 | j==2 | j==3 | j==11 %no data
% else
% if j==6;
% elseif j==14;
% elseif j==9 | j==10 | j==11 | j==12 | j==13 | j==14 | j==15 | j==16 % second half of data block
% else
% scatter(dist(undind3:undind4), SeaTemp_interp(undind3:undind4),25,Licor_Datetime((undind3:undind4)))
% %add date label
% text(dist(undind4),SeaTemp_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% end
% hold on
% if  j==4;
% elseif j==5 %mfc not working on way out
% elseif j==8 %not enough data
% elseif j==9 %water bad agreement
% elseif j==12 %water bad agreement
% elseif j==16 %no data
% elseif j==7 %no data
% elseif j==10 %no data
% elseif j==14 %no data
% elseif j==9 | j==10 | j==11 | j==12 | j==13 | j==14 | j==15 | j==16 % second half of data block
% else
% %underway out
% scatter(dist(undind1:undind2), SeaTemp_interp(undind1:undind2),25,Licor_Datetime((undind1:undind2)))
% %add date label
% text(dist(undind2),SeaTemp_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% end
% set(gca,'fontsize',12)
% set(gca,'fontsize',12)
% ylabel(['Temperature (',num2str(degree_symbol),'C)'],'FontSize',12);
% set(gca,'Xlim',[-1 11]);
% set(gca,'Ylim',[13 18]);
% u=colorbar; cbdate(u,'dd/mm');
% title(u,'Date(dd/mm)','FontSize',12);    hold on
% plot(zeros(length(13:1:18)),13:1:18,'k','LineWidth',2);
% plot(10*ones(length(13:1:18)),13:1:18,'k','LineWidth',2);
% set(gca,'Xticklabel',[]) %remove tick labels but keep marks
% end
% subplot(4,2,4)
% if j==1 | j==2 | j==3 | j==11 %no data
% else
% if j==6;
% elseif j==14;
% elseif j==1 | j==2 | j==3 | j==4 | j==5 | j==6 | j==7 | j==8; % second half of data block
% else
% scatter(dist(undind3:undind4), SeaTemp_interp(undind3:undind4),25,Licor_Datetime((undind3:undind4)))
% %add date label
% text(dist(undind4),SeaTemp_interp(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% end
% hold on
% if  j==4;
% elseif j==5 %mfc not working on way out
% elseif j==8 %not enough data
% elseif j==9 %water bad agreement
% elseif j==12 %water bad agreement
% elseif j==16 %no data
% elseif j==7 %no data
% elseif j==10 %no data
% elseif j==14 %no data
% elseif j==1 | j==2 | j==3 | j==4 | j==5 | j==6 | j==7 | j==8; % second half of data block
% else
% %underway out
% scatter(dist(undind1:undind2), SeaTemp_interp(undind1:undind2),25,Licor_Datetime((undind1:undind2)))
% %add date label
% text(dist(undind2),SeaTemp_interp(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% end
% set(gca,'fontsize',12)
% set(gca,'fontsize',12)
% %          ylabel(['Temperature (',num2str(degree_symbol),'C)'],'FontSize',12);
% set(gca,'Xlim',[-1 11]);
% set(gca,'Ylim',[13 18]);
% u=colorbar; cbdate(u,'dd/mm');
% title(u,'Date(dd/mm)','FontSize',12);    hold on
% plot(zeros(length(13:1:18)),13:1:18,'k','LineWidth',2);
% plot(10*ones(length(13:1:18)),13:1:18,'k','LineWidth',2);
% set(gca,'Xticklabel',[]) %remove tick labels but keep marks
% end
% subplot(4,2,5)
% if j==1 | j==2 | j==3 | j==11 %no data
% else
% if j==6;
% elseif j==14;
% elseif j==9 | j==10 | j==11 | j==12 | j==13 | j==14 | j==15 | j==16 % second half of data block
% else
% scatter(dist(undind3:undind4), Licor_fco2_combined(undind3:undind4),25,Licor_Datetime((undind3:undind4)))
% %add date label
% text(dist(undind4),Licor_fco2_combined(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% end
% hold on
% if  j==4;
% elseif j==5 %mfc not working on way out
% elseif j==8 %not enough data
% elseif j==9 %water bad agreement
% elseif j==12 %water bad agreement
% elseif j==16 %no data
% elseif j==7 %no data
% elseif j==10 %no data
% elseif j==14 %no data
% elseif j==9 | j==10 | j==11 | j==12 | j==13 | j==14 | j==15 | j==16 % second half of data block
% else
% %underway out
% scatter(dist(undind1:undind2), Licor_fco2_combined(undind1:undind2),25,Licor_Datetime((undind1:undind2)))
% %add date label
% text(dist(undind2),Licor_fco2_combined(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% end
% set(gca,'fontsize',12)
% set(gca,'fontsize',12)
% ylabel(['fCO_2_(_s_w_)' '(',num2str(micro_symbol),'atm)'],'FontSize',12);
% set(gca,'Xlim',[-1 11]);
% u=colorbar; cbdate(u,'dd/mm');
% title(u,'Date(dd/mm)','FontSize',12);    hold on
% plot(zeros(length(200:1:600)),200:1:600,'k','LineWidth',2)
% plot(10*ones(length(200:1:600)),200:1:600,'k','LineWidth',2)
% set(gca,'Xticklabel',[]); %remove tick labels but keep marks
% end
% subplot(4,2,6)
% if j==1 | j==2 | j==3 | j==11 %no data
% else
% if j==6;
% elseif j==14;
% elseif j==1 | j==2 | j==3 | j==4 | j==5 | j==6 | j==7 | j==8; % second half of data block
% else
% scatter(dist(undind3:undind4), Licor_fco2_combined(undind3:undind4),25,Licor_Datetime((undind3:undind4)))
% %add date label
% text(dist(undind4),Licor_fco2_combined(undind4),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% end
% hold on
% if  j==4;
% elseif j==5 %mfc not working on way out
% elseif j==8 %not enough data
% elseif j==9 %water bad agreement
% elseif j==12 %water bad agreement
% elseif j==16 %no data
% elseif j==7 %no data
% elseif j==10 %no data
% elseif j==14 %no data
% elseif j==1 | j==2 | j==3 | j==4 | j==5 | j==6 | j==7 | j==8; % second half of data block
% else
% %underway out
% scatter(dist(undind1:undind2), Licor_fco2_combined(undind1:undind2),25,Licor_Datetime((undind1:undind2)))
% %add date label
% text(dist(undind2),Licor_fco2_combined(undind2),({[fList(j,9:10) '/' fList(j,6:7)]}),'HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',18);
% end
% set(gca,'fontsize',12)
% set(gca,'fontsize',12)
% %          ylabel(['fCO_2_(_s_w_)' '(',num2str(micro_symbol),'atm)'],'FontSize',12);
% set(gca,'Xlim',[-1 11]);
% u=colorbar; cbdate(u,'dd/mm');
% title(u,'Date(dd/mm)','FontSize',12);    hold on
% plot(zeros(length(200:1:600)),200:1:600,'k','LineWidth',2)
% plot(10*ones(length(200:1:600)),200:1:600,'k','LineWidth',2)
% set(gca,'Xticklabel',[]); %remove tick labels but keep marks
% end
% subplot(4,2,7)
% scatter(underway_dist(underway_co2_pt1ss), underwayfCO2_sw(underway_co2_pt1ss),25,underway_DT(underway_co2_pt1ss))
% xlabel('Distance from L4 (km)');
% set(gca,'fontsize',12)
% set(gca,'fontsize',12)
% ylabel(['fCO_2_(_s_w_)' '(',num2str(micro_symbol),'atm)'],'FontSize',12);
% set(gca,'Xlim',[-1 11]);
% u=colorbar; cbdate(u,'dd/mm');
% title(u,'Date(dd/mm)','FontSize',12);    hold on
% plot(zeros(length(200:1:600)),200:1:600,'k','LineWidth',2)
% plot(10*ones(length(200:1:600)),200:1:600,'k','LineWidth',2)
% subplot(4,2,8)
% scatter(underway_dist(underway_co2_pt2ss), underwayfCO2_sw(underway_co2_pt2ss),25,underway_DT(underway_co2_pt2ss))
% xlabel('Distance from L4 (km)');
% set(gca,'fontsize',12)
% set(gca,'fontsize',12)
% %          ylabel(['fCO_2_(_s_w_)' '(',num2str(micro_symbol),'atm)'],'FontSize',12);
% set(gca,'Xlim',[-1 11]);
% u=colorbar; cbdate(u,'dd/mm');
% title(u,'Date(dd/mm)','FontSize',12);    hold on
% plot(zeros(length(200:1:600)),200:1:600,'k','LineWidth',2)
% plot(10*ones(length(200:1:600)),200:1:600,'k','LineWidth',2)
%       
%    
%    
%    
   
       %plotted in get_plot_seasonalstudy.m
% % 
% % %plot ALL TRACKS
% % figure(1170)
% % subsample=30;
% % m_proj('mercator','lon',[-4.25 -4.1],'lat',[50.23 50.35]);
% % m_grid('linestyle','none','tickdir','out','fontsize',16)
% % m_patch(coast(:,1),coast(:,2),[0.4660, 0.6740, 0.1880],'EdgeColor','none'); hold on;
% % m_ruler([0.8 0.95],0.1,[0 1000 2000])
% % if j==1 | j==2 | j==3 | j==11
% % else
% %     if j==6;
% %     elseif j==14;
% %     else
% %         undtimeavg=mean(Licor_Datetime(undind3:undind4));
% %         
% %         [ ind_temp_lowtide_und,~  ]=find(Plymouth_tidal_height_dt(ind_low_tide)<undtimeavg);
% %         nearest_lowtide_before_und=Plymouth_tidal_height_dt(ind_low_tide(max(ind_temp_lowtide_und)));
% %         dt_since_lowtide_und_ret=undtimeavg-nearest_lowtide_before_und;
% %         %         m_plot(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),'-ok')
% %         m_scatter(Longitude_interp((undind3:subsample:undind4))*-1,Latitude_interp((undind3:subsample:undind4)),25,Licor_Datetime((undind3:subsample:undind4)))     
% %     end
% %     hold on
% %     if  j==4;
% %     elseif j==5 %mfc not working on way out
% %     elseif j==8 %not enough data
% %     elseif j==9 %water bad agreement
% %     elseif j==12 %water bad agreement
% %     elseif j==16 %no data
% %     else
% %         %underway out
% %         %         m_plot(Longitude_interp(undind1:undind2)*-1,Latitude_interp((undind1:undind2)),'-ok')
% %         m_scatter(Longitude_interp(undind1:subsample:undind2)*-1,Latitude_interp((undind1:subsample:undind2)),25,Licor_Datetime((undind1:subsample:undind2)))
% % 
% %     end
% %      u=colorbar; cbdate(u,'dd/mm');
% %      ylabel(u,'Date(dd/mm)','FontSize',22);    hold on
% %     set(gca,'fontsize',22)
% %     hold on
% %     y=ylabel('Latitude');
% %     x=xlabel('Longitude');
% %     set(gca,'fontsize',24)
% %     set(gca,'fontsize',24)
% %     hold on
% %     %add L4
% %     [C,D]=m_ll2xy(-4.217,50.25);
% %     line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% %     [CC,D]=m_ll2xy(-4.217,50.25);
% %     text(CC,D,'L4','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',24);
% %     ring1=m_range_ring(-4.217,50.25,0.5);
% %     set(ring1,'color','r')
% %     ring2=m_range_ring(-4.217,50.25,1);
% %     set(ring2,'color','b')
% %     %add penlee
% %     [C,D]=m_ll2xy(-4.188801,50.318280);
% %     line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% %     [CC,D]=m_ll2xy(-4.188801,50.318280);
% %     text(CC,D,'Penlee','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',24);
% %     %add breakwater
% %     [C,D]=m_ll2xy(-4.158932,50.334499);
% %     line(C,D,'marker','square','markersize',8,'MarkerFaceColor','k','color','k');
% %     [CC,D]=m_ll2xy(-4.158932,50.334499);
% %     text(CC,D,'Breakwater','HorizontalAlignment', 'right','VerticalAlignment', 'bottom','fontsize',24);
% % end


     

    %%%%%%%%%%%%%%%%%%%%%%%%%%% Save ouput as loadable mat files
    if ~exist('Data', 'dir')
        mkdir('Data');
    end
    save([pwd '/Data/' fList(j,1:10) '.mat']);


    j; 
    telapsed=toc(tstart);
   %close all
   end




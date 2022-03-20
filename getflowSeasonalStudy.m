function [flowDT, MeasuredFlow_interp, lagsecs,respsecs,smoothflow] = getflowSeasonalStudy(fullpathflow, Licor_Datetime,j,fList);

if exist(fullpathflow, 'file')
    delimiter = ' ';formatSpec = '%s%f%s%f%s%s%s%[^\n\r]'; fileID = fopen(fullpathflow,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);
    fclose(fileID);
    Voltage = dataArray{:, 2};
    Precalcflow = dataArray{:, 4};
    flowymd = cell2mat(dataArray{:, 6});
    flowhms = cell2mat(dataArray{:, 7});
    flowymdhms=[flowymd,flowhms];
    flowmisecsDT=(str2num(flowymdhms(:,20:24)))*(1/(60*60*24))/100000; %need to add microsecs or not monotomically increasing
    flowDT= (datenum(flowymdhms,'yyyy-mm-ddHH:MM:SS'))+flowmisecsDT;
    
    if fList(j,1:10)==('2016-08-10');
    [findvoltrow,~]=find(Voltage(10000:47400)>2);
    Voltage(findvoltrow)=[];
    flowDT(findvoltrow)=[];
    end
    [flowDT1, ia, ic] = unique(flowDT);

    
    MeasuredFlow=(0.5196*Voltage)-0.0356; %values from sensor calibration (see excel file)
    MeasuredFlow1=MeasuredFlow(ia);
    MeasuredFlow_interp = interp1(flowDT1,MeasuredFlow1,Licor_Datetime);
    k = find(isnan(MeasuredFlow_interp ))'; MeasuredFlow_interp (k) = 0;

    %smooth flow rate data
    sF = 1; %sampling frequency for pump in Hz
    % spec = hspec_ss(Precalcflow_interp,Precalcflow_interp,sF);
    %figure(999); loglog(spec(:,1),spec(:,1).*spec(:,2),'.'); %plot the frequency response and use it to set cutoff frequency
    N=4; Wn = 0.02/(sF/2);  % Normalized cutoff frequency
    [B,A] = butter(N,Wn);   %butterworth filter
    %fvtool(B,A);          %plot magnitude response
    smoothflow = filtfilt(B,A,MeasuredFlow_interp); %smooth data using filtfilt
    %calculate the lag from  the flow rate
    %volume of tubing in m^3
    tubinglength=53.2;

    tubevol=tubinglength*pi*((0.0127/2)^2); %diameter tubing is 1.27cm and 54m long tubing
    tubevolL=tubevol*1000;
    % clearvars filename delimiter formatSpec fileID dataArray ans;
    lagsecs=tubevolL./((smoothflow(:,1))/60); %lagtime in seconds - add to CTD
    respsecs=24; %response time in seconds - add to first CO2 point

    else
    flowDT=Licor_Datetime;
    x=fullpathflow(end-22:end)
    filename = ['C:\Users\rps207\Documents\SSBjulyCruiseData\Pump Flow Timeseries from measurements\' x ];
    delimiter = ' '; formatSpec = '%f%[^\n\r]'; fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);
    fclose(fileID);
    MeasuredFlow_interp = dataArray{:, 1};
    %smooth flow rate data- as it is a fake interpolation no need to smooth
    smoothflow = MeasuredFlow_interp; %smooth data using filtfilt
    %calculate the lag from  the flow rate
    %volume of tubing in m^3
    tubevol=54*pi*((0.0127/2)^2); %diameter tubing is 1.27cm and 54m long tubing
    tubevolL=tubevol*1000;
    % clearvars filename delimiter formatSpec fileID dataArray ans;
    lagsecs=60*tubevolL./smoothflow(:,1); %lagtime in seconds - add to CTD
    lagsecs(~isfinite(lagsecs))=0;
    respsecs=24; %response time in seconds - add to first CO2 point
end
function [Temperature_probe_a , Temperature_probe_b , Datetime_probe_a, Datetime_probe_b] = getPT100seasonalstudy(fullpathPT100a,fullpathPT100b)

%importpt100 text docs into workspace be selecting folder location
delimiter = ' '; formatSpec = '%s%f%s%s%s%[^\n\r]'; fileID = fopen(fullpathPT100a,'r');fileID2 = fopen(fullpathPT100b,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);
dataArray2 = textscan(fileID2, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);
clearvars fileI fileID2 formatSpec path fnamePT100 fnamePT1002 delimiter; 
Temperature_probe_a_missing_temperatures = dataArray{:, 2};%define probe A Temperature
a=  dataArray{:, 4}; b=dataArray{:, 5};%define probe a time from ymd hms data columsn
Datetime_probe_a_missing_times= strcat(a,b);
Datetime_probe_a_missing_times = datenum(Datetime_probe_a_missing_times,'yyyy-mm-ddHH:MM:SS');
Datetime_probe_a=Datetime_probe_a_missing_times(1):datenum ([ 0 0 0 0 0 1]):Datetime_probe_a_missing_times(end);
Datetime_probe_a=Datetime_probe_a';
clearvars Datetime_seconds_probe_a Datetime_minutes_probe_a Datetime_hours_probe_a Datetime_Time_GMT_probe_a_char  Datetime_Time_GMT_probe_a Datetime_DOM_probe_a Datetime_month_probe_a Datetime_Year_probe_a Datetime_vector_probe_a
% account for the tinker forge python software not updating every second and multiple time points 
B = round(((Datetime_probe_a_missing_times(:) - Datetime_probe_a_missing_times(1)))*24*3600)+1; %convert from datenum to seconds- note (round is probably unnecessary)
[~,IA,~] = unique(B);%find only unique entries in time
A=round(B(1):1:B(end))'; %create a new time vector that has the desired time spacing
ind=find(ismember(A,B(IA))==1); %use ismember and find to find the indexes in A for which there is no value in B
Temperature_probe_a=ones(length(A),1)*nan; %Create a vector of Nans that is the length of A
Temperature_probe_a(ind)=Temperature_probe_a_missing_temperatures(IA); %Insert the Nan values into probe A temperatures
while sum(isnan(Temperature_probe_a))>0  %keep finding nans until the condition >0 is met
    vrow = find(isnan(Temperature_probe_a)); %find all the nans
    Temperature_probe_a(vrow) = Temperature_probe_a(vrow-1); %replace with the cell above
end
clearvars A B a b B_unique IA IC ind vrow Temperature_probe_a_missing_temperatures dataArray Datetime_probe_a_missing_times
Temperature_probe_b_missing_temperatures = dataArray2{:, 2};%define probe B Temperature
c=  dataArray2{:, 4}; d=dataArray2{:, 5};%define probe a time from ymd hms data columsn
Datetime_probe_b_missing_times= strcat(c,d);
Datetime_probe_b_missing_times = datenum(Datetime_probe_b_missing_times,'yyyy-mm-ddHH:MM:SS');
Datetime_probe_b=(Datetime_probe_b_missing_times(1):datenum ([ 0 0 0 0 0 1]):Datetime_probe_b_missing_times(end))';
clearvars Datetime_seconds_probe_a Datetime_minutes_probe_a Datetime_hours_probe_a Datetime_Time_GMT_probe_a_char  Datetime_Time_GMT_probe_a Datetime_DOM_probe_a Datetime_month_probe_a Datetime_Year_probe_a Datetime_vector_probe_a

% account for the tinker forge python software not updating every second and multiple time points 
B = round(((Datetime_probe_b_missing_times(:) - Datetime_probe_b_missing_times(1)))*24*3600)+1; %convert from datenum to seconds- note (round is probably unnecessary)
[~,IA,~] = unique(B);%find only unique entries in time
A=round(B(1):1:B(end))'; %create a new time vector that has the desired time spacing
ind=find(ismember(A,B(IA))==1); %use ismember and find to find the indexes in A for which there is no value in B
Temperature_probe_b=ones(length(A),1)*nan; %Create a vector of Nans that is the length of A
Temperature_probe_b(ind)=Temperature_probe_b_missing_temperatures(IA); %Insert the Nan values into probe A temperatures
while sum(isnan(Temperature_probe_b))>0  %keep finding nans until the condition >0 is met
    vrow = find(isnan(Temperature_probe_b)); %find all the nans
    Temperature_probe_b(vrow) = Temperature_probe_b(vrow-1); %replace with the cell above
end
%Use calibration coefficent to offset PT100 probes
Temperature_probe_a=(Temperature_probe_a*1.0052)-0.22;
Temperature_probe_b=(Temperature_probe_b*1.0089)-0.22;

end

% clearvars A B B_unique IA IC ind vrow Temperature_probe_b_missing_temperatures dataArray2 Datetime_probe_b_missing_times

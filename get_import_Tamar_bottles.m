%get the co2 from helens bottle river data
clear all

[~, ~, raw, dates] = xlsread('C:\Users\rps207\Documents\Seasonal Study Data\tamar_2014_carbon_Helen_pH.xlsx','Sheet2','B4:N16','',@convertSpreadsheetExcelDates);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[1,2,3,4,13]);
raw = raw(:,[6,7,8,9,10,11,12]);
dates = dates(:,5);
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),dates); % Find non-numeric cells
dates(R) = {NaN}; % Replace non-numeric Excel dates with NaN
data = reshape([raw{:}],size(raw));
Tamar_Stationname = cellVectors([1:3 5:13],1);
Tamar_Stationletter = cellVectors([1:3 5:13],2);
Tamar_Lat = cellVectors([1:3 5:13],3);
Tamar_Long = cellVectors([1:3 5:13],4);
Tamar_Date = datetime([dates{[1:3 5:13],1}].', 'ConvertFrom', 'Excel', 'Format', 'MM/dd/yyyy');
Tamar_LatN = data([1:3 5:13],1);
Tamar_LongW = data([1:3 5:13],2);
Tamar_Depth1 = data([1:3 5:13],3);
Tamar_DICumolkg = data([1:3 5:13],4);
Tamar_TAumolkg = data([1:3 5:13],5);
Tamar_Temp1 = data([1:3 5:13],6);
Tamar_Salinity = data([1:3 5:13],7);
Tamar_Notes = cellVectors([1:3 5:13],5);
clearvars data raw dates cellVectors R;

A=[]%pre allocate matrix
for i=1:length(Tamar_TAumolkg)
A(i,:)=CO2SYS(Tamar_TAumolkg(i),Tamar_DICumolkg(i),1,2,Tamar_Salinity(i),Tamar_Temp1(i),25,1,0,10,2,1,4,1);
end
Tamar_fco2_calc_ta_dic=A(:,4);
clearvars i A Tamar_Lat Tamar_Long

save('Data/Tamar_bottles.mat')

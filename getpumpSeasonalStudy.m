function [pumponind, pumpoffind] = getpumpSeasonalStudy(fullpathpump,Licor_Datetime);
s={}; fid=fopen(fullpathpump,'r'); tline = fgetl(fid);
while ischar(tline)
   s=[s;tline];
   tline = fgetl(fid);
end
Pumpon=s(1,:);  pumpDT_on= datenum(Pumpon,'yyyy-mm-dd HH:MM:SS'); 
Pumpoff=s(2,:);pumpDT_off= datenum(Pumpoff,'yyyy-mm-dd HH:MM:SS'); 

%takes a while for pump to start up so add a precautionary 100 seconds
pumpondelay=100;

pumponind=find(Licor_Datetime>(pumpDT_on +(pumpondelay/86400)));
pumponind=pumponind(1);
pumpoffind=find(Licor_Datetime>(pumpDT_off));
pumpoffind=pumpoffind(1);
%clearvars s tline path fname fid clear Pumpon Pumpoff
end
function[undind1, undind2,undind3, undind4] = getunderwaySeasonalStudy(fullpathund,Licor_Datetime,fList,j);

delimiter = '\t'; formatSpec = '%*s%s%[^\n\r]'; fileID = fopen(fullpathund,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false); fclose(fileID);
VarName4 = dataArray{:, 1};
clearvars filename delimiter formatSpec fileID dataArray ans;

VarName4=VarName4([2 3 5 6],:)
VarName4= (datenum(VarName4,'yyyy-mm-dd HH:MM:SS'));



if  fList(j,1:10)==('2016-04-27'); 
undind1=[];
undind2=[];
undind3=VarName4(3,:);
undind4=VarName4(4,:);
elseif fList(j,1:10)==('2016-05-11'); 
undind1=[];
undind2=[];
undind3=VarName4(3,:);
undind4=VarName4(4,:);
elseif fList(j,1:10)==('2016-06-10'); 
undind1=[];
undind2=[];
undind3=VarName4(3,:);
undind4=VarName4(4,:);
elseif fList(j,1:10)==('2016-06-10'); 
undind1=[];
undind2=[];
undind3=VarName4(3,:);
undind4=VarName4(4,:);
elseif fList(j,1:10)==('2016-09-21'); 
undind1=[];
undind2=[];
undind3=VarName4(3,:);
undind4=VarName4(4,:);
elseif fList(j,1:10)==('2016-05-25'); 
undind1=VarName4(1,:);
undind2=VarName4(2,:);
undind3=[];
undind4=[];
elseif fList(j,1:10)==('2016-09-02'); 
undind1=VarName4(1,:);
undind2=VarName4(2,:);
undind3=[];
undind4=[];
elseif fList(j,1:10)==('2016-06-22'); 
undind1=VarName4(1,:);
undind2=VarName4(2,:);
undind3=[];
undind4=[];
else  
undind1=VarName4(1,:);
undind2=VarName4(2,:);
undind3=VarName4(3,:);
undind4=VarName4(4,:);

end



if isempty(undind1)
    
else
undind1=find(Licor_Datetime>(undind1));
undind1=undind1(1);
end

if isempty(undind2)
    
else
undind2=find(Licor_Datetime>(undind2));
undind2=undind2(1);
end

if isempty(undind3)
    
else
undind3=find(Licor_Datetime>(undind3));
undind3=undind3(1);
end

if isempty(undind4)
    
else
undind4=find(Licor_Datetime>(undind4));
undind4=undind4(1);
end



   
   
   
   
   
   
   
   

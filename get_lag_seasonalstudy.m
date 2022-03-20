%use flist 10 from seasonal study as the data to assess ques delay
load('C:\Users\rps207\Documents\MATLAB\CO2 NSOP output analysis\Data\2016-06-10.mat')


%load in data from 2016-06-10 - average of [77;69;91] -mean =79
%using previous calculation can work out spatial offset. assuming 9 knot
%speed of ship
% 48 second =220 m
% 79 second =
% 
% (220/48)*79=360m

XXX = Licor_Datetime(undind3:undind4);   %XXX = replace with your time vector (secs)
YYY = Equilibrator_Temperature(undind3:undind4);   %YYY = replace with your temp sensor 1
ZZZ = SeaTemp_interp(undind3:undind4);   %ZZZ = replace with your temp sensor 2
AAA = underwayT_sea_interp(undind3:undind4);   %ZZZ = replace with your temp sensor 2

plot(XXX,YYY);hold on; plot(XXX,ZZZ,'r'); plot(XXX,AAA,'g');dynamicDateTicks([], [], ' HH:MM'); %plot vectors to check they work

bestLags=[];  %preallocate matrix for speed
int = 300;   %define interval at which to assess lag (secs) eg- 120 means the lag is determined every 120 seconds
numInt = fix((length(XXX)./int))-1; %identify number of intervals within dataset to look at the lag eg if data is 8000 long and int is 100, there are 80 intervals
maxL = 200;        %define maximum lag (secs), theoretically suggest what the lag should be pick an upper measurement
intSt = 1+maxL;     %interval start from where to start the assement of the lag
lags = -maxL:maxL;  %range of lags (from negative to positive), the reason for this is lag may be negative
SIZE=size(XXX);
intend= SIZE(1)-int-max(lags)




for step=intSt:int:intend % using the start point (intst) and interval(int) to find the indexes for the step to work through for lag 
    c = zeros(1,length(lags));  %reset correlation output to zeros
    for ilag = 1:length(lags);      %loop through lags- this works by working out the correlation coefficent for every value between -maxL:maxL
        cmatrix = corrcoef(YYY(step+lags(ilag):step+int+lags(ilag)),ZZZ(step:step+int));     %perform correlation analysis
       c(ilag) = cmatrix(1,2);     %build correlation vector, combining all the correlation stats for each of the lags
    end
    bestLag = lags(c==max(c));      %identify maximum correlation within vector,  this should correspond to the point where the lags are
    figure;subplot(2,1,1); 
    plot(lags,c,'.'); text(10,0.95,sprintf(['Best Lag = ' num2str(bestLag) ' secs']), 'fontsize',16);
    xlabel('Time Lag (mins)', 'fontsize',16); ylabel('Corr Coefficient', 'fontsize',16);
    set(gca,'fontsize',18, 'LineWidth',1.5);
    bestLags = [bestLags ; XXX(step) , bestLag];    %build output array of time and bestlag
    
    subplot(2,1,2); plot(YYY(step:step+int), 'c', 'LineWidth',1.5); hold on; 
plot(ZZZ(step:step+int), 'LineWidth',1.5); 
plot(YYY(step+bestLag:step+bestLag+int),'--k', 'LineWidth',1.5);
legend('YYY','ZZZ','ZZZ+bestLag');
set(gca,'fontsize',18, 'LineWidth',1.5);

end

figure(100000);
plot(bestLags(:,1),bestLags(:,2),'r*')
hold on
plot(Licor_Datetime,lagsecs)
ylim([0 360])
        dynamicDateTicks([], [], ' HH:MM');



mean(bestLags(:,2))


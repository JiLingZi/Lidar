%% Î²¼£±ä»¯Í¼
TimeCap = TimeXCap(:,2075:2085)';
TimeCa = TimeXCa(:,2085:2091)';
TimeNumCap = datenum(TimeCap);
TimeNumCa = datenum(TimeCa);
DenTrailCap = [17 41 72 82 111 106 111 111 90 26 13];
DenTrailCa = [160 366 594 766 773 419 74];
% DenTrailCap = [40 60 34];
% DenTrailCa = [78 141 100];
Profx = 2075;
SEC = 1/(24*60*60);
TimeIntCap = (TimeNumCap(1,:):SEC:TimeNumCap(end,:))';
TimeIntCa = (TimeNumCa(1,:):SEC:TimeNumCa(end,:))';
FitCap = interp1(TimeNumCap,DenTrailCap,TimeIntCap,'linear');
FitCa = interp1(TimeNumCa,DenTrailCa,TimeIntCa,'linear');

% ¼ì²éÍ¼Ïñ
figure('Name','Check')
plot(TimeIntCa,FitCa,'-b','linewidth',1.5)
hold on;
plot(TimeIntCap,FitCap,'--b','linewidth',1.5)
datetick('x','HH:MM:SS')
% ylim([0 600])
ylabel('Density of Trail (cm^{-3})')
grid on
TDstr = 'starttime: '+string(Profx)+' | '+TimeList(:,Profx);
title(TDstr)
xlabel('Time')
datetick('x','HH:MM:SS','keeplimits','keepticks');
legend('Original Ca','Original Ca^+','location','northwest')

%% ÐÞ¸ÄÓë²Ã¼ô

% Ca ×ó

% StdT = round((TimeIntCap(1,:)-TimeIntCa(1,:))./SEC);
% EndT = round((TimeIntCap(end,:)-TimeIntCa(end,:))./SEC);
% 
% IntCap = FitCap(1:end-EndT,:);
% IntCa = FitCa(StdT+1:end,:);
% TimeCut = TimeIntCap(1:end-EndT,:);

% Cap ×ó

% StdT = round((TimeIntCa(1,:)-TimeIntCap(1,:))./SEC);
% EndT = round((TimeIntCa(end,:)-TimeIntCap(end,:))./SEC);
% 
% IntCa = FitCa(1:end-EndT,:);
% IntCap = FitCap(StdT+1:end,:);
% TimeCut = TimeIntCap(StdT+1:end,:);

% ¶¨ÖÆ

StdT = round((TimeIntCa(1,:)-TimeIntCap(1,:))./SEC);
EndT = round((TimeIntCap(end,:)-TimeIntCa(end,:))./SEC);

IntCa = FitCa(1:end,:);
IntCap = FitCap(StdT+1:end-EndT,:);
TimeCut = TimeIntCap(StdT+1:end-EndT,:);

SunCaCap = IntCap+IntCa;
RatioCap = IntCap./SunCaCap;
RatioCa = IntCa./SunCaCap;

RatioSum = zeros(size(RatioCa,1)-1,1);
for jr = 2:size(RatioCa,1)
    RatioSum(jr-1,:) = RatioCa(jr,:)-RatioCa(jr-1,:);
end
RatioCaAVG = mean(RatioSum)
RatioCaSUM = sum(RatioSum)


%% »æÍ¼
figure('Name','Origin','position',[50,100,1450,735])
subplot(2,2,1)
plot(TimeIntCa,FitCa,'-b','linewidth',1.5)
hold on;
plot(TimeIntCap,FitCap,'--b','linewidth',1.5)
datetick('x','HH:MM:SS')
% ylim([0 600])
ylabel('Density of Trail (cm^{-3})')
grid on
TDstr = 'starttime: '+string(Profx)+' | '+TimeList(:,Profx);
title(TDstr)
xlabel('Time')
datetick('x','HH:MM:SS','keeplimits','keepticks');
legend('Original Ca','Original Ca^+','location','northwest')

% figure('Name','Cut')
subplot(2,2,2)
plot(TimeCut,IntCa,'-b','linewidth',1.5)
hold on;
plot(TimeCut,IntCap,'--b','linewidth',1.5)
datetick('x','HH:MM:SS')
% ylim([0 600])
ylabel('Density of Trail (cm^{-3})')
grid on
TDstr = 'starttime: '+string(Profx)+' | '+TimeList(:,Profx);
title(TDstr)
xlabel('Time')
datetick('x','HH:MM:SS','keeplimits','keepticks');
legend('Fit Ca','Fit Ca^+','location','northwest')

% figure('Name','Add')
subplot(2,2,3)
plot(TimeCut,SunCaCap,'-r','linewidth',1.5)
datetick('x','HH:MM:SS')
% ylim([0 600])
ylabel('Density of Trail (cm^{-3})')
grid on
TDstr = 'starttime: '+string(Profx)+' | '+TimeList(:,Profx);
title(TDstr)
xlabel('Time')
datetick('x','HH:MM:SS','keeplimits','keepticks');
legend('Sum Density','location','northwest')

% figure('Name','Ratio')
subplot(2,2,4)
plot(TimeCut,RatioCa.*100,'-b','linewidth',1.5)
hold on;
plot(TimeCut,RatioCap.*100,'--b','linewidth',1.5)
datetick('x','HH:MM:SS')
% ylim([0 600])
ylabel('Ratio of Trail (%)')
grid on
TDstr = 'starttime: '+string(Profx)+' | '+TimeList(:,Profx);
title(TDstr)
xlabel('Time')
datetick('x','HH:MM:SS','keeplimits','keepticks');
legend('Ratio Ca','Ratio Ca^+','location','northwest')

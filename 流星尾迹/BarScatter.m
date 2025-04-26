XlsName ='D:\KYBF\TrailPPT\Datas\�ָ߶�ʱ��ֲ�.xlsx';
Table = readtable(XlsName,'sheet','CapAlt');
AltCeil = table2array(Table(:,4));
AltCeil = AltCeil(~isnan(AltCeil));
DenSum = table2array(Table(:,5));
DenSum = DenSum(~isnan(DenSum));
AlCnt = table2array(Table(:,6));
AlCnt = AlCnt(~isnan(AlCnt));

fig = openfig('D:\File\20240927ppt\MeanCap.fig'); 
lines = findobj(fig, 'Type', 'line'); 
for i = 1
    xData = get(lines(i), 'XData'); 
    yData = get(lines(i), 'YData'); 
end
AltitudeCap = xData';
DenMean = yData';
close;

SumAltCnt = sum(AltCeil.*AlCnt,1);
SumCnt = sum(AlCnt,1);
CoeAlt = SumAltCnt./SumCnt;

SumAltDen = sum(AltCeil.*DenSum,1);
SumDen2 = sum(DenSum,1);
CoeDenAlt = SumAltDen./SumDen2;

CoeAltx = AltitudeCap(246:376,:);
CoeDenx = DenMean(246:376,:);
CoeXxX = CoeDenx.*CoeAltx;
AltDenx = sum(CoeXxX,1)./sum(CoeDenx,1);

% cftool(AltCeil,AlCnt)

GausF = fit(AltCeil,AlCnt,'gauss1');
GausMix = GausF(CoeAltx);
GausMax = max(GausMix);
GausAltidx = find(GausMix == GausMax);
GausAlt = CoeAltx(GausAltidx);

figure('name','��״ͼ')
yyaxis left
bar(AltCeil,AlCnt,'b');
xlabel('Altitude (km)');
ylabel('Number');
title('Number of Ca^+ Meteor Trails');
% xlim([0 20]);
xlim([75 115]);
grid on; hold on;
plot([CoeAlt CoeAlt],[0 100],'--b','linewidth',1.5);
hold on;
plot([AltDenx AltDenx],[0 100],'--r','linewidth',1.5);
yyaxis right
plot(AltitudeCap,DenMean,'-r','linewidth',1.5);
ylim([0 120])
ylabel('Density (cm^{-3})');
xlabel('Altitude (km)');
disp([CoeAlt,AltDenx,GausAlt])
legend('Trail number','Median altitude','Centroid altitude','Mean density')
set(gca,'fontsize',12)

% figure('name','�ܶ�ͼ')
% plot(DenSum,AltCeil,'-b','linewidth',1.5);
% xlabel('Density (cm^{-3})');
% ylabel('Altitude (km)');
% title('Density Distribution of Ca ion Meteor Trails');
% % xlim([0 20]);
% ylim([75 115]);
% grid on;

%% ʱ��
Table2 = readtable(XlsName,'sheet','Fe');
DateStr = char(Table2{:,1});
DateTrail = datetime(DateStr,'InputFormat','yyyyMMdd''T');
TimeSum = table2array(Table2(:,2));
TimeSum = TimeSum(~isnan(TimeSum));
Counts = table2array(Table2(:,3));
Counts = Counts(~isnan(Counts));
Density = table2array(Table2(:,4));
Density = Density(~isnan(Density));
CR = table2array(Table2(:,5));
CR = CR(~isnan(CR));
DR = table2array(Table2(:,6));
DR = DR(~isnan(DR));

figure('name','��ʱ��','position',[100,50,1280,830])
subplot(3,1,1)
bar(DateTrail,TimeSum,'b')
title('�۲�ʱ��')
xlabel('Date');
ylabel('hours');
grid on;
subplot(3,1,2)
yyaxis left;
bar(DateTrail,Counts,'b')
title('β������')
xlabel('Date');
ylabel('Number');
grid on;
yyaxis right;
plot(DateTrail,CR,'-r','linewidth',1.5)
ylabel('Number / h');
subplot(3,1,3)
yyaxis left;
bar(DateTrail,Density,'b')
title('β���ܶ�')
xlabel('Date');
ylabel('Density (cm^{-3})');
grid on;
yyaxis right;
plot(DateTrail,DR,'-r','linewidth',1.5)
ylabel('Density (cm^{-3}/h)');

%% �¶�ɢ��ͼ
startDate = datetime(2023, 10, 1);
endDate = datetime(2024, 9, 1);
stlim = datetime(2023, 9, 15);
edlim = datetime(2024, 9, 15);
Yue = startDate:calmonths(1):endDate;

XlsNameN ='D:\KYBF\TrailPPT\Datas\ʱ��������״ͼ.xlsx';
TableN = readtable(XlsNameN,'sheet','Ca');
AryN = table2array(TableN);
Num = AryN(:,2);%����
Num2 = AryN(:,3);%ʱ��
Num3 = Num./Num2;
NaNto0 = isnan(Num3);
Num3(NaNto0) = 0;
figure('name','�¶�ɢ��-����')
yyaxis left
bar(Yue,Num,0.7,'b','linewidth',1)
grid on;
box on;
xlim([stlim,edlim]);
datetick('x','mmm','keeplimits')
xlabel('Time');
ylabel('Ca Meteor Trail Number')
yyaxis right
plot(Yue,Num2,'--r','linewidth',2)
ylabel('Observation Duration (h)')
title('Ca Meteor Trail Monthly Distribution')
set(gca, 'fontsize',12);

figure('name','�¶�ɢ��-ƽ��ÿСʱ����')
yyaxis left
bar(Yue,Num,0.7,'b','linewidth',1)
grid on;
box on;
xlim([stlim,edlim]);
datetick('x','mmm','keeplimits')
xlabel('Time');
ylabel('Fe Meteor Trail Number')
yyaxis right
plot(Yue,Num3,'--r','linewidth',2)
ylabel('Average Number of Trail per Hour')
title('Na Meteor Trail Monthly Distribution')
set(gca, 'fontsize',12);

% Den = [7179 17177 0 5643 0 53 0 1519];
% figure('name','�¶�ɢ��-�ܶ�')
% bar(Yue,Den,0.7,'linewidth',2)
% grid on;
% box on;
% xlim([stlim,edlim]);
% datetick('x','mmm','keeplimits')
% xlabel('Time');
% ylabel('Fe Meteor Trail Density (cm^{-3})')
% set(gca, 'fontsize',12)




%%

% �� figure �ļ� 
fig = openfig('D:\File\20240927ppt\MeanCap.fig'); 
% ��ȡ figure �����е� line �������ߣ� 
lines = findobj(fig, 'Type', 'line'); 
% ��ȡ figure �����е� bar ������״ͼ�� 
bars = findobj(fig, 'Type', 'bar'); 
% ����ÿ�� line ������ȡ���� 
for i = 1:length(lines) 
    xData = get(lines(i), 'XData'); 
    yData = get(lines(i), 'YData'); 
    % ��ӡ���ݻ������������ 
    fprintf('Line %d:\n', i); 
    disp('X Data:'); 
    disp(xData); 
    disp('Y Data:'); 
    disp(yData); 
end

for i = 1
    xData = get(lines(i), 'XData'); 
    yData = get(lines(i), 'YData'); 
    % ��ӡ���ݻ������������ 
    fprintf('Line %d:\n', i); 
    disp('X Data:'); 
    disp(xData); 
    disp('Y Data:'); 
    disp(yData); 
end













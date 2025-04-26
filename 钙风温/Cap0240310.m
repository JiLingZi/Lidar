% CaIon�ܶ�
% ��ȡ����������,��ԭ�ӡ���������4��ͨ��
ReadDate = '20231105';
fldstart = ['F:\RawData\ZWDATA\MOHEnew\Cap\',ReadDate,'\'];
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1;% Ϊ��ͳһ�ļ�������һ
% ��һ����Ϊ�Ǹ��ƻ�����������һ���ļ������²������ĸ�ͨ���ļ�����ͬ
% ��ʽ��Ҫ�ٸ��ˣ�̽��߶�����Ҳ��Ҫ�ٸ��ˣ�����ʵ���Զ���������������

% bin = 26432;
bin = 16384;
% bin = 13700;
lveHang = 20;
tt = 9;
ttstr = 'yyyy-MM-dd''T''HH:mm:ss.SSS';

DenPh1 = zeros(bin, num_files);
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
for j = 1:num_files
    filename = fullfile(folder,files(j).name);
    fid = fopen(filename);
    PhData1 = textscan(fid, '%f %f','HeaderLines',lveHang);
    fclose(fid);
    DenPh1(:,j) = PhData1{1,2}(1:bin,:);
end


% DenPh2 = zeros(bin, num_files);
% folder = [fldstart,'CH2\'];
% files = dir(fullfile(folder, '*.txt'));
% for j = 1:num_files
%     filename = fullfile(folder,files(j).name);
%     fid = fopen(filename);
%     PhData2 = textscan(fid, '%f %f','HeaderLines',lveHang);
%     fclose(fid);
%     DenPh2(:,j) = PhData2{1,2}(1:bin,:);
% end
% 
% DenPh3 = zeros(bin, num_files);
% folder = [fldstart,'CH3\'];
% files = dir(fullfile(folder, '*.txt'));
% for j = 1:num_files
%     filename = fullfile(folder,files(j).name);
%     fid = fopen(filename);
%     PhData3 = textscan(fid, '%f %f','HeaderLines',lveHang);
%     fclose(fid);
%     DenPh3(:,j) = PhData3{1,2}(1:bin,:);
% end
% 
% DenPh4 = zeros(bin, num_files);
% folder = [fldstart,'CH4\'];
% files = dir(fullfile(folder, '*.txt'));
% for j = 1:num_files
%     filename = fullfile(folder,files(j).name);
%     fid = fopen(filename);
%     PhData4 = textscan(fid, '%f %f','HeaderLines',lveHang);
%     fclose(fid);
%     DenPh4(:,j) = PhData4{1,2}(1:bin,:);
% end

% DenPh = DenPh1+DenPh2+DenPh3+DenPh4;
DenPh = DenPh1;

% ��ȡʱ������
TimeValues = cell(1,num_files);
for j = 1:num_files
    filename = fullfile(folder, files(j).name);
    fid = fopen(filename);
    TXT = textscan(fid, '%s', 'Delimiter', '\n');
    RowTime = TXT{1}{tt};
    TimeData = strsplit(RowTime);
    fclose(fid);
    TimeValues(:,j) = cellstr(TimeData{2:end});
end

% �ϲ�ʱ���ļ���
% jNum = 54;% 1h
jNum = 9;% 10min

% �ϲ��ռ�ֱ���
iNum = 3;% 100m
% iNum = 32;% 1km
% iNum = 65;% 2km
% iNum = 163;% 5km
% iNum = 325;% 10km
% iNum = 651;% 20km

TimeList = cell(1,floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
end

% ��ȡԭʼ�߶Ⱦ���
height_num_origin = PhData1{1,1}(1:bin,:);

% �ϲ��У�3�кϲ���96��
PhiSum = zeros(floor(bin/iNum),num_files);
for i = 1:floor(bin/iNum)
    PhiSum(i,:) = sum(DenPh(iNum*(i-1)+1:iNum*i,:),1);
end

% �ϲ��У�15�кϲ���15����
PhjSum = zeros(floor(bin/iNum),floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
    PhjSum(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
end

% ����120-150��������������߶�
Altitude = height_num_origin(1:floor(bin/iNum));
for i = 1:floor(bin/iNum)
    Altitude(i,:) = height_num_origin(i*iNum-floor(iNum/2),:);
end
% Altitude = height_num_origin(1:floor(26432/iNum))*iNum;
KM30 = size(Altitude(Altitude<30),1)+1;
KM75 = size(Altitude(Altitude<75),1)+1;
KM80 = size(Altitude(Altitude<80),1)+1;
KM90 = size(Altitude(Altitude<90),1)+1;
KM100 = size(Altitude(Altitude<100),1)+1;
KM110 = size(Altitude(Altitude<110),1)+1;
KM115 = size(Altitude(Altitude<115),1)+1;
KM120 = size(Altitude(Altitude<120),1)+1;
KM125 = size(Altitude(Altitude<125),1)+1;
KM150 = size(Altitude(Altitude<150),1)+1;
KM200 = size(Altitude(Altitude<200),1)+1;
KM300 = size(Altitude(Altitude<300),1)+1;
KM350 = size(Altitude(Altitude<350),1)+1;
KM400 = size(Altitude(Altitude<400),1)+1;

% ��һ�����������
CheckH = 50;% ����������һ���߶ȣ�KM
KMCheck = size(Altitude(Altitude<CheckH),1);
CheckPh = PhjSum./PhjSum(KMCheck,:);
figure('Name','���������','position',[10 100 1350 500])
for j = 1:size(CheckPh,2)
    plot(CheckPh(:,j)+j,Altitude,'-','linewidth',1.5);
    hold on;
end
xlim([0 31]);
ylim([75 115]);
% ylim([150 350]);
xInterval = 1;  % x��̶ȼ��

% % ����������ͼ
% JCheck = 3;
% plot(Altitude,PhjSum(:,JCheck),'-k','linewidth',1.5);
% grid on;
% % xlim([75 140]);
% xlim([200 300]);
% xlabel('Altitude (km)');
% ylabel('Photon Counts');
% TitleStr = string(TimeList(:,JCheck)) + " " + 'UTC';
% title(TitleStr)
% % ylim([75 115]);

mean(max(PhjSum(KM80:KM110,:)))
noise = mean(mean(PhjSum(KM350:KM400,:)))

% % �������
% peak = 12887735;
% noise = mean(PhjSum(KM120:KM125,JCheck),1)
% Err = sqrt(peak)./(peak-noise)

%{
temporal resolution: 1 h
spatial resolution: 5 km
peak photon counts: 17077
%}

% ����
Noise = mean(PhjSum(KM350:KM400,:),1);
SumPh = PhjSum - Noise;

SNR = SumPh./Noise;



STD = std(PhjSum(KM350:KM400,:),1);
DetcRange = Noise + STD;



% ��ȡ��������ɢ��������Ч����ɢ�����
ScatterTable = readtable('RayEffScatter.txt');
Scatter = table2array(ScatterTable(:,2:3));
RayScatter = Scatter(3,1);
EffScatter = Scatter(3,2);

% ��ȡ30 km ������ģ���¶����ܶ�
[TemRay,DenRay] = atmosnrlmsise00(30000,53.3,122.7,2023,201,17);

% �����ܶ�
Z = Altitude;                       % �߶Ⱦ���
ZR = 30;                            % �ο��߶�
SigmaRay = RayScatter;              % ��������ɢ�����
SigmaEff = EffScatter;              % ��Ч����ɢ�����
NZ = PhjSum;                        % ����������
NB = Noise;                         % ��������
NZR = PhjSum(KM30,:);               % �ο��߶ȴ�������
NumberRay = sum(DenRay)*1e-6;       % �ο��߶ȴ�����ģ�����ܶ�
NumberZ = (((Z.^2).*SigmaRay.*(NZ-NB))./((ZR^2)*SigmaEff*(NZR-NB)))*NumberRay;

% % �ܶ����ߴּ��
% NaNJ = [27:30];
% NumberZ(:,NaNJ) = NaN;

figure(2)
Density = NumberZ(KM75:KM125,:);
Density(Density<0) = NaN;
TimeX = datetime(TimeList, 'InputFormat', ttstr);
TimeY = 1:size(TimeX,2);
hours = hour(TimeX);
minutes = minute(TimeX);
Time = (((hours(1,:)+8)*3600+minutes(1,:)*60)+((((hours(end,:)+8)*3600+minutes(end,:)*60)-((hours(1,:)+8)*3600+minutes(1,:)*60)).*(TimeY./TimeY(end))))./3600;
AltitudeY = Altitude(KM75:KM125,:);
plot(AltitudeY,Density);
xlim([75 125]);

% ���ܶ�ʱ��ͼ
figure(3)
[Map,Line] = contourf(Time,AltitudeY,Density,30);
set(Line,'LineColor','none')
set(gca,'FontSize',12,'FontName','Times New Roman');
title(fldstart);
% set(gca, 'ColorScale', 'log');
colormap(jet); 
colorbar;
ylabel(colorbar,'Density (cm^{-3})');
xlabel('Local Time');
ylabel('Altitude (km)');
ylim([75 125]);
% xlim([18.2 23])

figure(4)
Density = NumberZ(KM150:KM400,:);
Density(Density<0) = NaN;
TimeX = datetime(TimeList, 'InputFormat', ttstr);
TimeY = 1:size(TimeX,2);
hours = hour(TimeX);
minutes = minute(TimeX);
Time = (((hours(1,:)+8)*3600+minutes(1,:)*60)+((((hours(end,:)+8)*3600+minutes(end,:)*60)-((hours(1,:)+8)*3600+minutes(1,:)*60)).*(TimeY./TimeY(end))))./3600;
AltitudeY = Altitude(KM150:KM400,:);
plot(AltitudeY,Density);
xlim([115 400]);

% ���ܶ�ʱ��ͼ
figure(5)
[Map,Line] = contourf(Time,AltitudeY,Density,30);
set(Line,'LineColor','none')
set(gca,'FontSize',12,'FontName','Times New Roman');
title(fldstart);
set(gca, 'ColorScale', 'log');
colormap(jet); 
colorbar;
ylabel(colorbar,'Density (cm^{-3})');
xlabel('Local Time');
ylabel('Altitude (km)');
ylim([150 400]);
% xlim([18.2 23])


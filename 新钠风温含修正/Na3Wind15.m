%% 读取卫星数据
% load '卫星数据\Altitude.mat';
% load '卫星数据\Date.mat';
% load '卫星数据\Latitude.mat';
% load '卫星数据\Longitude.mat';
% load '卫星数据\Temperature.mat';
% load '卫星数据\Time.mat';
% 
% SbAlt = double(Saber_Altitude);
% SbDate = Saber_date';
% SbLat = double(Saber_Latitude);
% SbLon = double(Saber_Longtitude);
% SbTem = double(Saber_Temperture);
% SbTime = saber_yymmdd;
% 
% RunRead = '卫星数据读取OK了'

%% 读取光子数矩阵
ReadDate = '20240120';
fldstart = ['F:\RawData\ZWDATA\MOHEnew\Na\',ReadDate,'\'];
folder = [fldstart,'Na\'];
files = dir(fullfile(folder, '*.dat'));
num_files = length(files)-1;% 为了统一文件数，减一

V_F_0 = zeros(8192, num_files);
V_F_R = zeros(8192, num_files);
V_F_L = zeros(8192, num_files);
N_F_0 = zeros(8192, num_files);
N_F_R = zeros(8192, num_files);
N_F_L = zeros(8192, num_files);
E_F_0 = zeros(8192, num_files);
E_F_R = zeros(8192, num_files);
E_F_L = zeros(8192, num_files);
for j = 1:num_files
    filename = fullfile(folder, files(j).name);
    Na_table = readtable(filename);
    Na_data = table2array(Na_table);
    V_F_0(:,j) = Na_data(1:8192,2);
    V_F_R(:,j) = Na_data(1:8192,3);
    V_F_L(:,j) = Na_data(1:8192,4);
    N_F_0(:,j) = Na_data(1:8192,5);
    N_F_R(:,j) = Na_data(1:8192,6);
    N_F_L(:,j) = Na_data(1:8192,7);
    E_F_0(:,j) = Na_data(1:8192,8);
    E_F_R(:,j) = Na_data(1:8192,9);
    E_F_L(:,j) = Na_data(1:8192,10);
end

% 读取时间序列
TimeValues = cell(1,num_files);
for j = 1:num_files
    filename = fullfile(folder, files(j).name);
    fid = fopen(filename);
    TXT = textscan(fid, '%s', 'Delimiter', '\n');
    RowTime = TXT{1}{4};
    TimeData = strsplit(RowTime);
    fclose(fid);
    TimeValues(:,j) = cellstr(TimeData{4});
end
RunRead = '数据读取OK了'

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

%% 筛选数据
% NoCh = [];
% V_F_0(1140,NoCh);

%% 合并时间文件数
% jNum = 1;
% iNum = 16;
% 
% TimeList = cell(1,floor(num_files/jNum));
% for j = 1:floor(num_files/jNum)
% 	TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
% end
% 
% TimeX = datetime(TimeList, 'InputFormat', ':yyyy-MM-dd''T''HH:mm:ss.SSS');
% 
% % 获取原始高度矩阵
% height_num_origin = Na_data(1:8192,1);
% 
% % 合并行
% VF_0i = zeros(floor(8192/iNum),num_files);
% VF_Ri = zeros(floor(8192/iNum),num_files);
% VF_Li = zeros(floor(8192/iNum),num_files);
% NF_0i = zeros(floor(8192/iNum),num_files);
% NF_Ri = zeros(floor(8192/iNum),num_files);
% NF_Li = zeros(floor(8192/iNum),num_files);
% EF_0i = zeros(floor(8192/iNum),num_files);
% EF_Ri = zeros(floor(8192/iNum),num_files);
% EF_Li = zeros(floor(8192/iNum),num_files);
% for i = 1:floor(8192/iNum)
%     VF_0i(i,:) = sum(V_F_0(iNum*(i-1)+1:iNum*i,:),1);
%     VF_Ri(i,:) = sum(V_F_R(iNum*(i-1)+1:iNum*i,:),1);
%     VF_Li(i,:) = sum(V_F_L(iNum*(i-1)+1:iNum*i,:),1);
%     NF_0i(i,:) = sum(N_F_0(iNum*(i-1)+1:iNum*i,:),1);
%     NF_Ri(i,:) = sum(N_F_R(iNum*(i-1)+1:iNum*i,:),1);
%     NF_Li(i,:) = sum(N_F_L(iNum*(i-1)+1:iNum*i,:),1);
%     EF_0i(i,:) = sum(E_F_0(iNum*(i-1)+1:iNum*i,:),1);
%     EF_Ri(i,:) = sum(E_F_R(iNum*(i-1)+1:iNum*i,:),1);
%     EF_Li(i,:) = sum(E_F_L(iNum*(i-1)+1:iNum*i,:),1);
% end
% 
% % 合并列
% VF_0 = zeros(floor(8192/iNum),floor(num_files/jNum));
% VF_R = zeros(floor(8192/iNum),floor(num_files/jNum));
% VF_L = zeros(floor(8192/iNum),floor(num_files/jNum));
% NF_0 = zeros(floor(8192/iNum),floor(num_files/jNum));
% NF_R = zeros(floor(8192/iNum),floor(num_files/jNum));
% NF_L = zeros(floor(8192/iNum),floor(num_files/jNum));
% EF_0 = zeros(floor(8192/iNum),floor(num_files/jNum));
% EF_R = zeros(floor(8192/iNum),floor(num_files/jNum));
% EF_L = zeros(floor(8192/iNum),floor(num_files/jNum));
% for j = 1:floor(num_files/jNum)
%     VF_0(:,j) = sum(VF_0i(:,jNum*(j-1)+1:jNum*j),2);
%     VF_R(:,j) = sum(VF_Ri(:,jNum*(j-1)+1:jNum*j),2);
%     VF_L(:,j) = sum(VF_Li(:,jNum*(j-1)+1:jNum*j),2);
%     NF_0(:,j) = sum(NF_0i(:,jNum*(j-1)+1:jNum*j),2);
%     NF_R(:,j) = sum(NF_Ri(:,jNum*(j-1)+1:jNum*j),2);
%     NF_L(:,j) = sum(NF_Li(:,jNum*(j-1)+1:jNum*j),2);
%     EF_0(:,j) = sum(EF_0i(:,jNum*(j-1)+1:jNum*j),2);
%     EF_R(:,j) = sum(EF_Ri(:,jNum*(j-1)+1:jNum*j),2);
%     EF_L(:,j) = sum(EF_Li(:,jNum*(j-1)+1:jNum*j),2);
% end
% 
% % 噪声120-150，计算光子数，高度
% Altitude = height_num_origin(1:floor(8192/iNum))*iNum;
% KM30 = size(Altitude(Altitude<30),1)+1;
% KM35 = size(Altitude(Altitude<35),1)+1;
% KM40 = size(Altitude(Altitude<40),1)+1;
% KM75 = size(Altitude(Altitude<75),1)+1;
% KM80 = size(Altitude(Altitude<80),1)+1;
% KM85 = size(Altitude(Altitude<85),1)+1;
% KM90 = size(Altitude(Altitude<90),1)+1;
% KM95 = size(Altitude(Altitude<95),1)+1;
% KM98 = size(Altitude(Altitude<98),1)+1;
% KM100 = size(Altitude(Altitude<100),1)+1;
% KM105 = size(Altitude(Altitude<105),1)+1;
% KM115 = size(Altitude(Altitude<115),1)+1;
% KM120 = size(Altitude(Altitude<120),1)+1;
% KM150 = size(Altitude(Altitude<150),1)+1;
% Alt30 = Altitude*(sqrt(3)/2);
% KMx30 = size(Alt30(Alt30<30),1)+1;
% KMx35 = size(Alt30(Alt30<35),1)+1;
% KMx40 = size(Alt30(Alt30<40),1)+1;
% KMx75 = size(Alt30(Alt30<75),1)+1;
% KMx80 = size(Alt30(Alt30<80),1)+1;
% KMx85 = size(Alt30(Alt30<85),1)+1;
% KMx90 = size(Alt30(Alt30<90),1)+1;
% KMx95 = size(Alt30(Alt30<95),1)+1;
% KMx98 = size(Alt30(Alt30<98),1)+1;
% KMx100 = size(Alt30(Alt30<100),1)+1;
% KMx105 = size(Alt30(Alt30<105),1)+1;
% KMx115 = size(Alt30(Alt30<115),1)+1;
% KMx120 = size(Alt30(Alt30<120),1)+1;
% KMx150 = size(Alt30(Alt30<150),1)+1;
% 
% % 去噪声
% NoiseV0 = mean(VF_0(KM120:KM150,:),1);
% NoiseVR = mean(VF_R(KM120:KM150,:),1);
% NoiseVL = mean(VF_L(KM120:KM150,:),1);
% NoiseN0 = mean(NF_0(KM120:KM150,:),1);
% NoiseNR = mean(NF_R(KM120:KM150,:),1);
% NoiseNL = mean(NF_L(KM120:KM150,:),1);
% NoiseE0 = mean(EF_0(KM120:KM150,:),1);
% NoiseER = mean(EF_R(KM120:KM150,:),1);
% NoiseEL = mean(EF_L(KM120:KM150,:),1);
% 
% VF0 = VF_0 - NoiseV0;
% VFR = VF_R - NoiseVR;
% VFL = VF_L - NoiseVL;
% NF0 = NF_0 - NoiseN0;
% NFR = NF_R - NoiseNR;
% NFL = NF_L - NoiseNR;
% EF0 = EF_0 - NoiseE0;
% EFR = EF_R - NoiseER;
% EFL = EF_L - NoiseEL;
% 
% figure('name','交替噪声检查','position',[400 50 800 750])
% subplot(3,1,1)
% plot(TimeX,NoiseV0, '-k', 'linewidth',1.)
% hold on;
% plot(TimeX,NoiseVR, '-b', 'linewidth',1.)
% hold on;
% plot(TimeX,NoiseVL, '-r', 'linewidth',1.)
% grid on;
% datetick('x','HH:MM')
% % xlim([datetime(2024, 1, 19, 22, 30, 0), datetime(2024, 1, 19, 23, 0, 0)]);
% xlabel('Time (UT)')
% ylabel('Vertical Noise')
% legend('\nu_0','\nu_+','\nu_-')
% subplot(3,1,2)
% plot(TimeX,NoiseN0, '-k', 'linewidth',1.)
% hold on;
% plot(TimeX,NoiseNR, '-b', 'linewidth',1.)
% hold on;
% plot(TimeX,NoiseNL, '-r', 'linewidth',1.)
% grid on;
% datetick('x','HH:MM')
% % xlim([datetime(2024, 1, 19, 22, 30, 0), datetime(2024, 1, 19, 23, 0, 0)]);
% xlabel('Time (UT)')
% ylabel('North Noise')
% legend('\nu_0','\nu_+','\nu_-')
% subplot(3,1,3)
% plot(TimeX,NoiseE0, '-k', 'linewidth',1.)
% hold on;
% plot(TimeX,NoiseER, '-b', 'linewidth',1.)
% hold on;
% plot(TimeX,NoiseEL, '-r', 'linewidth',1.)
% grid on;
% datetick('x','HH:MM')
% % xlim([datetime(2024, 1, 19, 22, 30, 0), datetime(2024, 1, 19, 23, 0, 0)]);
% xlabel('Time (UT)')
% ylabel('East Noise')
% legend('\nu_0','\nu_+','\nu_-')
% 
% 
% 
% 
% RunRead = '合并去噪OK了'
% 
% 
% % 归一化
% V0 = VF0./VF0(KM35,:);
% VR = VFR./VFR(KM35,:);
% VL = VFL./VFL(KM35,:);
% N0 = NF0./NF0(KM35,:);
% NR = NFR./NFR(KM35,:);
% NL = NFL./NFL(KM35,:);
% E0 = EF0./EF0(KM35,:);
% ER = EFR./EFR(KM35,:);
% EL = EFL./EFL(KM35,:);
% 
% RunRead = '归一化OK了'

%% 三频检查
% figure('name','V三频检查','position',[600 50 350 800])
% pf = size(V0,2)-6;
% % pf = 300;
% plot(E0(:,pf), Altitude, '-k', 'linewidth',1.5)
% hold on
% plot(ER(:,pf), Altitude, '-b', 'linewidth',1.5)
% hold on
% plot(EL(:,pf), Altitude, '-r', 'linewidth',1.5)
% grid on
% % set(gca, 'XScale', 'log');
% ylim([30 125])
% % ttstr = string(pf) + ' | ' + TimeList(pf);
% Tits = 'Vertical ';
% ttstr = "Mohe Na " + TimeList(pf);
% title('(c) East')
% set(gca,'FontSize',12.5);
% % ylim([3000 7000])
% xlabel('Photon Counts')
% ylabel('Altitude (km)')
% legend('\nu_{0}','\nu_{+}','\nu_{-}')

%% 归一化光子数检查
% MaxBh = 1000;
% YBc = 0.2;
% CheckV0 = NF0./MaxBh;
% CheckVR = NFR./MaxBh;
% CheckVL = NFL./MaxBh;
% figure('Name','V光子数检查','position',[100 300 1350 500])
% for j = 1:size(CheckV0,2)
%     plot(CheckV0(:,j)*YBc+j,Altitude.*(sqrt(3)/2),'-k','linewidth',1.5);
%     hold on;
%     plot(CheckVR(:,j)*YBc+j,Altitude.*(sqrt(3)/2),'-b','linewidth',1.5);
%     hold on;
%     plot(CheckVL(:,j)*YBc+j,Altitude.*(sqrt(3)/2),'-r','linewidth',1.5);
%     hold on;
%     grid on
% end
% set(gca,'FontSize',12,'FontName','Times New Roman');
% tstr = "Na Photon Profiles "+ReadDate;
% title(tstr);
% xlabel('Local Time');
% ylabel('Altitude (km)');
% xlim([-1 11]);
% ylim([40 125]);

%% 分割时间
% Tfenge = 66;
% TT = TimeList(Tfenge)

% Ratio0 = (VF_0(:,63)+VF_0(:,64)+VF_0(:,65))./(VF_0(:,67)+VF_0(:,68)+VF_0(:,69));
% RatioR = (VF_R(:,63)+VF_R(:,64)+VF_R(:,65))./(VF_R(:,67)+VF_R(:,68)+VF_R(:,69));
% RatioL = (VF_L(:,63)+VF_L(:,64)+VF_L(:,65))./(VF_L(:,67)+VF_L(:,68)+VF_L(:,69));

% Ratio0 = (NF_0(:,63)+NF_0(:,64)+NF_0(:,65))./(NF_0(:,67)+NF_0(:,68)+NF_0(:,69));
% RatioR = (NF_R(:,63)+NF_R(:,64)+NF_R(:,65))./(NF_R(:,67)+NF_R(:,68)+NF_R(:,69));
% RatioL = (NF_L(:,63)+NF_L(:,64)+NF_L(:,65))./(NF_L(:,67)+NF_L(:,68)+NF_L(:,69));

% Ratio0 = (EF_0(:,63)+EF_0(:,64)+EF_0(:,65))./(EF_0(:,67)+EF_0(:,68)+EF_0(:,69));
% RatioR = (EF_R(:,63)+EF_R(:,64)+EF_R(:,65))./(EF_R(:,67)+EF_R(:,68)+EF_R(:,69));
% RatioL = (EF_L(:,63)+EF_L(:,64)+EF_L(:,65))./(EF_L(:,67)+EF_L(:,68)+EF_L(:,69));
% 
% 
% figure('name','Ratio3','position',[200 50 1050 800])
% subplot(1,3,1)
% plot(Ratio0(KM35:KM40,:),Altitude(KM35:KM40,:),'-k','linewidth',1.5)
% hold on;
% plot(RatioR(KM35:KM40,:),Altitude(KM35:KM40,:),'-b','linewidth',1.5)
% hold on;
% plot(RatioL(KM35:KM40,:),Altitude(KM35:KM40,:),'-r','linewidth',1.5)
% grid on;
% title('Ratio of Ray')
% xlabel('Ratio')
% ylabel('Altitude (km)')
% % legend('\nu_0','\nu_+','\nu_-')
% set(gca,'fontname','times new roman','fontsize',15)
% 
% subplot(1,3,2)
% plot(Ratio0(KM85:KM95,:),Altitude(KM85:KM95,:),'-k','linewidth',1.5)
% hold on;
% plot(RatioR(KM85:KM95,:),Altitude(KM85:KM95,:),'-b','linewidth',1.5)
% hold on;
% plot(RatioL(KM85:KM95,:),Altitude(KM85:KM95,:),'-r','linewidth',1.5)
% grid on;
% title('Ratio of Na')
% xlabel('Ratio')
% ylabel('Altitude (km)')
% % legend('\nu_0','\nu_+','\nu_-')
% set(gca,'fontname','times new roman','fontsize',15)
% 
% subplot(1,3,3)
% plot(Ratio0(KM120:KM150,:),Altitude(KM120:KM150,:),'-k','linewidth',1.5)
% hold on;
% plot(RatioR(KM120:KM150,:),Altitude(KM120:KM150,:),'-b','linewidth',1.5)
% hold on;
% plot(RatioL(KM120:KM150,:),Altitude(KM120:KM150,:),'-r','linewidth',1.5)
% grid on;
% title('Ratio of Noise')
% xlabel('Ratio')
% ylabel('Altitude (km)')
% legend('\nu_0','\nu_+','\nu_-')
% set(gca,'fontname','times new roman','fontsize',15)





%% 滤光器系数,Ray选区
% LN1 = Tfenge-4;LN2 = Tfenge-3;LN3 = Tfenge-2;
% RN1 = Tfenge+2;RN2 = Tfenge+3;RN3 = Tfenge+4;
% 
% Ray10 = VF0(:,LN1)+VF0(:,LN2)+VF0(:,LN3);
% Ray20 = VF0(:,RN1)+VF0(:,RN2)+VF0(:,RN3);
% Ray1R = VFR(:,LN1)+VFR(:,LN2)+VFR(:,LN3);
% Ray2R = VFR(:,RN1)+VFR(:,RN2)+VFR(:,RN3);
% Ray1L = VFL(:,LN1)+VFL(:,LN2)+VFL(:,LN3);
% Ray2L = VFL(:,RN1)+VFL(:,RN2)+VFL(:,RN3);
% NaPh10 = VF0(KM85,LN1)+VF0(KM85,LN2)+VF0(KM85,LN3);
% NaPh20 = VF0(KM85,RN1)+VF0(KM85,RN2)+VF0(KM85,RN3);
% NaPh1R = VFR(KM85,LN1)+VFR(KM85,LN2)+VFR(KM85,LN3);
% NaPh2R = VFR(KM85,RN1)+VFR(KM85,RN2)+VFR(KM85,RN3);
% NaPh1L = VFL(KM85,LN1)+VFL(KM85,LN2)+VFL(KM85,LN3);
% NaPh2L = VFL(KM85,RN1)+VFL(KM85,RN2)+VFL(KM85,RN3);
% 
% 
% Tray0 = Ray10./Ray20;
% Tna0 = NaPh10./NaPh20;
% Rf0 = Tray0./Tna0
% TrayR = Ray1R./Ray2R;
% TnaR = NaPh1R./NaPh2R;
% RfR = TrayR./TnaR
% TrayL = Ray1L./Ray2L;
% TnaL = NaPh1L./NaPh2L;
% RfL = TrayL./TnaL
% 
% RunRead = '系数计算OK了'
% 
% plot([Rf0,RfR,RfL],Altitude,'linewidth',2);
% ylim([25 45]);
% xlabel('R_{\nu_i} i = (0,R,L)')
% ylabel('Altitude (km)')
% title('Na @85 km peak, Ray [25,45] km')
% grid on;
% set(gca,'fontname','times new roman')
% legend('Rf0','RfR','RfL');


%% 滤光器系数,Na选区

% Ray10 = VF0(KM35,LN1)+VF0(KM35,LN2)+VF0(KM35,LN3);
% Ray20 = VF0(KM35,RN1)+VF0(KM35,RN2)+VF0(KM35,RN3);
% Ray1R = VFR(KM35,LN1)+VFR(KM35,LN2)+VFR(KM35,LN3);
% Ray2R = VFR(KM35,RN1)+VFR(KM35,RN2)+VFR(KM35,RN3);
% Ray1L = VFL(KM35,LN1)+VFL(KM35,LN2)+VFL(KM35,LN3);
% Ray2L = VFL(KM35,RN1)+VFL(KM35,RN2)+VFL(KM35,RN3);
% NaPh10 = VF0(:,LN1)+VF0(:,LN2)+VF0(:,LN3);
% NaPh20 = VF0(:,RN1)+VF0(:,RN2)+VF0(:,RN3);
% NaPh1R = VFR(:,LN1)+VFR(:,LN2)+VFR(:,LN3);
% NaPh2R = VFR(:,RN1)+VFR(:,RN2)+VFR(:,RN3);
% NaPh1L = VFL(:,LN1)+VFL(:,LN2)+VFL(:,LN3);
% NaPh2L = VFL(:,RN1)+VFL(:,RN2)+VFL(:,RN3);
% 
% 
% Tray0 = Ray10./Ray20;
% Tna0 = NaPh10./NaPh20;
% Rf0 = Tray0./Tna0
% TrayR = Ray1R./Ray2R;
% TnaR = NaPh1R./NaPh2R;
% RfR = TrayR./TnaR
% TrayL = Ray1L./Ray2L;
% TnaL = NaPh1L./NaPh2L;
% RfL = TrayL./TnaL
% 
% RunRead = '系数计算OK了'
% 
% plot([Rf0,RfR,RfL],Altitude,'linewidth',2);
% ylim([80 105]);
% xlabel('R_{\nu_i} i = (0,R,L)')
% ylabel('Altitude (km)')
% title('Ray @35 km, Na [80,105] km')
% grid on;
% set(gca,'fontname','times new roman')
% legend('Rf0','RfR','RfL');


%% 标定系数计算
% Ray10 = mean(VF0(KM35:KM40,LN1)+VF0(KM35:KM40,LN2)+VF0(KM35:KM40,LN3),1);
% Ray20 = mean(VF0(KM35:KM40,RN1)+VF0(KM35:KM40,RN2)+VF0(KM35:KM40,RN3),1);
% Ray1R = mean(VFR(KM35:KM40,LN1)+VFR(KM35:KM40,LN2)+VFR(KM35:KM40,LN3),1);
% Ray2R = mean(VFR(KM35:KM40,RN1)+VFR(KM35:KM40,RN2)+VFR(KM35:KM40,RN3),1);
% Ray1L = mean(VFL(KM35:KM40,LN1)+VFL(KM35:KM40,LN2)+VFL(KM35:KM40,LN3),1);
% Ray2L = mean(VFL(KM35:KM40,RN1)+VFL(KM35:KM40,RN2)+VFL(KM35:KM40,RN3),1);
% NaPh10 = mean(VF0(KM85:KM95,LN1)+VF0(KM85:KM95,LN2)+VF0(KM85:KM95,LN3),1);
% NaPh20 = mean(VF0(KM85:KM95,RN1)+VF0(KM85:KM95,RN2)+VF0(KM85:KM95,RN3),1);
% NaPh1R = mean(VFR(KM85:KM95,LN1)+VFR(KM85:KM95,LN2)+VFR(KM85:KM95,LN3),1);
% NaPh2R = mean(VFR(KM85:KM95,RN1)+VFR(KM85:KM95,RN2)+VFR(KM85:KM95,RN3),1);
% NaPh1L = mean(VFL(KM85:KM95,LN1)+VFL(KM85:KM95,LN2)+VFL(KM85:KM95,LN3),1);
% NaPh2L = mean(VFL(KM85:KM95,RN1)+VFL(KM85:KM95,RN2)+VFL(KM85:KM95,RN3),1);
% 
% Tray0 = Ray10./Ray20;
% Tna0 = NaPh10./NaPh20;
% Rf0 = Tray0./Tna0
% TrayR = Ray1R./Ray2R;
% TnaR = NaPh1R./NaPh2R;
% RfR = TrayR./TnaR
% TrayL = Ray1L./Ray2L;
% TnaL = NaPh1L./NaPh2L;
% RfL = TrayL./TnaL

check = aa;

%% 标定系数
Rf0 = 1.1026;
RfL = 1.2123;
RfR = 0.9621;

RunRead = '系数标定OK了'

% 合并时间文件数
jNum = 30;
iNum = 32;
iNumx = 37;

TimeList = cell(1,floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
	TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
end

TimeX = datetime(TimeList, 'InputFormat', ':yyyy-MM-dd''T''HH:mm:ss.SSS');

% 获取原始高度矩阵
height_num_origin = Na_data(1:8192,1);

% 合并行
VF_0i = zeros(floor(8192/iNum),num_files);
VF_Ri = zeros(floor(8192/iNum),num_files);
VF_Li = zeros(floor(8192/iNum),num_files);

NF_0i = zeros(floor(8192/iNumx),num_files);
NF_Ri = zeros(floor(8192/iNumx),num_files);
NF_Li = zeros(floor(8192/iNumx),num_files);
EF_0i = zeros(floor(8192/iNumx),num_files);
EF_Ri = zeros(floor(8192/iNumx),num_files);
EF_Li = zeros(floor(8192/iNumx),num_files);

for i = 1:floor(8192/iNum)
    VF_0i(i,:) = sum(V_F_0(iNum*(i-1)+1:iNum*i,:),1);
    VF_Ri(i,:) = sum(V_F_R(iNum*(i-1)+1:iNum*i,:),1);
    VF_Li(i,:) = sum(V_F_L(iNum*(i-1)+1:iNum*i,:),1);
end
    
for i = 1:floor(8192/iNumx)
    NF_0i(i,:) = sum(N_F_0(iNumx*(i-1)+1:iNumx*i,:),1);
    NF_Ri(i,:) = sum(N_F_R(iNumx*(i-1)+1:iNumx*i,:),1);
    NF_Li(i,:) = sum(N_F_L(iNumx*(i-1)+1:iNumx*i,:),1);
    EF_0i(i,:) = sum(E_F_0(iNumx*(i-1)+1:iNumx*i,:),1);
    EF_Ri(i,:) = sum(E_F_R(iNumx*(i-1)+1:iNumx*i,:),1);
    EF_Li(i,:) = sum(E_F_L(iNumx*(i-1)+1:iNumx*i,:),1);
end

% 合并列
VF_0 = zeros(floor(8192/iNum),floor(num_files/jNum));
VF_R = zeros(floor(8192/iNum),floor(num_files/jNum));
VF_L = zeros(floor(8192/iNum),floor(num_files/jNum));

NF_0 = zeros(floor(8192/iNumx),floor(num_files/jNum));
NF_R = zeros(floor(8192/iNumx),floor(num_files/jNum));
NF_L = zeros(floor(8192/iNumx),floor(num_files/jNum));
EF_0 = zeros(floor(8192/iNumx),floor(num_files/jNum));
EF_R = zeros(floor(8192/iNumx),floor(num_files/jNum));
EF_L = zeros(floor(8192/iNumx),floor(num_files/jNum));

for j = 1:floor(num_files/jNum)
    VF_0(:,j) = sum(VF_0i(:,jNum*(j-1)+1:jNum*j),2);
    VF_R(:,j) = sum(VF_Ri(:,jNum*(j-1)+1:jNum*j),2);
    VF_L(:,j) = sum(VF_Li(:,jNum*(j-1)+1:jNum*j),2);
    
    NF_0(:,j) = sum(NF_0i(:,jNum*(j-1)+1:jNum*j),2);
    NF_R(:,j) = sum(NF_Ri(:,jNum*(j-1)+1:jNum*j),2);
    NF_L(:,j) = sum(NF_Li(:,jNum*(j-1)+1:jNum*j),2);
    EF_0(:,j) = sum(EF_0i(:,jNum*(j-1)+1:jNum*j),2);
    EF_R(:,j) = sum(EF_Ri(:,jNum*(j-1)+1:jNum*j),2);
    EF_L(:,j) = sum(EF_Li(:,jNum*(j-1)+1:jNum*j),2);
end

% 噪声120-150，计算光子数，高度
Altitude = height_num_origin(1:floor(8192/iNum))*iNum;
KM30 = size(Altitude(Altitude<30),1)+1;
KM35 = size(Altitude(Altitude<35),1)+1;
KM40 = size(Altitude(Altitude<40),1)+1;
KM75 = size(Altitude(Altitude<75),1)+1;
KM80 = size(Altitude(Altitude<80),1)+1;
KM85 = size(Altitude(Altitude<85),1)+1;
KM90 = size(Altitude(Altitude<90),1)+1;
KM95 = size(Altitude(Altitude<95),1)+1;
KM98 = size(Altitude(Altitude<98),1)+1;
KM100 = size(Altitude(Altitude<100),1)+1;
KM105 = size(Altitude(Altitude<105),1)+1;
KM115 = size(Altitude(Altitude<115),1)+1;
KM120 = size(Altitude(Altitude<120),1)+1;
KM150 = size(Altitude(Altitude<150),1)+1;
Alt30 = (height_num_origin(1:floor(8192/iNumx))*iNumx)*(sqrt(3)/2);
KMx30 = size(Alt30(Alt30<30),1)+1;
KMx35 = size(Alt30(Alt30<35),1)+1;
KMx40 = size(Alt30(Alt30<40),1)+1;
KMx75 = size(Alt30(Alt30<75),1)+1;
KMx80 = size(Alt30(Alt30<80),1)+1;
KMx85 = size(Alt30(Alt30<85),1)+1;
KMx90 = size(Alt30(Alt30<90),1)+1;
KMx95 = size(Alt30(Alt30<95),1)+1;
KMx98 = size(Alt30(Alt30<98),1)+1;
KMx100 = size(Alt30(Alt30<100),1)+1;
KMx105 = size(Alt30(Alt30<105),1)+1;
KMx115 = size(Alt30(Alt30<115),1)+1;
KMx120 = size(Alt30(Alt30<120),1)+1;
KMx150 = size(Alt30(Alt30<150),1)+1;

% 去噪声
VF0 = VF_0 - mean(VF_0(KM120:KM150,:),1);
VFR = VF_R - mean(VF_R(KM120:KM150,:),1);
VFL = VF_L - mean(VF_L(KM120:KM150,:),1);
NF0 = NF_0 - mean(NF_0(KMx120:KMx150,:),1);
NFR = NF_R - mean(NF_R(KMx120:KMx150,:),1);
NFL = NF_L - mean(NF_L(KMx120:KMx150,:),1);
EF0 = EF_0 - mean(EF_0(KMx120:KMx150,:),1);
EFR = EF_R - mean(EF_R(KMx120:KMx150,:),1);
EFL = EF_L - mean(EF_L(KMx120:KMx150,:),1);

RunRead = '第二次合并去噪OK了'

% 归一化
FL = 1;
FR = size(VF0,2);
V0 = VF0(:,FL:FR)./VF0(KM35,FL:FR);
VR = VFR(:,FL:FR)./VFR(KM35,FL:FR);
VL = VFL(:,FL:FR)./VFL(KM35,FL:FR);
N0 = NF0(:,FL:FR)./NF0(KMx35,FL:FR);
NR = NFR(:,FL:FR)./NFR(KMx35,FL:FR);
NL = NFL(:,FL:FR)./NFL(KMx35,FL:FR);
E0 = EF0(:,FL:FR)./EF0(KMx35,FL:FR);
ER = EFR(:,FL:FR)./EFR(KMx35,FL:FR);
EL = EFL(:,FL:FR)./EFL(KMx35,FL:FR);

RaySV0 = sum(VF_0(KM35:KM40,:),1);
RaySVR = sum(VF_R(KM35:KM40,:),1);
RaySVL = sum(VF_L(KM35:KM40,:),1);
RaySN0 = sum(NF_0(KMx35:KMx40,:),1);
RaySNR = sum(NF_R(KMx35:KMx40,:),1);
RaySNL = sum(NF_L(KMx35:KMx40,:),1);
RaySE0 = sum(EF_0(KMx35:KMx40,:),1);
RaySER = sum(EF_R(KMx35:KMx40,:),1);
RaySEL = sum(EF_L(KMx35:KMx40,:),1);

Tlist1 = datestr(TimeX','yyyymmdd HH:MM')
Tlist = Tlist1(FL:FR,:);
RayV0 = RaySV0(:,FL:FR);
RayVR = RaySVR(:,FL:FR);
RayVL = RaySVL(:,FL:FR);
RayN0 = RaySN0(:,FL:FR);
RayNR = RaySNR(:,FL:FR);
RayNL = RaySNL(:,FL:FR);
RayE0 = RaySE0(:,FL:FR);
RayER = RaySER(:,FL:FR);
RayEL = RaySEL(:,FL:FR);

ErayVR = RayV0./RayVR;
ErayVL = RayV0./RayVL;
ErayNR = RayN0./RayNR;
ErayNL = RayN0./RayNL;
ErayER = RayE0./RayER;
ErayEL = RayE0./RayEL;

NoiseV0 = VF_0 - VF0;
NoiseVR = VF_R - VFR;
NoiseVL = VF_L - VFL;

RunRead = '第二次归一化OK了'

%% 实验性白天投射系数
% ErayVR = ones(1,size(RayN0,2));
% ErayVL = ones(1,size(RayN0,2));
% ErayNR = ones(1,size(RayN0,2));
% ErayNL = ones(1,size(RayN0,2));
% ErayER = ones(1,size(RayN0,2));
% ErayEL = ones(1,size(RayN0,2));


%% RV RT 夜间
% 计算RV
R_V_real_V = (VR - VL) ./ (V0);
R_V_real_N = (NR - NL) ./ (N0);
R_V_real_E = (ER - EL) ./ (E0);

% 计算RT
R_T_real_V = (VR + VL) ./ (2 * V0);
R_T_real_N = (NR + NL) ./ (2 * N0);
R_T_real_E = (ER + EL) ./ (2 * E0);

RunRead = '夜间RTRV OK了'

%% RV RT 白天
% % 计算RV
% R_V_real_V = (VR.*RfR - VL.*RfL) ./ (V0.*Rf0);
% R_V_real_N = (NR.*RfR - NL.*RfL) ./ (N0.*Rf0);
% R_V_real_E = (ER.*RfR - EL.*RfL) ./ (E0.*Rf0);
% 
% % 计算RT
% R_T_real_V = (VR.*RfR + VL.*RfL) ./ (2 * (V0.*Rf0));
% R_T_real_N = (NR.*RfR + NL.*RfL) ./ (2 * (N0.*Rf0));
% R_T_real_E = (ER.*RfR + EL.*RfL) ./ (2 * (E0.*Rf0));
% 
% RunRead = '白天RTRV OK了'

%% 垂直风矫正
load PxT.mat
load PxV.mat
PxV = PxV;
PxT = PxT;
VmdfV = zeros(1,601);
RankMin = zeros(1,size(V0,2));
for jt = 1:size(V0,2)
    for jj = 1:601
        V_real_V = FV(R_V_real_V(:,jt),R_T_real_V(:,jt),PxV,jj);
        VmdfV(:,jj) = abs(mean(V_real_V(KM90:KM95,:),1));
        
    end
    minValue = min(VmdfV);
    RankMin(:,jt) = find(VmdfV == minValue);
end
RunRead = '垂直风矫正OK了'

% 记录真实频移量
TrueDnu = RankMin-300
TruePxT = zeros(21,size(V0,2));
TruePxV = zeros(21,size(V0,2));
for jTrue = 1:size(V0,2)
    TruePxT(:,jTrue) = PxT(:,RankMin(:,jTrue));
    TruePxV(:,jTrue) = PxV(:,RankMin(:,jTrue));
end
RunRead = '记录频移量OK了'

%% 不修正移频量
load PxT.mat
load PxV.mat
TruePxT = zeros(21,size(V0,2));
TruePxV = zeros(21,size(V0,2));
for jTrue = 1:size(V0,2)
    TruePxT(:,jTrue) = PxT(:,301);
    TruePxV(:,jTrue) = PxV(:,301);
end
RunRead = '不修正频移量OK了'


%% 温度彩图
Ttk = 1/24;
MapTime = datenum(TimeX)+(8*Ttk);
Temp_V = V0;

LdtT = datetime(2023, 11, 3, 18, 0, 0);
RdtT = datetime(2023, 11, 3, 24, 0, 0);
Ldt = datenum(LdtT);
Rdt = datenum(RdtT);

for jvt = 1:size(V0,2)
    Temp_V(:,jvt) = FT(R_V_real_V(:,jvt),R_T_real_V(:,jvt),TruePxT,jvt);
end
[Map, Line] = contourf(MapTime,Altitude(KM85:KM100,:),Temp_V(KM85:KM100,:));
set(Line,'LineColor','none')
ylabel(colorbar,'Temperature (K)');
colormap(jet);
xlim([Ldt,Rdt]);
Ttk = 1/24;
x_ticks = Ldt:Ttk:Rdt;
set(gca, 'XTick', x_ticks);
datetick('x','HH:MM','keeplimits','keepticks');
xlabel('Time (LT)');
ylabel('Altitude (km)')


%% 单温度廓线

Tlist(135,:)

for jvt = 1:size(V0,2)
    Temp_V(:,jvt) = FT(R_V_real_V(:,jvt),R_T_real_V(:,jvt),TruePxT,jvt);
end

syms RV RT;

PreT = Altitude;
PreV = Altitude;
PreTN = Alt30;
PreVN = Alt30;
PreTE = Alt30;
PreVE = Alt30;

Jnight = 135

for jpre = 1:size(Altitude,1)
    RV_real = R_V_real_V(jpre,Jnight);
    RT_real = R_T_real_V(jpre,Jnight);
    Ph0 = VF_0(jpre,Jnight);
    PhR = VF_R(jpre,Jnight).*ErayVR(:,Jnight);
    PhL = VF_L(jpre,Jnight).*ErayVL(:,Jnight);

    % 对 fTL 的 RT 和 RV 分别求偏导
    dT_dRTL = diff(FTerr(RV,RT),RT);
    dT_dRVL = diff(FTerr(RV,RT),RV);
    dvalueTRTL = subs(dT_dRTL,[RV RT],[RV_real RT_real]);
    dvalueTRTL = double(dvalueTRTL);
    dvalueTRVL = subs(dT_dRVL,[RV RT],[RV_real RT_real]);
    dvalueTRVL = double(dvalueTRVL);

    % 对 fVL 的 RT 和 RV 分别求偏导
    dVL_dRT = diff(FVerr(RV,RT),RT);
    dVL_dRV = diff(FVerr(RV,RT),RV);
    dvalueVRTL = subs(dVL_dRT,[RV RT],[RV_real RT_real]);
    dvalueVRTL = double(dvalueVRTL);
    dvalueVRVL = subs(dVL_dRV,[RV RT],[RV_real RT_real]);
    dvalueVRVL = double(dvalueVRVL);

    % 定义并计算Pre_TL
    Pre_TL = sqrt((PhR.*(0.5.*dvalueTRTL+dvalueTRVL).^2+PhL.*(0.5.*dvalueTRTL-dvalueTRVL).^2+Ph0.*(RT_real.*dvalueTRTL+RV_real.*dvalueTRVL).^2)./(Ph0.^2));
    Pre_TL = real(double(Pre_TL));

    % 定义并计算Pre_VL
    Pre_VL = sqrt((PhR.*(0.5.*dvalueVRTL+dvalueVRVL).^2+PhL.*(0.5.*dvalueVRTL-dvalueVRVL).^2+Ph0.*(RT_real.*dvalueVRTL+RV_real.*dvalueVRVL).^2)./(Ph0.^2));
    Pre_VL = real(double(Pre_VL));

    PreT(jpre,:) = Pre_TL;
    PreV(jpre,:) = Pre_VL;

end

fTV = figure('name','温度');
PlotT = Temp_V(:,Jnight);
PlotPre = PreT;
% plot(PlotT,Altitude,'-k','linewidth',1.5);

e=errorbar(PlotT(KM80:KM105,:),Altitude(KM80:KM105,:),PlotPre(KM80:KM105,:),'horizontal');
e.Color = 'black';
e.Marker = 'none';
e.MarkerSize = 5;
e.CapSize = 5;
e.LineWidth=1.5;

ylim([80 105]);
xlim([180 300]);
ttstr = "Temperature (K) " + Tlist(Jnight,:);
title(ttstr);
xlabel('T (K)');
ylabel('Altitude (km)')
set(gca,'FontName','arial','FontSize',10)
grid on;
legend('MoheLidar','Nrlmsise','Saber')

%% 选择夜间数据群
FNum = size(V0,2)

PreXX = zeros(FNum*3,2);

DenUp = 2000;
FigUp = 102;
FigDown = 83;
FgTL = 120;
FgTR = 320;
FgNL = -200;
FgNR = 200;
FgEL = -200;
FgER = 200;
tkN = 40;
tkE = 40;

Jnight = 1:FNum;

TlistSTR = datestr(TimeX','HHMM');


for Jnight = Jnight
    
syms RV RT;

PreT = Altitude;
PreV = Altitude;
PreTN = Alt30;
PreVN = Alt30;
PreTE = Alt30;
PreVE = Alt30;
    
    for jpre = 1:size(Alt30,1)
        RV_real = R_V_real_N(jpre,Jnight);
        RT_real = R_T_real_N(jpre,Jnight);
        Ph0 = NF_0(jpre,Jnight);
        PhR = NF_R(jpre,Jnight).*ErayNR(:,Jnight);
        PhL = NF_L(jpre,Jnight).*ErayNL(:,Jnight);

        % 对 fTL 的 RT 和 RV 分别求偏导
        dT_dRTL = diff(FTerr(RV,RT),RT);
        dT_dRVL = diff(FTerr(RV,RT),RV);
        dvalueTRTL = subs(dT_dRTL,[RV RT],[RV_real RT_real]);
        dvalueTRTL = double(dvalueTRTL);
        dvalueTRVL = subs(dT_dRVL,[RV RT],[RV_real RT_real]);
        dvalueTRVL = double(dvalueTRVL);

        % 对 fVL 的 RT 和 RV 分别求偏导
        dVL_dRT = diff(FVerr(RV,RT),RT);
        dVL_dRV = diff(FVerr(RV,RT),RV);
        dvalueVRTL = subs(dVL_dRT,[RV RT],[RV_real RT_real]);
        dvalueVRTL = double(dvalueVRTL);
        dvalueVRVL = subs(dVL_dRV,[RV RT],[RV_real RT_real]);
        dvalueVRVL = double(dvalueVRVL);

        % 定义并计算Pre_TL
        Pre_TL = sqrt((PhR.*(0.5.*dvalueTRTL+dvalueTRVL).^2+PhL.*(0.5.*dvalueTRTL-dvalueTRVL).^2+Ph0.*(RT_real.*dvalueTRTL+RV_real.*dvalueTRVL).^2)./(Ph0.^2));
        Pre_TL = real(double(Pre_TL));

        % 定义并计算Pre_VL
        Pre_VL = sqrt((PhR.*(0.5.*dvalueVRTL+dvalueVRVL).^2+PhL.*(0.5.*dvalueVRTL-dvalueVRVL).^2+Ph0.*(RT_real.*dvalueVRTL+RV_real.*dvalueVRVL).^2)./(Ph0.^2));
        Pre_VL = real(double(Pre_VL));

        PreTN(jpre,:) = Pre_TL;
        PreVN(jpre,:) = Pre_VL;

    end
    
    PreXX((Jnight-1)*3+2,1) = PreTN(KMx90,:)
    PreXX((Jnight-1)*3+2,2) = PreVN(KMx90,:)
    
    
    % 画北向风速
    Wind_N = N0;
    for jvt = Jnight
        Wind_N(:,jvt) = FV(R_V_real_N(:,jvt),R_T_real_N(:,jvt),TruePxV,jvt);
    end
    fWN = figure('name','北向风');
%     plot(Wind_N(:,Jnight).*2,Altitude*(sqrt(3)/2),'-k','linewidth',1.5);
    
    e=errorbar(Wind_N(:,Jnight).*2,Alt30,PreVN,'horizontal');
    e.Color = 'black';
    e.Marker = 'none';
    e.MarkerSize = 5;
    e.CapSize = 5;
    e.LineWidth=1.5;
    
    ylim([FigDown FigUp]);
    xlim([FgNL FgNR]);
    xtkN = xlim;
    xticks(xtkN(1):tkN:xtkN(end))
    tnstr = "Meridional Wind Velocity (m/s) " + Tlist(Jnight,:);
    title(tnstr);
    xlabel('V (m/s)');
    ylabel('Altitude (km)')
    set(gca,'FontName','Times New Roman','FontSize',12)
    grid on;
    
    
    DataWN = [floor(Alt30(KMx80:KMx105,:)),Wind_N(KMx80:KMx105,Jnight),PreVN(KMx80:KMx105,:)];
    
    PassTXT = "D:\File\20240710T\MoheNaWind\MeridionalWind\";
    mkdir(PassTXT)
    PassWN = PassTXT + "MoheNaMeridionalWind" + ReadDate + TlistSTR(Jnight,:) + ".TXT";
    FileWN = fopen(char(PassWN),'w');
    fprintf(FileWN, '# Mohe Na Meridional Wind\n');
    fprintf(FileWN, '%-7s  %-9s  %-9s\n', 'Elev', 'Velocity', 'Err(m/s)');
    fprintf(FileWN, '%7.3f  %9.1f  %9.2f\n', DataWN.');
    fclose(FileWN);

    % figure('name','北向风map');
    % MapTime = datenum(TimeX(:,Jnight));
    % MapAlt = Alt30(KMx85:KMx100,:);
    % MapWind_N = Wind_N(KMx85:KMx100,Jnight).*2;
    % [Map, Line] = contourf(MapTime,MapAlt,MapWind_N,10);
    % set(Line,'LineColor','none');
    % datetick('x','HH:MM','keepticks');
    % colormap jet
    % colorbar
    % title('Merdional Wind Velocity (m/s)');
    % xlabel('Time (UT)');
    % ylabel('Altitude (km)')
    % set(gca,'FontName','Times New Roman','FontSize',12)

    
    for jpre = 1:size(Alt30,1)
        RV_real = R_V_real_E(jpre,Jnight);
        RT_real = R_T_real_E(jpre,Jnight);
        Ph0 = EF_0(jpre,Jnight);
        PhR = EF_R(jpre,Jnight).*ErayER(:,Jnight);
        PhL = EF_L(jpre,Jnight).*ErayEL(:,Jnight);

        % 对 fTL 的 RT 和 RV 分别求偏导
        dT_dRTL = diff(FTerr(RV,RT),RT);
        dT_dRVL = diff(FTerr(RV,RT),RV);
        dvalueTRTL = subs(dT_dRTL,[RV RT],[RV_real RT_real]);
        dvalueTRTL = double(dvalueTRTL);
        dvalueTRVL = subs(dT_dRVL,[RV RT],[RV_real RT_real]);
        dvalueTRVL = double(dvalueTRVL);

        % 对 fVL 的 RT 和 RV 分别求偏导
        dVL_dRT = diff(FVerr(RV,RT),RT);
        dVL_dRV = diff(FVerr(RV,RT),RV);
        dvalueVRTL = subs(dVL_dRT,[RV RT],[RV_real RT_real]);
        dvalueVRTL = double(dvalueVRTL);
        dvalueVRVL = subs(dVL_dRV,[RV RT],[RV_real RT_real]);
        dvalueVRVL = double(dvalueVRVL);

        % 定义并计算Pre_TL
        Pre_TL = sqrt((PhR.*(0.5.*dvalueTRTL+dvalueTRVL).^2+PhL.*(0.5.*dvalueTRTL-dvalueTRVL).^2+Ph0.*(RT_real.*dvalueTRTL+RV_real.*dvalueTRVL).^2)./(Ph0.^2));
        Pre_TL = real(double(Pre_TL));

        % 定义并计算Pre_VL
        Pre_VL = sqrt((PhR.*(0.5.*dvalueVRTL+dvalueVRVL).^2+PhL.*(0.5.*dvalueVRTL-dvalueVRVL).^2+Ph0.*(RT_real.*dvalueVRTL+RV_real.*dvalueVRVL).^2)./(Ph0.^2));
        Pre_VL = real(double(Pre_VL));

        PreTE(jpre,:) = Pre_TL;
        PreVE(jpre,:) = Pre_VL;

    end
    
    PreXX((Jnight-1)*3+3,1) = PreTE(KMx90,:)
    PreXX((Jnight-1)*3+3,2) = PreVE(KMx90,:)
    
    
    % 画东向风速
    Wind_E = E0;
    for jvt = Jnight
        Wind_E(:,jvt) = FV(R_V_real_E(:,jvt),R_T_real_E(:,jvt),TruePxV,jvt);
    end
    fWE = figure('name','东向风');
%     plot(Wind_E(:,Jnight).*2,Altitude*(sqrt(3)/2),'-k','linewidth',1.5);
    
    e=errorbar(Wind_E(:,Jnight).*2,Alt30,PreVE,'horizontal');
    e.Color = 'black';
    e.Marker = 'none';
    e.MarkerSize = 5;
    e.CapSize = 5;
    e.LineWidth=1.5;
    
    ylim([FigDown FigUp]);
    xlim([FgEL FgER]);
    xtkE = xlim;
    xticks(xtkE(1):tkE:xtkE(end))
    testr = "Zonal Wind Velocity (m/s) " + Tlist(Jnight,:);
    title(testr);
    xlabel('V (m/s)');
    ylabel('Altitude (km)')
    set(gca,'FontName','Times New Roman','FontSize',12)
    grid on;
    
    
    DataWE = [floor(Alt30(KMx80:KMx105,:)),Wind_E(KMx80:KMx105,Jnight),PreVE(KMx80:KMx105,:)];
    
    PassTXT = "D:\File\20240710T\MoheNaWind\ZonalWind\";
    mkdir(PassTXT)
    PassWE = PassTXT + "MoheNaZonalWind" + ReadDate + TlistSTR(Jnight,:) + ".TXT";
    FileWE = fopen(char(PassWE),'w');
    fprintf(FileWE, '# Mohe Na Zonal Wind\n');
    fprintf(FileWE, '%-7s  %-9s  %-9s\n', 'Elev', 'Velocity', 'Err(m/s)');
    fprintf(FileWE, '%7.3f  %9.1f  %9.2f\n', DataWE.');
    fclose(FileWE);
    
close all;

end

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);





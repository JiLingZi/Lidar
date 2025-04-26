% 读取光子数矩阵
ReadDate = '20240412';
fldstart = ['F:\RawData\ZWDATA\MOHEnew\Na\',ReadDate,'\'];
% folder = [fldstart,'Na\other\'];
folder = [fldstart,'Na\'];
files = dir(fullfile(folder, '*.dat'));
num_files = length(files);% 为了统一文件数，减一

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
rd1 = 'OK了'

%% 合并时间文件数
jNum = 60;
iNum = 32;

TimeList = cell(1,floor(num_files/jNum));
for j = 1:floor(num_files/jNum)
	TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
end

TimeX = datetime(TimeList, 'InputFormat', ':yyyy-MM-dd''T''HH:mm:ss.SSS');

% 获取原始高度矩阵
height_num_origin = Na_data(1:8192,1);
rd2 = 'OK了'

%% 合并行
VF_0i = zeros(floor(8192/iNum),num_files);
VF_Ri = zeros(floor(8192/iNum),num_files);
VF_Li = zeros(floor(8192/iNum),num_files);
NF_0i = zeros(floor(8192/iNum),num_files);
NF_Ri = zeros(floor(8192/iNum),num_files);
NF_Li = zeros(floor(8192/iNum),num_files);
EF_0i = zeros(floor(8192/iNum),num_files);
EF_Ri = zeros(floor(8192/iNum),num_files);
EF_Li = zeros(floor(8192/iNum),num_files);
for i = 1:floor(8192/iNum)
    VF_0i(i,:) = sum(V_F_0(iNum*(i-1)+1:iNum*i,:),1);
    VF_Ri(i,:) = sum(V_F_R(iNum*(i-1)+1:iNum*i,:),1);
    VF_Li(i,:) = sum(V_F_L(iNum*(i-1)+1:iNum*i,:),1);
    NF_0i(i,:) = sum(N_F_0(iNum*(i-1)+1:iNum*i,:),1);
    NF_Ri(i,:) = sum(N_F_R(iNum*(i-1)+1:iNum*i,:),1);
    NF_Li(i,:) = sum(N_F_L(iNum*(i-1)+1:iNum*i,:),1);
    EF_0i(i,:) = sum(E_F_0(iNum*(i-1)+1:iNum*i,:),1);
    EF_Ri(i,:) = sum(E_F_R(iNum*(i-1)+1:iNum*i,:),1);
    EF_Li(i,:) = sum(E_F_L(iNum*(i-1)+1:iNum*i,:),1);
end
rd3 = 'OK了'

%% 合并列
VF_0 = zeros(floor(8192/iNum),floor(num_files/jNum));
VF_R = zeros(floor(8192/iNum),floor(num_files/jNum));
VF_L = zeros(floor(8192/iNum),floor(num_files/jNum));
NF_0 = zeros(floor(8192/iNum),floor(num_files/jNum));
NF_R = zeros(floor(8192/iNum),floor(num_files/jNum));
NF_L = zeros(floor(8192/iNum),floor(num_files/jNum));
EF_0 = zeros(floor(8192/iNum),floor(num_files/jNum));
EF_R = zeros(floor(8192/iNum),floor(num_files/jNum));
EF_L = zeros(floor(8192/iNum),floor(num_files/jNum));
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
rd4 = 'OK了'

%% 噪声120-150，计算光子数，高度
Altitude = height_num_origin(1:floor(8192/iNum))*iNum;
KM30 = size(Altitude(Altitude<30),1)+1;
KM35 = size(Altitude(Altitude<35),1)+1;
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
Alt30 = Altitude*(sqrt(3)/2);
KMx30 = size(Alt30(Alt30<30),1)+1;
KMx35 = size(Alt30(Alt30<35),1)+1;
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

%% 去噪声
VF0 = VF_0 - mean(VF_0(KM120:KM150,:),1);
VFR = VF_R - mean(VF_R(KM120:KM150,:),1);
VFL = VF_L - mean(VF_L(KM120:KM150,:),1);
NF0 = NF_0 - mean(NF_0(KM120:KM150,:),1);
NFR = NF_R - mean(NF_R(KM120:KM150,:),1);
NFL = NF_L - mean(NF_L(KM120:KM150,:),1);
EF0 = EF_0 - mean(EF_0(KM120:KM150,:),1);
EFR = EF_R - mean(EF_R(KM120:KM150,:),1);
EFL = EF_L - mean(EF_L(KM120:KM150,:),1);

%% 归一化
V0 = VF0./VF0(KM35,:);
VR = VFR./VFR(KM35,:);
VL = VFL./VFL(KM35,:);
N0 = NF0./NF0(KM35,:);
NR = NFR./NFR(KM35,:);
NL = NFL./NFL(KM35,:);
E0 = EF0./EF0(KM35,:);
ER = EFR./EFR(KM35,:);
EL = EFL./EFL(KM35,:);

%% 
figure('name','E三频检查','position',[600 50 350 800])
pf = 1;
plot(V0(:,pf), Altitude, '-k', 'linewidth',1.5)
hold on
plot(VR(:,pf), Altitude, '-b', 'linewidth',1.5)
hold on
plot(VL(:,pf), Altitude, '-r', 'linewidth',1.5)
grid on
set(gca, 'XScale', 'log');
ylim([30 125])
% ttstr = string(pf) + ' | ' + TimeList(pf);
Tits = 'Vertical ';
ttstr = "Mohe Na " + TimeList(pf);
title('(c) East')
set(gca,'FontSize',12.5);
% ylim([3000 7000])
xlabel('Photon Counts')
ylabel('Altitude (km)')
legend('\nu_{0}','\nu_{+}','\nu_{-}')

%% RV RT
% 计算RV
R_V_real_V = (VR - VL) ./ (V0);
R_V_real_N = (NR - NL) ./ (N0);
R_V_real_E = (ER - EL) ./ (E0);

% 计算RT
R_T_real_V = (VR + VL) ./ (2 * V0);
R_T_real_N = (NR + NL) ./ (2 * N0);
R_T_real_E = (ER + EL) ./ (2 * E0);

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
rd5 = 'OK了'

%% 记录真实频移量
TrueDnu = RankMin-300
TruePxT = zeros(21,size(V0,2));
TruePxV = zeros(21,size(V0,2));
for jTrue = 1:size(V0,2)
    TruePxT(:,jTrue) = PxT(:,RankMin(:,jTrue));
    TruePxV(:,jTrue) = PxV(:,RankMin(:,jTrue));
end
rd6 = 'OK了'

%%
% CnuMin = 0.69;
% R_V_real_V = (VR.*CnuMin - VL.*CnuMin) ./ (V0);
% R_V_real_N = (NR.*CnuMin - NL.*CnuMin) ./ (N0);
% R_V_real_E = (ER.*CnuMin - EL.*CnuMin) ./ (E0);
% R_T_real_V = (VR.*CnuMin + VL.*CnuMin) ./ (2 * V0);
% R_T_real_N = (NR.*CnuMin + NL.*CnuMin) ./ (2 * N0);
% R_T_real_E = (ER.*CnuMin + EL.*CnuMin) ./ (2 * E0);

%% 选择夜间数据群
Jnight = 1:size(V0,2);
Jnight = 1:3;

%% 画垂直风速
Wind_V = V0;
for jvt = Jnight
    Wind_V(:,jvt) = FV(R_V_real_V(:,jvt),R_T_real_V(:,jvt),TruePxV,jvt);
end
figure('name','垂直风');
plot(Wind_V(:,Jnight),Altitude,'-k','linewidth',2);
ylim([85 100]);
xlim([-20 20]);
title('Vertical Wind Velocity (m/s)');
xlabel('V (m/s)');
ylabel('Altitude (km)')
set(gca,'FontName','黑体','FontSize',12)
grid on;

figure('name','垂直风map');
MapTime = datenum(TimeX(:,Jnight));
MapAlt = Altitude(KM85:KM100,:);
MapWind_V = Wind_V(KM85:KM100,Jnight);
[Map, Line] = contourf(MapTime,MapAlt,MapWind_V,10);
set(Line,'LineColor','none');
datetick('x','HH:MM','keepticks');
colormap jet
colorbar
title('Vertical Wind Velocity (m/s)');
xlabel('Time (UT)');
ylabel('Altitude (km)')
set(gca,'FontName','黑体','FontSize',12)

%% 画北向风速
Wind_N = V0;
for jvt = Jnight
    Wind_N(:,jvt) = FV(R_V_real_N(:,jvt),R_T_real_N(:,jvt),TruePxV,jvt);
end
figure('name','北向风');
plot(Wind_N(:,Jnight).*2,Altitude*(sqrt(3)/2),'-k','linewidth',2);
ylim([85 100]);
% xlim([-20 20]);
title('Meridional Wind Velocity (m/s)');
xlabel('V (m/s)');
ylabel('Altitude (km)')
set(gca,'FontName','黑体','FontSize',12)
grid on;

figure('name','北向风map');
MapTime = datenum(TimeX(:,Jnight));
MapAlt = Alt30(KMx85:KMx100,:);
MapWind_N = Wind_N(KMx85:KMx100,Jnight).*2;
[Map, Line] = contourf(MapTime,MapAlt,MapWind_N,10);
set(Line,'LineColor','none');
datetick('x','HH:MM','keepticks');
colormap jet
colorbar
title('Merdional Wind Velocity (m/s)');
xlabel('Time (UT)');
ylabel('Altitude (km)')
set(gca,'FontName','黑体','FontSize',12)

%% 画东向风速
Wind_E = V0;
for jvt = Jnight
    Wind_E(:,jvt) = FV(R_V_real_E(:,jvt),R_T_real_E(:,jvt),TruePxV,jvt);
end
figure('name','东向风');
plot(Wind_E(:,Jnight).*2,Altitude*(sqrt(3)/2),'-k','linewidth',2);
ylim([85 100]);
% xlim([-20 20]);
title('Zonal Wind Velocity (m/s)');
xlabel('V (m/s)');
ylabel('Altitude (km)')
set(gca,'FontName','黑体','FontSize',12)
grid on;

figure('name','东向风map');
MapTime = datenum(TimeX(:,Jnight));
MapAlt = Alt30(KMx85:KMx100,:);
MapWind_E = Wind_E(KMx85:KMx100,Jnight).*2;
[Map, Line] = contourf(MapTime,MapAlt,MapWind_E,10);
set(Line,'LineColor','none');
datetick('x','HH:MM','keepticks');
colormap jet
colorbar
title('Zonal Wind Velocity (m/s)');
xlabel('Time (UT)');
ylabel('Altitude (km)')
set(gca,'FontName','黑体','FontSize',12)

%% 画温度
Temp_V = V0;
for jvt = Jnight
    Temp_V(:,jvt) = FT(R_V_real_V(:,jvt),R_T_real_V(:,jvt),TruePxT,jvt);
end
figure('name','温度');
plot(Temp_V(:,Jnight),Altitude,'-k','linewidth',2);
ylim([84 105]);
% xlim([-20 20]);
title('Temperature (K)');
xlabel('T (K)');
ylabel('Altitude (km)')
set(gca,'FontName','黑体','FontSize',12)
grid on;

figure('name','温度map');
MapTime = datenum(TimeX(:,Jnight));
MapAlt = Altitude(KM85:KM100,:);
MapTemp_V = Temp_V(KM85:KM100,Jnight);
[Map, Line] = contourf(MapTime,MapAlt,MapTemp_V,10);
set(Line,'LineColor','none');
datetick('x','HH:MM','keepticks');
colormap jet
colorbar
title('Temperature (K)');
xlabel('Time (UT)');
ylabel('Altitude (km)')
set(gca,'FontName','黑体','FontSize',12)







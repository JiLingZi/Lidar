PF = 12;
PhOrigin = VF0(:,PF);
PhKeep = VF0(:,PF);
T_ZCk = T_ZCheck(1:end,:)+10;

ins = ones(size(PhOrigin,1),1);
OnDleT = zeros(size(PhOrigin,1),1);
sns = ones(1,201);
RankMin = ins;
jnsk = 0.5:0.01:2.5;

for jc = 45:-1:10
    for jns = 0.5:0.01:2.5
%         if abs(T_ZCk(jc,:)-T_ZNa(jc,:))>3
            Ph = [PhKeep(1:jc-1,:);PhKeep(jc,:).*jns;PhOrigin(jc+1:end,:)];
            Z = AltVX;                           % 高度矩阵
            Z_0 = Z_0_NRL*1e-3;                     % 密度参考高度
            N_Z = Ph;                               % 30-60 km光子数
            N_B = mean(Ph(i_120:i_140,:));          % 光子数平均噪声水平
            N_Z_0 = Ph(i_40,:);                     % 密度参考高度处光子数
            Rho_ZNa = (Rho_Z_0 .* ((Z.^2)./(Z_0^2)) .* ((N_Z-N_B)./(N_Z_0-N_B)));

            Curve = Rho_ZNa.*G;

            Rho_Z_top = Rho_ZNa(i_70,:);
            T_ZNa = (1:(i_70-1))';
            for i = 1:(i_70-1)
                T_ZNa(i,:) = (Rho_Z_top/Rho_ZNa(i,:))*T_Z_top + ((M/R)/Rho_ZNa(i,:))*((sum(Curve(i:i_70,:)))*Alt_Res);
            end
            PhOrigin= Ph;
%         end
            DleteT = abs(T_ZCk(jc,:)-T_ZNa(jc,:));
            
            x1 = 0.5;
            x2 = 2.5;
            y1 = 1;
            y2 = 201;
            a = (y2 - y1) / (x2 - x1);
            b = y1 - a * x1;
            y = round(a * jns + b);
            sns(:,y) = DleteT;
            Checksns = sns';
            
            minValue = min(sns);  
    end
    OnDleT(jc,:) = minValue;
    RankMin(jc,:) = find(sns == minValue);
    ts = find(sns == minValue);
    ins(jc,:) = jnsk(:,RankMin(jc,:));
    PhOrigin(jc,:)= PhKeep(jc,:).*ins(jc,:);
plot(T_ZNa,AltTZ,'-b','linewidth',2)
hold on;
grid on;
plot(T_ZCk,AltTZ,'-r','linewidth',2)
ylim([25 55])
xlabel('Temperature (K)');
ylabel('Altitude (km)');

end

% ins = PhOrigin./VF0(:,PF);

%%
PhFinish = VF0(:,PF).*ins;
Ph = PhFinish;

Z = AltVX;                           % 高度矩阵
Z_0 = Z_0_NRL*1e-3;                     % 密度参考高度
N_Z = Ph;                               % 30-60 km光子数
N_B = mean(Ph(i_120:i_140,:));          % 光子数平均噪声水平
N_Z_0 = Ph(i_40,:);                     % 密度参考高度处光子数
Rho_ZNa = (Rho_Z_0 .* ((Z.^2)./(Z_0^2)) .* ((N_Z-N_B)./(N_Z_0-N_B)));

Curve = Rho_ZNa.*G;

Rho_Z_top = Rho_ZNa(i_70,:);
T_ZNa = (1:(i_70-1))';
for i = 1:(i_70-1)
    T_ZNa(i,:) = (Rho_Z_top/Rho_ZNa(i,:))*T_Z_top + ((M/R)/Rho_ZNa(i,:))*((sum(Curve(i:i_70,:)))*Alt_Res);
end

plot(T_ZNa,AltTZ,'-b','linewidth',2)
hold on;
grid on;
plot(T_ZCk,AltTZ,'-r','linewidth',2)
hold on;

Ph = PhKeep;
Z = AltVX;                           % 高度矩阵
Z_0 = Z_0_NRL*1e-3;                     % 密度参考高度
N_Z = Ph;                               % 30-60 km光子数
N_B = mean(Ph(i_120:i_140,:));          % 光子数平均噪声水平
N_Z_0 = Ph(i_40,:);                     % 密度参考高度处光子数
Rho_ZNa = (Rho_Z_0 .* ((Z.^2)./(Z_0^2)) .* ((N_Z-N_B)./(N_Z_0-N_B)));
Curve = Rho_ZNa.*G;
Rho_Z_top = Rho_ZNa(i_70,:);
T_ZNa = (1:(i_70-1))';
for i = 1:(i_70-1)
    T_ZNa(i,:) = (Rho_Z_top/Rho_ZNa(i,:))*T_Z_top + ((M/R)/Rho_ZNa(i,:))*((sum(Curve(i:i_70,:)))*Alt_Res);
end
plot(T_ZNa,AltTZ,'--b','linewidth',2)

ylim([25 65])
xlabel('Temperature (K)');
ylabel('Altitude (km)');
legend('Na After','K','Na Before')
ttstr = "Rayleight Temperature Modify " + string(TlistRay(pf));
title(ttstr);

%% 
p3f = 300;
AltBin = 8192;

Hins = repelem(ins,32);

AltV = height_num_origin;
AltX1 = height_num_origin*(sqrt(3)/2);
AltX = AltV(1:AltBin,:);

AltXVoj = AltX;

figure('name','Hins')
plot(AltXVoj,Hins,'-b','linewidth',1.5)
hold on;
HinsSm = smoothdata(Hins, 'movmean', 64);
plot(AltXVoj,HinsSm,'-r','linewidth',1.5)
xlim([25 55])
xlabel('Altitude (km)');
ylabel('Correction factor');
ttstr = "Correction Factor " + string(TlistRay(pf));
title(ttstr);

figure('name','F0RL','position',[300 300 1200 350])
subplot(1,3,1)
phoj0 = V_F_0(:,p3f);
phsm1 = phoj0(1:AltBin,:);
plot(AltXVoj,phsm1,'-b','linewidth',1.5)
hold on;
phsm2 = smoothdata(phsm1, 'movmean', 128);
plot(AltXVoj,phsm2,'-r','linewidth',1.5)
xlim([25 55])
xlabel('Altitude (km)');
ylabel('Counts');
ttstr = "f0 " + string(TlistRay(pf));
title(ttstr);
subplot(1,3,2)
phojR = V_F_R(:,p3f);
phsm1 = phojR(1:AltBin,:);
plot(AltXVoj,phsm1,'-b','linewidth',1.5)
hold on;
phsm2 = smoothdata(phsm1, 'movmean', 128);
plot(AltXVoj,phsm2,'-r','linewidth',1.5)
xlim([25 55])
xlabel('Altitude (km)');
ylabel('Counts');
ttstr = "f+ " + string(TlistRay(pf));
title(ttstr);
subplot(1,3,3)
phojR = V_F_L(:,p3f);
phsm1 = phojR(1:AltBin,:);
plot(AltXVoj,phsm1,'-b','linewidth',1.5)
hold on;
phsm2 = smoothdata(phsm1, 'movmean', 128);
plot(AltXVoj,phsm2,'-r','linewidth',1.5)
xlim([25 55])
xlabel('Altitude (km)');
ylabel('Counts');
ttstr = "f- " + string(TlistRay(pf));
title(ttstr);

PhinsTk = [AltXVoj,phsm2,HinsSm];
PhsTkUsf = PhinsTk(652:1280,:);
PhHins2 = flip(PhsTkUsf(:,2));
PhHins3 = flip(PhsTkUsf(:,3));

Modify_V_F_0 = V_F_0;
Modify_V_F_R = V_F_R;
Modify_V_F_L = V_F_L;
Modify_N_F_0 = N_F_0;
Modify_N_F_R = N_F_R;
Modify_N_F_L = N_F_L;
Modify_E_F_0 = E_F_0;
Modify_E_F_R = E_F_R;
Modify_E_F_L = E_F_L;

indexAbove = ones(8192,1);

for jj = 2820:8192
    if V_F_0(jj,p3f)>PhHins2(1,:)
    CheckVal = V_F_0(jj,p3f);
    indexAbove(jj,:) = find(PhHins2 >= CheckVal, 1);
    Modify_V_F_0(jj,p3f) = V_F_0(jj,p3f).*PhHins3(indexAbove(jj,:),:);
    end
    if V_F_R(jj,p3f)>PhHins2(1,:)
    CheckVal = V_F_R(jj,p3f);
    indexAbove(jj,:) = find(PhHins2 >= CheckVal, 1);
    Modify_V_F_R(jj,p3f) = V_F_R(jj,p3f).*PhHins3(indexAbove(jj,:),:);
    end
    if V_F_L(jj,p3f)>PhHins2(1,:)
    CheckVal = V_F_L(jj,p3f);
    indexAbove(jj,:) = find(PhHins2 >= CheckVal, 1);
    Modify_V_F_L(jj,p3f) = V_F_L(jj,p3f).*PhHins3(indexAbove(jj,:),:);
    end
end

% Af = [1,2,3,4,5,6,7];
% Bf = [2.5 3.9 4.6]';
% indexAbove = find(Af >= Bf(1), 1)

figure('name','修正后比较','position',[300 300 1200 350])
subplot(1,3,1)
plot(AltX1,Modify_V_F_0(:,p3f),'-r','linewidth',1.5);
hold on;
plot(AltX1,V_F_0(:,p3f),'-b','linewidth',1.5);
xlim([25 115])
grid on;
legend('After','Before')
subplot(1,3,2)
plot(AltX1,Modify_V_F_R(:,p3f),'-r','linewidth',1.5);
hold on;
plot(AltX1,V_F_R(:,p3f),'-b','linewidth',1.5);
xlim([25 115])
grid on;
legend('After','Before')
subplot(1,3,3)
plot(AltX1,Modify_V_F_L(:,p3f),'-r','linewidth',1.5);
hold on;
plot(AltX1,V_F_L(:,p3f),'-b','linewidth',1.5);
xlim([25 115])
grid on;
legend('After','Before')

%% 修改初始光子数，谨慎运行

V_F_0 = Modify_V_F_0;
V_F_R = Modify_V_F_R;
V_F_L = Modify_V_F_L;
N_F_0 = Modify_N_F_0;
N_F_R = Modify_N_F_R;
N_F_L = Modify_N_F_L;
E_F_0 = Modify_E_F_0;
E_F_R = Modify_E_F_R;
E_F_L = Modify_E_F_L;



%% 标定系数
Rf0 = 1.1026;
RfL = 1.2123;
RfR = 0.9621;

RunRead = '系数标定OK了'

% 合并时间文件数，15、32、37
jNum = 1;
iNum = 1;
iNumx = 1;

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

Tlist1 = datestr(TimeX','yyyymmdd HH:MM');
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



%% 单温度廓线
Tfind = TimeX';

% Tlist(300,:)

for jvt = 1:size(V0,2)
    Temp_Plot(:,jvt) = FT(R_V_real_E(:,jvt),R_T_real_E(:,jvt),TruePxT,jvt);
end

syms RV RT;

PreT = Altitude;
PreV = Altitude;
PreTN = Alt30;
PreVN = Alt30;
PreTE = Alt30;
PreVE = Alt30;

Jnight = 300

PhF0 = EF0;
PhFR = EFR;
PhFL = EFL;

AltF = Alt30;

figure('name','origin','position',[300 300 1200 350])
subplot(1,3,1)
plot(PhF0(:,Jnight),AltF)
ylim([75 115]);
grid on;
xlabel('Photon Counts');
ylabel('Altitude (km)');
title('(a) f0');
subplot(1,3,2)
plot(PhFR(:,Jnight),AltF)
ylim([75 115]);
grid on;
xlabel('Photon Counts');
ylabel('Altitude (km)');
title('(a) f+');
subplot(1,3,3)
plot(PhFL(:,Jnight),AltF)
ylim([75 115]);
grid on;
xlabel('Photon Counts');
ylabel('Altitude (km)');
title('(a) f-');

% for jpre = 1:size(Altitude,1)
%     RV_real = R_V_real_V(jpre,Jnight);
%     RT_real = R_T_real_V(jpre,Jnight);
%     Ph0 = VF_0(jpre,Jnight);
%     PhR = VF_R(jpre,Jnight).*ErayVR(:,Jnight);
%     PhL = VF_L(jpre,Jnight).*ErayVL(:,Jnight);
% 
%     % 对 fTL 的 RT 和 RV 分别求偏导
%     dT_dRTL = diff(FTerr(RV,RT),RT);
%     dT_dRVL = diff(FTerr(RV,RT),RV);
%     dvalueTRTL = subs(dT_dRTL,[RV RT],[RV_real RT_real]);
%     dvalueTRTL = double(dvalueTRTL);
%     dvalueTRVL = subs(dT_dRVL,[RV RT],[RV_real RT_real]);
%     dvalueTRVL = double(dvalueTRVL);
% 
%     % 对 fVL 的 RT 和 RV 分别求偏导
%     dVL_dRT = diff(FVerr(RV,RT),RT);
%     dVL_dRV = diff(FVerr(RV,RT),RV);
%     dvalueVRTL = subs(dVL_dRT,[RV RT],[RV_real RT_real]);
%     dvalueVRTL = double(dvalueVRTL);
%     dvalueVRVL = subs(dVL_dRV,[RV RT],[RV_real RT_real]);
%     dvalueVRVL = double(dvalueVRVL);
% 
%     % 定义并计算Pre_TL
%     Pre_TL = sqrt((PhR.*(0.5.*dvalueTRTL+dvalueTRVL).^2+PhL.*(0.5.*dvalueTRTL-dvalueTRVL).^2+Ph0.*(RT_real.*dvalueTRTL+RV_real.*dvalueTRVL).^2)./(Ph0.^2));
%     Pre_TL = real(double(Pre_TL));
% 
%     % 定义并计算Pre_VL
%     Pre_VL = sqrt((PhR.*(0.5.*dvalueVRTL+dvalueVRVL).^2+PhL.*(0.5.*dvalueVRTL-dvalueVRVL).^2+Ph0.*(RT_real.*dvalueVRTL+RV_real.*dvalueVRVL).^2)./(Ph0.^2));
%     Pre_VL = real(double(Pre_VL));
% 
%     PreT(jpre,:) = Pre_TL;
%     PreV(jpre,:) = Pre_VL;
% 
% end

fTV = figure('name','温度');
PlotT = Temp_Plot(:,Jnight);
PlotPre = PreT;

PlotT(PlotT<0)=200;
PlotT(PlotT>5000)=200;
plot(PlotT,AltF,'-k','linewidth',1.5);

% e=errorbar(PlotT(KM80:KM105,:),Altitude(KM80:KM105,:),PlotPre(KM80:KM105,:),'horizontal');
% e.Color = 'black';
% e.Marker = 'none';
% e.MarkerSize = 5;
% e.CapSize = 5;
% e.LineWidth=1.5;

ylim([80 105]);
% xlim([180 300]);
ttstr = "Temperature (K) " + Tlist(Jnight,:);
title(ttstr);
xlabel('T (K)');
ylabel('Altitude (km)')
set(gca,'FontName','arial','FontSize',10)
grid on;
legend('MoheLidar','Nrlmsise','Saber')


%% 单点温度图
% 计算RV
R_V_real_V = (VR - VL) ./ (V0);
R_V_real_N = (NR - NL) ./ (N0);
R_V_real_E = (ER(3458,:) - EL(3456,:)) ./ (E0(3457,:));

% 计算RT
R_T_real_V = (VR + VL) ./ (2 * V0);
R_T_real_N = (NR + NL) ./ (2 * N0);
R_T_real_E = (ER(3458,:) + EL(3456,:)) ./ (2 * E0(3457,:));


for jvt = 1:size(V0,2)
    Temp_Single(:,jvt) = FT(R_V_real_E(:,jvt),R_T_real_E(:,jvt),TruePxT,jvt);
end

T_Sig = Temp_Single(:,300)

%% 饱和验证

figure('name','origin','position',[300 50 900 800])
subplot(1,3,1)
plot(PhF0(:,Jnight),AltF)
ylim([25 115]);
grid on;
xlabel('Photon Counts');
ylabel('Altitude (km)');
title('(a) f0');
subplot(1,3,2)
plot(PhFR(:,Jnight),AltF)
ylim([25 115]);
grid on;
xlabel('Photon Counts');
ylabel('Altitude (km)');
title('(a) f+');
subplot(1,3,3)
plot(PhFL(:,Jnight),AltF)
ylim([25 115]);
grid on;
xlabel('Photon Counts');
ylabel('Altitude (km)');
title('(a) f-');
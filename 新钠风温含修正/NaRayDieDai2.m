PF = 1;
PhMove = VF0(:,PF);
PhMove2 = smoothdata(PhMove, 'movmean', 5);
PhOrigin = PhMove2;
PhKeep = PhMove2;
T_ZCk = T_ZCheck(1:end,:)+2;

ins = ones(size(PhOrigin,1),1);
OnDleT = zeros(size(PhOrigin,1),1);
sns = ones(1,201);
RankMin = ins;
jnsk = 0.5:0.01:2.5;

for jc = 35:-1:10
    for jns = 0.5:0.01:2.5
%         if abs(T_ZCk(jc,:)-T_ZNa(jc,:))>3
            Ph = [PhKeep(1:jc-1,:);PhKeep(jc,:).*jns;PhOrigin(jc+1:end,:)];
%             Ph = smoothdata(Ph, 'movmean', 5);
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

%% 直接光子数拟合
PhOtXog = PhOt(19:35,:);
PhNaY = PhNa(19:35,:);
FktNaOt = PhOtXog(end,:)./PhNaY(end,:);
PhOtX = PhOtXog./FktNaOt;
plot(PhOtX,PhNaY,'-m','linewidth',1.5);
xlabel('Reference photon number');
ylabel('Na photon number')
title('18.7 km to 34.4 km')
xlim([0 3e6]);ylim([0 3e6]);
grid on; box on; hold on;

% cftool(PhOtX,PhNaY)

FexpPh = fit(PhOtX,PhNaY, 'exp2');
plot(FexpPh,'-g')


%%
PhFinish = PhMove2.*ins;
Ph = PhFinish;

Z = AltVX;                              % 高度矩阵
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
ttstr = "Rayleigh Temperature Modify " + string(TlistRay(pf));
title(ttstr);

%% 
p3f = 1;
AltBin = 256;

Hins = ins;

AltV = AltVX;
% AltX1 = AltVX*(sqrt(3)/2);
AltX = AltV(1:AltBin,:);

AltXVoj = AltX;

figure('name','Hins')
plot(AltXVoj,Hins,'-b','linewidth',1.5)
hold on;
HinsSm = smoothdata(Hins, 'movmean', 5);
plot(AltXVoj,HinsSm,'-r','linewidth',1.5)
xlim([25 55])
xlabel('Altitude (km)');
ylabel('Correction factor');
ttstr = "Correction Factor " + string(TlistRay(pf));
title(ttstr);

figure('name','F0RL','position',[300 300 1200 350])
subplot(1,3,1)
phoj0 = VF0(:,p3f);
phsm1 = phoj0(1:AltBin,:);
plot(AltXVoj,phsm1,'-b','linewidth',1.5)
hold on;
phsm2 = smoothdata(phsm1, 'movmean', 5);
plot(AltXVoj,phsm2,'-r','linewidth',1.5)
xlim([25 55])
xlabel('Altitude (km)');
ylabel('Counts');
ttstr = "f0 " + string(TlistRay(pf));
title(ttstr);
subplot(1,3,2)
phojR = VFR(:,p3f);
phsm1 = phojR(1:AltBin,:);
plot(AltXVoj,phsm1,'-b','linewidth',1.5)
hold on;
phsm2 = smoothdata(phsm1, 'movmean', 5);
plot(AltXVoj,phsm2,'-r','linewidth',1.5)
xlim([25 55])
xlabel('Altitude (km)');
ylabel('Counts');
ttstr = "f+ " + string(TlistRay(pf));
title(ttstr);
subplot(1,3,3)
phojR = VFL(:,p3f);
phsm1 = phojR(1:AltBin,:);
plot(AltXVoj,phsm1,'-b','linewidth',1.5)
hold on;
phsm2 = smoothdata(phsm1, 'movmean', 5);
plot(AltXVoj,phsm2,'-r','linewidth',1.5)
xlim([25 55])
xlabel('Altitude (km)');
ylabel('Counts');
ttstr = "f- " + string(TlistRay(pf));
title(ttstr);

PhinsTk = [AltXVoj,phsm2,HinsSm];
PhsTkUsf = PhinsTk(19:35,:);
PhHins2 = flip(PhsTkUsf(:,2));
PhHins3 = flip(PhsTkUsf(:,3));

figure('name','光子数与修正系数关系')
plot(PhHins2,PhHins3,'-k','linewidth',2);
grid on;hold on;
Fexp2 = fit(PhHins2,PhHins3, 'exp2');
plot(Fexp2,'-r')
xlabel('Photon Number');
ylabel('Correction factor')
legend('Correction factor','fitted curve')

PhHinsOrigin = PhHins2./(28*32);
Fexp2Ogn = fit(PhHinsOrigin,PhHins3, 'exp2');

figure('name','原始光子数与修正系数关系')
plot(PhHinsOrigin,PhHins3,'-k','linewidth',2);
grid on;hold on;
plot(Fexp2Ogn,'-r')
xlabel('Photon Number');
ylabel('Correction factor')
legend('Correction factor','fitted curve')

% Factor = 4.946*exp(( -9.645e-08)*CntX) + (-3.954)*exp((-1.568e-07)*CntX);

Modify_V_F_0 = VF0;
Modify_V_F_R = VFR;
Modify_V_F_L = VFL;
Modify_N_F_0 = NF0;
Modify_N_F_R = NFR;
Modify_N_F_L = NFL;
Modify_E_F_0 = EF0;
Modify_E_F_R = EFR;
Modify_E_F_L = EFL;

indexAbove = ones(265,1);

for jj = 20:117
    if 	VF0(jj,p3f)>PhHins2(1,:)
    CntX = VF0(jj,p3f);
    FitFactor = 4.946*exp(( -9.645e-08)*CntX) + (-3.954)*exp((-1.568e-07)*CntX);
    Modify_V_F_0(jj,p3f) = VF0(jj,p3f).*FitFactor;
    end
    if VFR(jj,p3f)>PhHins2(1,:)
    CntX = VFR(jj,p3f);
    FitFactor = 4.946*exp(( -9.645e-08)*CntX) + (-3.954)*exp((-1.568e-07)*CntX);
    Modify_V_F_R(jj,p3f) = VFR(jj,p3f).*FitFactor;
    end
    if VFL(jj,p3f)>PhHins2(1,:)
    CntX = VFL(jj,p3f);
    FitFactor = 4.946*exp(( -9.645e-08)*CntX) + (-3.954)*exp((-1.568e-07)*CntX);
    Modify_V_F_L(jj,p3f) = VFL(jj,p3f).*FitFactor;
    end
end

% for jj = 20:117
%     if 	VF0(jj,p3f)>PhHins2(1,:)
%     CheckVal = VF0(jj,p3f);
%     indexAbove(jj,:) = find(PhHins2 >= CheckVal, 1);
%     Modify_V_F_0(jj,p3f) = VF0(jj,p3f).*PhHins3(indexAbove(jj,:),:);
%     end
%     if VFR(jj,p3f)>PhHins2(1,:)
%     CheckVal = VFR(jj,p3f);
%     indexAbove(jj,:) = find(PhHins2 >= CheckVal, 1);
%     Modify_V_F_R(jj,p3f) = VFR(jj,p3f).*PhHins3(indexAbove(jj,:),:);
%     end
%     if VFL(jj,p3f)>PhHins2(1,:)
%     CheckVal = VFL(jj,p3f);
%     indexAbove(jj,:) = find(PhHins2 >= CheckVal, 1);
%     Modify_V_F_L(jj,p3f) = VFL(jj,p3f).*PhHins3(indexAbove(jj,:),:);
%     end
% end

% Af = [1,2,3,4,5,6,7];
% Bf = [2.5 3.9 4.6]';
% indexAbove = find(Af >= Bf(1), 1)

figure('name','修正后比较','position',[50 50 1200 800])
subplot(3,1,1)
plot(AltX,Modify_V_F_0(:,p3f),'-r','linewidth',1.5);
hold on;
plot(AltX,VF0(:,p3f),'-b','linewidth',1.5);
hold on;
plot([20 115],[1.4240e+05 1.4240e+05],'--g','linewidth',1);
xlim([20 115])
grid on;
% set(gca, 'YScale', 'log');
legend('After','Before')
subplot(3,1,2)
plot(AltX,Modify_V_F_R(:,p3f),'-r','linewidth',1.5);
hold on;
plot(AltX,VFR(:,p3f),'-b','linewidth',1.5);
hold on;
plot([20 115],[1.4240e+05 1.4240e+05],'--g','linewidth',1);
xlim([20 115])
grid on;
legend('After','Before')
subplot(3,1,3)
plot(AltX,Modify_V_F_L(:,p3f),'-r','linewidth',1.5);
hold on;
plot(AltX,VFL(:,p3f),'-b','linewidth',1.5);
hold on;
plot([20 115],[1.4240e+05 1.4240e+05],'--g','linewidth',1);
xlim([20 115])
grid on;
legend('After','Before')

%% 修改初始光子数，谨慎运行

VF0 = Modify_V_F_0;
VFR = Modify_V_F_R;
VFL = Modify_V_F_L;
NF0 = Modify_N_F_0;
NFR = Modify_N_F_R;
NFL = Modify_N_F_L;
EF0 = Modify_E_F_0;
EFR = Modify_E_F_R;
EFL = Modify_E_F_L;


%% 归一化
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
Temp_Plot = V0;
for jvt = 1:size(V0,2)
    Temp_Plot(:,jvt) = FT(R_V_real_V(:,jvt),R_T_real_V(:,jvt),TruePxT,jvt);
end

syms RV RT;

PreT = Altitude;
PreV = Altitude;
PreTN = Alt30;
PreVN = Alt30;
PreTE = Alt30;
PreVE = Alt30;

Jnight = 1

PhF0 = VF0;
PhFR = VFR;
PhFL = VFL;

AltF = Altitude;

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
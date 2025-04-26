% CapErrorWindTemperature

% 基本参数
TT = 200;                   % 温度，K
VV = 0;                     % 风速，m/s
% DnuC = (600:100:3000)';      % 频移量，MHz
DnuC = 630;                % 频移量，MHz

Plot_TL = DnuC;
Plot_TH = DnuC;
Plot_VL = DnuC;
Plot_VH = DnuC;

for i = 1:size(DnuC,1)

% 不同空间分辨率光子数确定，时间分辨率均为1 h
[N0E1, NRE1, NLE1] = DefCounts(DnuC(i), TT, VV, 12887735);     % E区 1 km
[N0E2, NRE2, NLE2] = DefCounts(DnuC(i), TT, VV, 17270188);     % E区 2 km
[N0E5, NRE5, NLE5] = DefCounts(DnuC(i), TT, VV, 27745289);     % E区 5 km
[N0F5, NRF5, NLF5] = DefCounts(DnuC(i), TT, VV, 17077);        % F区 5 km
[N0F10, NRF10, NLF10] = DefCounts(DnuC(i), TT, VV, 34119);     % F区 10 km
[N0F20, NRF20, NLF20] = DefCounts(DnuC(i), TT, VV, 64382);     % F区 20 km

Hang = 1;
% 计算误差
N0arr = [N0E1;N0E2;N0E5;N0F5;N0F10;N0F20];
NRarr = [NRE1;NRE2;NRE5;NRF5;NRF10;NRF20];
NLarr = [NLE1;NLE2;NLE5;NLF5;NLF10;NLF20];
N0 = N0arr(Hang);
NR = NRarr(Hang);
NL = NLarr(Hang);
RT_real = (NR+NL)./(2*N0);
RV_real = (NR-NL)./N0;

syms RT RV;

% 对 fTL 的 RT 和 RV 分别求偏导
dT_dRTL = diff(fTL(RT,RV),RT);
dT_dRVL = diff(fTL(RT,RV),RV);
dvalueTRTL = subs(dT_dRTL,[RT RV],[RV_real RT_real]);
dvalueTRTL = double(dvalueTRTL);
dvalueTRVL = subs(dT_dRVL,[RT RV],[RV_real RT_real]);
dvalueTRVL = double(dvalueTRVL);

% 对 fTH 的 RT 和 RV 分别求偏导
dT_dRTH = diff(fTH(RT,RV),RT);
dT_dRVH = diff(fTH(RT,RV),RV);
dvalueTRTH = subs(dT_dRTH,[RT RV],[RV_real RT_real]);
dvalueTRTH = double(dvalueTRTH);
dvalueTRVH = subs(dT_dRVH,[RT RV],[RV_real RT_real]);
dvalueTRVH = double(dvalueTRVH);

% 对 fVL 的 RT 和 RV 分别求偏导
dVL_dRT = diff(fVL(RT,RV),RT);
dVL_dRV = diff(fVL(RT,RV),RV);
dvalueVRTL = subs(dVL_dRT,[RT RV],[RV_real RT_real]);
dvalueVRTL = double(dvalueVRTL);
dvalueVRVL = subs(dVL_dRV,[RT RV],[RV_real RT_real]);
dvalueVRVL = double(dvalueVRVL);

% 对 fVH 的 RT 和 RV 分别求偏导
dVH_dRT = diff(fVH(RT,RV),RT);
dVH_dRV = diff(fVH(RT,RV),RV);
dvalueVRTH = subs(dVH_dRT,[RT RV],[RV_real RT_real]);
dvalueVRTH = double(dvalueVRTH);
dvalueVRVH = subs(dVH_dRV,[RT RV],[RV_real RT_real]);
dvalueVRVH = double(dvalueVRVH);

% 定义并计算Pre_TL
Pre_TL = sqrt((NR.*(0.5.*dvalueTRTL+dvalueTRVL).^2+NL.*(0.5.*dvalueTRTL-dvalueTRVL).^2+N0.*(RT_real.*dvalueTRTL+RV_real.*dvalueTRVL).^2)./(N0.^2));
Pre_TL = real(double(Pre_TL));

% 定义并计算Pre_TH
Pre_TH = sqrt((NR.*(0.5.*dvalueTRTH+dvalueTRVH).^2+NL.*(0.5.*dvalueTRTH-dvalueTRVH).^2+N0.*(RT_real.*dvalueTRTH+RV_real.*dvalueTRVH).^2)./(N0.^2));
Pre_TH = real(double(Pre_TH));

% 定义并计算Pre_VL
Pre_VL = sqrt((NR.*(0.5.*dvalueVRTL+dvalueVRVL).^2+NL.*(0.5.*dvalueVRTL-dvalueVRVL).^2+N0.*(RT_real.*dvalueVRTL+RV_real.*dvalueVRVL).^2)./(N0.^2));
Pre_VL = real(double(Pre_VL));

% 定义并计算Pre_VH
Pre_VH = sqrt((NR.*(0.5.*dvalueVRTH+dvalueVRVH).^2+NL.*(0.5.*dvalueVRTH-dvalueVRVH).^2+N0.*(RT_real.*dvalueVRTH+RV_real.*dvalueVRVH).^2)./(N0.^2));
Pre_VH = real(double(Pre_VH));

Plot_TL(i) = Pre_TL;
Plot_TH(i) = Pre_TH;
Plot_VL(i) = Pre_VL;
Plot_VH(i) = Pre_VH;

stri = string(i);
disp(['现在是第',stri]);

end
disp(['TL=',string(Pre_TL)]);
disp(['TH=',string(Pre_TH)]);
disp(['VL=',string(Pre_VL)]);
disp(['VH=',string(Pre_VH)]);

% % 总输出±
% figure('Name','ErrorTem')
% plot(DnuC,Plot_TL,'-r','linewidth',2);
% grid on;
% xlabel('Frequency Offset (MHz)');
% ylabel('Temperature Error (K)');
% TitStr = 'Temperature Error: '+string(TT)+' K, '+string(VV)+' m/s,'+' 1 h, 5 km';
% title(TitStr);
% 
% % 总输出
% figure('Name','ErrorWind')
% plot(DnuC,Plot_VL,'-b','linewidth',2);
% grid on;
% xlabel('Frequency Offset (MHz)');
% ylabel('Wind Velocity Error (m/s)');
% TitStr = 'Wind Velocity Error: '+string(TT)+' K, '+string(VV)+' m/s,'+' 1 h, 5 km';
% title(TitStr);






















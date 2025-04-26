% CapErrorWindTemperature

% ��������
TT = 200;                   % �¶ȣ�K
VV = 0;                     % ���٣�m/s
% DnuC = (600:100:3000)';      % Ƶ������MHz
DnuC = 630;                % Ƶ������MHz

Plot_TL = DnuC;
Plot_TH = DnuC;
Plot_VL = DnuC;
Plot_VH = DnuC;

for i = 1:size(DnuC,1)

% ��ͬ�ռ�ֱ��ʹ�����ȷ����ʱ��ֱ��ʾ�Ϊ1 h
[N0E1, NRE1, NLE1] = DefCounts(DnuC(i), TT, VV, 12887735);     % E�� 1 km
[N0E2, NRE2, NLE2] = DefCounts(DnuC(i), TT, VV, 17270188);     % E�� 2 km
[N0E5, NRE5, NLE5] = DefCounts(DnuC(i), TT, VV, 27745289);     % E�� 5 km
[N0F5, NRF5, NLF5] = DefCounts(DnuC(i), TT, VV, 17077);        % F�� 5 km
[N0F10, NRF10, NLF10] = DefCounts(DnuC(i), TT, VV, 34119);     % F�� 10 km
[N0F20, NRF20, NLF20] = DefCounts(DnuC(i), TT, VV, 64382);     % F�� 20 km

Hang = 1;
% �������
N0arr = [N0E1;N0E2;N0E5;N0F5;N0F10;N0F20];
NRarr = [NRE1;NRE2;NRE5;NRF5;NRF10;NRF20];
NLarr = [NLE1;NLE2;NLE5;NLF5;NLF10;NLF20];
N0 = N0arr(Hang);
NR = NRarr(Hang);
NL = NLarr(Hang);
RT_real = (NR+NL)./(2*N0);
RV_real = (NR-NL)./N0;

syms RT RV;

% �� fTL �� RT �� RV �ֱ���ƫ��
dT_dRTL = diff(fTL(RT,RV),RT);
dT_dRVL = diff(fTL(RT,RV),RV);
dvalueTRTL = subs(dT_dRTL,[RT RV],[RV_real RT_real]);
dvalueTRTL = double(dvalueTRTL);
dvalueTRVL = subs(dT_dRVL,[RT RV],[RV_real RT_real]);
dvalueTRVL = double(dvalueTRVL);

% �� fTH �� RT �� RV �ֱ���ƫ��
dT_dRTH = diff(fTH(RT,RV),RT);
dT_dRVH = diff(fTH(RT,RV),RV);
dvalueTRTH = subs(dT_dRTH,[RT RV],[RV_real RT_real]);
dvalueTRTH = double(dvalueTRTH);
dvalueTRVH = subs(dT_dRVH,[RT RV],[RV_real RT_real]);
dvalueTRVH = double(dvalueTRVH);

% �� fVL �� RT �� RV �ֱ���ƫ��
dVL_dRT = diff(fVL(RT,RV),RT);
dVL_dRV = diff(fVL(RT,RV),RV);
dvalueVRTL = subs(dVL_dRT,[RT RV],[RV_real RT_real]);
dvalueVRTL = double(dvalueVRTL);
dvalueVRVL = subs(dVL_dRV,[RT RV],[RV_real RT_real]);
dvalueVRVL = double(dvalueVRVL);

% �� fVH �� RT �� RV �ֱ���ƫ��
dVH_dRT = diff(fVH(RT,RV),RT);
dVH_dRV = diff(fVH(RT,RV),RV);
dvalueVRTH = subs(dVH_dRT,[RT RV],[RV_real RT_real]);
dvalueVRTH = double(dvalueVRTH);
dvalueVRVH = subs(dVH_dRV,[RT RV],[RV_real RT_real]);
dvalueVRVH = double(dvalueVRVH);

% ���岢����Pre_TL
Pre_TL = sqrt((NR.*(0.5.*dvalueTRTL+dvalueTRVL).^2+NL.*(0.5.*dvalueTRTL-dvalueTRVL).^2+N0.*(RT_real.*dvalueTRTL+RV_real.*dvalueTRVL).^2)./(N0.^2));
Pre_TL = real(double(Pre_TL));

% ���岢����Pre_TH
Pre_TH = sqrt((NR.*(0.5.*dvalueTRTH+dvalueTRVH).^2+NL.*(0.5.*dvalueTRTH-dvalueTRVH).^2+N0.*(RT_real.*dvalueTRTH+RV_real.*dvalueTRVH).^2)./(N0.^2));
Pre_TH = real(double(Pre_TH));

% ���岢����Pre_VL
Pre_VL = sqrt((NR.*(0.5.*dvalueVRTL+dvalueVRVL).^2+NL.*(0.5.*dvalueVRTL-dvalueVRVL).^2+N0.*(RT_real.*dvalueVRTL+RV_real.*dvalueVRVL).^2)./(N0.^2));
Pre_VL = real(double(Pre_VL));

% ���岢����Pre_VH
Pre_VH = sqrt((NR.*(0.5.*dvalueVRTH+dvalueVRVH).^2+NL.*(0.5.*dvalueVRTH-dvalueVRVH).^2+N0.*(RT_real.*dvalueVRTH+RV_real.*dvalueVRVH).^2)./(N0.^2));
Pre_VH = real(double(Pre_VH));

Plot_TL(i) = Pre_TL;
Plot_TH(i) = Pre_TH;
Plot_VL(i) = Pre_VL;
Plot_VH(i) = Pre_VH;

stri = string(i);
disp(['�����ǵ�',stri]);

end
disp(['TL=',string(Pre_TL)]);
disp(['TH=',string(Pre_TH)]);
disp(['VL=',string(Pre_VL)]);
disp(['VH=',string(Pre_VH)]);

% % �������
% figure('Name','ErrorTem')
% plot(DnuC,Plot_TL,'-r','linewidth',2);
% grid on;
% xlabel('Frequency Offset (MHz)');
% ylabel('Temperature Error (K)');
% TitStr = 'Temperature Error: '+string(TT)+' K, '+string(VV)+' m/s,'+' 1 h, 5 km';
% title(TitStr);
% 
% % �����
% figure('Name','ErrorWind')
% plot(DnuC,Plot_VL,'-b','linewidth',2);
% grid on;
% xlabel('Frequency Offset (MHz)');
% ylabel('Wind Velocity Error (m/s)');
% TitStr = 'Wind Velocity Error: '+string(TT)+' K, '+string(VV)+' m/s,'+' 1 h, 5 km';
% title(TitStr);






















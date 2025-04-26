%% �����ӷ���������

% E���¶� 200 K ; F���¶� 900 K
syms RT RV x y

% ��������
TT = 200;                   % �¶ȣ�K
VV = 0;                     % ���٣�m/s
% DnuC = (100:100:3000)';      % Ƶ������MHz
DnuC = 648;                % Ƶ������MHz

Plot_TL = DnuC;
Plot_TH = DnuC;
Plot_VL = DnuC;
Plot_VH = DnuC;

for DD = 1:size(DnuC,1)
    %% E��ɢ�����

    % ͨ�ó���
    c = 2.99792458e8;                           %����(m/s)
    k_B = 1.3806505e-23;                        %������������(J/K)
    e = 1.60217662e-19;                         %���ӵ���(C)
    m_e = 9.10938215e-31;                       %��������(kg)
    epsilon_0 = 8.854187817e-12;                %��յ�����(F/m)

    % ���ⳣ��
    lambda_L = 393.3663e-9;                     %���Ⲩ��(m)
    nu_L = c / lambda_L;                        %��������Ƶ��(Hz)
    sigma_FWHM = 180e6;                         %�����߿�
    sigma_L = sigma_FWHM / (2*sqrt(2*log(2)));  %��˹���͵�RMS���

    % ƫ����
%     DnuR = DnuC(DD);
%     DnuL = -DnuC(DD);
    DnuR = 650;
    DnuL = -650;

    % �����ӳ���
    M = 6.665e-26;                              %�����Ӿ�������(kg)
    fD2 = 0.69;                                 %����������ǿ��

    % �Ϳշ�������
    T = TT-5:0.1:TT+5;
    V = VV-5:0.1:VV+5;

    % ��Ч����ɢ�����ѭ������
    Data_eff_0 = zeros(size(T,2),size(V,2));
    Data_eff_R = zeros(size(T,2),size(V,2));
    Data_eff_L = zeros(size(T,2),size(V,2));

    % �������ѭ��
    for j = 1:length(V)

        nu_0_1 = nu_L * ((c+V(j))/c);
        f1 = 1;
        nu = nu_L-3e9:1e6:nu_L+3e9;
        g_L = (1./(sqrt(2.*pi).*sigma_L)) .* exp(-(nu-nu_L).^2 ./ (2.*sigma_L.^2));

        % �����¶�ѭ��
        for i = 1:length(T)
            sigma_D_1 = nu_0_1 * sqrt(k_B*T(i) / (M*c^2));
            sigma_0_1 = (1 / (sqrt(2*pi)*sigma_D_1)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
            sigma_abs_1 = sigma_0_1 * exp(-(nu-nu_0_1).^2 / (2*sigma_D_1^2));
            sigma_abs_1_fi = sigma_abs_1*f1;
            sigma_abs = sigma_abs_1_fi;
            sigma_eff_fx = conv(sigma_abs,g_L);
            sigma_eff_fx = sigma_eff_fx * 1e6;
            Data_eff_0(i,j) = sigma_eff_fx(6001)/(4*pi);
            Data_eff_R(i,j) = sigma_eff_fx(6001+DnuR)/(4*pi);
            Data_eff_L(i,j) = sigma_eff_fx(6001+DnuL)/(4*pi);
        end
    end

    % ������Ϸ��ٺ��¶�����
    T = T';
    T = repmat(T,1,size(V,2));
    V = repmat(V,size(T,1),1);

    % ��Ϸ��±��ʾ���
    R_T = (Data_eff_R + Data_eff_L) ./ (2 * Data_eff_0);
    R_V = (Data_eff_R - Data_eff_L) ./ (Data_eff_0);

    % ���¶Ƚ��ж�ά�������
    ft = fittype( 'poly55' );
    [fit_T, gof_T] = fit( [R_V(:), R_T(:)], T(:), ft );
    coeTL=coeffvalues(fit_T);
    f_TL = fit_T;

    p00 = coeTL(1);
    p10 = coeTL(2);
    p01 = coeTL(3);
    p20 = coeTL(4);
    p11 = coeTL(5);
    p02 = coeTL(6);
    p30 = coeTL(7);
    p21 = coeTL(8);
    p12 = coeTL(9);
    p03 = coeTL(10);
    p40 = coeTL(11);
    p31 = coeTL(12);
    p22 = coeTL(13);
    p13 = coeTL(14);
    p04 = coeTL(15);
    p50 = coeTL(16);
    p41 = coeTL(17);
    p32 = coeTL(18);
    p23 = coeTL(19);
    p14 = coeTL(20);
    p05 = coeTL(21);

    fFTL = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 ...
        + p21*x^2*y + p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y ...
        + p22*x^2*y^2 + p13*x*y^3 + p04*y^4 + p50*x^5 + p41*x^4*y ...
        + p32*x^3*y^2 + p23*x^2*y^3 + p14*x*y^4 + p05*y^5;
    fFTL = subs(fFTL, {'x', 'y'}, {RT, RV});    % �滻����

    % �Է��ٽ��ж�ά�������
    ft = fittype( 'poly55' );
    [fit_V, gof_V] = fit( [R_V(:), R_T(:)], V(:), ft );
    coeVL=coeffvalues(fit_V);
    f_VL = fit_V;

    p00 = coeVL(1);
    p10 = coeVL(2);
    p01 = coeVL(3);
    p20 = coeVL(4);
    p11 = coeVL(5);
    p02 = coeVL(6);
    p30 = coeVL(7);
    p21 = coeVL(8);
    p12 = coeVL(9);
    p03 = coeVL(10);
    p40 = coeVL(11);
    p31 = coeVL(12);
    p22 = coeVL(13);
    p13 = coeVL(14);
    p04 = coeVL(15);
    p50 = coeVL(16);
    p41 = coeVL(17);
    p32 = coeVL(18);
    p23 = coeVL(19);
    p14 = coeVL(20);
    p05 = coeVL(21);

    fFVL = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 ...
        + p21*x^2*y + p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y ...
        + p22*x^2*y^2 + p13*x*y^3 + p04*y^4 + p50*x^5 + p41*x^4*y ...
        + p32*x^3*y^2 + p23*x^2*y^3 + p14*x*y^4 + p05*y^5;
    fFVL = subs(fFVL, {'x', 'y'}, {RT, RV});    % �滻����

    %% F��ɢ�����-�����������

    % ͨ�ó���
    c = 2.99792458e8;                           %����(m/s)
    k_B = 1.3806505e-23;                        %������������(J/K)
    e = 1.60217662e-19;                         %���ӵ���(C)
    m_e = 9.10938215e-31;                       %��������(kg)
    epsilon_0 = 8.854187817e-12;                %��յ�����(F/m)

    % ���ⳣ��
    lambda_L = 393.477469e-9;                     %���Ⲩ��(m)
    nu_L = c / lambda_L;                        %��������Ƶ��(Hz)
    sigma_FWHM = 180e6;                         %�����߿�
    sigma_L = sigma_FWHM / (2*sqrt(2*log(2)));  %��˹���͵�RMS���

    % ƫ����
    DnuR = 1500;
    DnuL = -1500;

    % �����ӳ���
    M = 6.665e-26;                              %�����Ӿ�������(kg)
    fD2 = 0.69;                                 %����������ǿ��

    % �߿շ�������
    T = TT-5:0.1:TT+5;
    V = VV-5:0.1:VV+5;

    % ��Ч����ɢ�����ѭ������
    Data_eff_0 = zeros(size(T,2),size(V,2));
    Data_eff_R = zeros(size(T,2),size(V,2));
    Data_eff_L = zeros(size(T,2),size(V,2));

    % �������ѭ��
    for j = 1:length(V)

        nu_0_1 = nu_L * ((c+V(j))/c);
        f1 = 1;
        nu = nu_L-3e9:1e6:nu_L+3e9;
        g_L = (1./(sqrt(2.*pi).*sigma_L)) .* exp(-(nu-nu_L).^2 ./ (2.*sigma_L.^2));

        % �����¶�ѭ��
        for i = 1:length(T)
            sigma_D_1 = nu_0_1 * sqrt(k_B*T(i) / (M*c^2));
            sigma_0_1 = (1 / (sqrt(2*pi)*sigma_D_1)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
            sigma_abs_1 = sigma_0_1 * exp(-(nu-nu_0_1).^2 / (2*sigma_D_1^2));
            sigma_abs_1_fi = sigma_abs_1*f1;
            sigma_abs = sigma_abs_1_fi;
            sigma_eff_fx = conv(sigma_abs,g_L);
            sigma_eff_fx = sigma_eff_fx * 1e6;
            Data_eff_0(i,j) = sigma_eff_fx(6001)/(4*pi);
            Data_eff_R(i,j) = sigma_eff_fx(6001+DnuR)/(4*pi);
            Data_eff_L(i,j) = sigma_eff_fx(6001+DnuL)/(4*pi);
        end
    end

    % ������Ϸ��ٺ��¶�����
    T = T';
    T = repmat(T,1,size(V,2));
    V = repmat(V,size(T,1),1);

    % ��Ϸ��±��ʾ���
    R_T = (Data_eff_R + Data_eff_L) ./ (2 * Data_eff_0);
    R_V = (Data_eff_R - Data_eff_L) ./ (Data_eff_0);

    % ���¶Ƚ��ж�ά�������
    ft = fittype( 'poly55' );
    [fit_T, gof_T] = fit( [R_V(:), R_T(:)], T(:), ft );
    coeTH=coeffvalues(fit_T);
    f_TH = fit_T;

    p00 = coeTH(1);
    p10 = coeTH(2);
    p01 = coeTH(3);
    p20 = coeTH(4);
    p11 = coeTH(5);
    p02 = coeTH(6);
    p30 = coeTH(7);
    p21 = coeTH(8);
    p12 = coeTH(9);
    p03 = coeTH(10);
    p40 = coeTH(11);
    p31 = coeTH(12);
    p22 = coeTH(13);
    p13 = coeTH(14);
    p04 = coeTH(15);
    p50 = coeTH(16);
    p41 = coeTH(17);
    p32 = coeTH(18);
    p23 = coeTH(19);
    p14 = coeTH(20);
    p05 = coeTH(21);

    fFTH = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 ...
        + p21*x^2*y + p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y ...
        + p22*x^2*y^2 + p13*x*y^3 + p04*y^4 + p50*x^5 + p41*x^4*y ...
        + p32*x^3*y^2 + p23*x^2*y^3 + p14*x*y^4 + p05*y^5;
    fFTH = subs(fFTH, {'x', 'y'}, {RT, RV});    % �滻����

    % �Է��ٽ��ж�ά�������
    ft = fittype( 'poly55' );
    [fit_V, gof_V] = fit( [R_V(:), R_T(:)], V(:), ft );
    coeVH=coeffvalues(fit_V);
    f_VH = fit_V;

    p00 = coeVH(1);
    p10 = coeVH(2);
    p01 = coeVH(3);
    p20 = coeVH(4);
    p11 = coeVH(5);
    p02 = coeVH(6);
    p30 = coeVH(7);
    p21 = coeVH(8);
    p12 = coeVH(9);
    p03 = coeVH(10);
    p40 = coeVH(11);
    p31 = coeVH(12);
    p22 = coeVH(13);
    p13 = coeVH(14);
    p04 = coeVH(15);
    p50 = coeVH(16);
    p41 = coeVH(17);
    p32 = coeVH(18);
    p23 = coeVH(19);
    p14 = coeVH(20);
    p05 = coeVH(21);

    fFVH = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 ...
        + p21*x^2*y + p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y ...
        + p22*x^2*y^2 + p13*x*y^3 + p04*y^4 + p50*x^5 + p41*x^4*y ...
        + p32*x^3*y^2 + p23*x^2*y^3 + p14*x*y^4 + p05*y^5;
    fFVH = subs(fFVH, {'x', 'y'}, {RT, RV});    % �滻����


    %% CapErrorWindTemperature

    % ��ͬ�ռ�ֱ��ʹ�����ȷ����ʱ��ֱ��ʾ�Ϊ1 h
    
    [N0E1, NRE1, NLE1] = DefCounts(DnuC(DD), TT, VV, 150000);     % E�� 1 km
    
%     [N0E1, NRE1, NLE1] = DefCounts(DnuC(DD), TT, VV, 12887735);     % E�� 1 km
    [N0E2, NRE2, NLE2] = DefCounts(DnuC(DD), TT, VV, 17270188);     % E�� 2 km
    [N0E5, NRE5, NLE5] = DefCounts(DnuC(DD), TT, VV, 27745289);     % E�� 5 km
    [N0F5, NRF5, NLF5] = DefCounts(DnuC(DD), TT, VV, 17077);        % F�� 5 km
    [N0F10, NRF10, NLF10] = DefCounts(DnuC(DD), TT, VV, 34119);     % F�� 10 km
    [N0F20, NRF20, NLF20] = DefCounts(DnuC(DD), TT, VV, 64382);     % F�� 20 km

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

    % �� fTL �� RT �� RV �ֱ���ƫ��
    dfTL_dRT = diff(fFTL, RT);
    dfTL_dRV = diff(fFTL, RV);
    dvalueTRTL = subs(dfTL_dRT, {RT, RV}, {RT_real,RV_real});
    dvalueTRTL = double(dvalueTRTL);
    dvalueTRVL = subs(dfTL_dRV, {RT, RV}, {RT_real,RV_real});
    dvalueTRVL = double(dvalueTRVL);

    % �� fTH �� RT �� RV �ֱ���ƫ��
    dfTH_dRT = diff(fFTH, RT);
    dfTH_dRV = diff(fFTH, RV);
    dvalueTRTH = subs(dfTH_dRT, {RT, RV}, {RT_real,RV_real});
    dvalueTRTH = double(dvalueTRTH);
    dvalueTRVH = subs(dfTH_dRV, {RT, RV}, {RT_real,RV_real});
    dvalueTRVH = double(dvalueTRVH);

    % �� fVL �� RT �� RV �ֱ���ƫ��
    dfVL_dRT = diff(fFVL, RT);
    dfVL_dRV = diff(fFVL, RV);
    dvalueVRTL = subs(dfVL_dRT, {RT, RV}, {RT_real,RV_real});
    dvalueVRTL = double(dvalueVRTL);
    dvalueVRVL = subs(dfVL_dRV, {RT, RV}, {RT_real,RV_real});
    dvalueVRVL = double(dvalueVRVL);

    % �� fVH �� RT �� RV �ֱ���ƫ��
    dfVH_dRT = diff(fFVH, RT);
    dfVH_dRV = diff(fFVH, RV);
    dvalueVRTH = subs(dfVH_dRT, {RT, RV}, {RT_real,RV_real});
    dvalueVRTH = double(dvalueVRTH);
    dvalueVRVH = subs(dfVH_dRV, {RT, RV}, {RT_real,RV_real});
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

    Plot_TL(DD) = Pre_TL;
    Plot_TH(DD) = Pre_TH;
    Plot_VL(DD) = Pre_VL;
    Plot_VH(DD) = Pre_VH;

    strDD = string(DD);
    disp(['�����ǵ�',strDD]);

end

disp(['TL=',string(Pre_TL)]);
disp(['TH=',string(Pre_TH)]);
disp(['VL=',string(Pre_VL)]);
disp(['VH=',string(Pre_VH)]);

% % �������
% figure('Name','ErrorTem')
% plot(DnuC,Plot_TH,'-r','linewidth',2);
% grid on;
% xlabel('Frequency Offset (MHz)');
% ylabel('Temperature Error (K)');
% TitStr = 'Temperature Error: '+string(TT)+' K, '+string(VV)+' m/s,'+' 1 h, 1 km';
% title(TitStr);
% xlim([0 3000])
% ylim([0 30])
% 
% % �����
% figure('Name','ErrorWind')
% plot(DnuC,Plot_VH,'-b','linewidth',2);
% grid on;
% xlabel('Frequency Offset (MHz)');
% ylabel('Wind Velocity Error (m/s)');
% TitStr = 'Wind Velocity Error: '+string(TT)+' K, '+string(VV)+' m/s,'+' 1 h, 1 km';
% title(TitStr);
% xlim([0 3000])
% ylim([0 30])













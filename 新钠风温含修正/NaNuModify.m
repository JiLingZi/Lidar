%% ���ݷ糡
clear all

% ͨ�ó���
c = 2.99792458e8;                           %����(m/s)
k_B = 1.3806505e-23;                        %������������(J/K)
e = 1.60217662e-19;                         %���ӵ���(C)
m_e = 9.10938215e-31;                       %��������(kg)
epsilon_0 = 8.854187817e-12;                %��յ�����(F/m)
% ���ⳣ��
lambda_L = 589.1583e-9;                     %���Ⲩ��(m)
nu_L = c / lambda_L;                        %��������Ƶ��(Hz)
sigma_FWHM = 80e6;                          %�����߿�
sigma_L = sigma_FWHM / (2*sqrt(2*log(2)));  %��˹���͵�RMS���
% ��ԭ�ӳ���
M = 3.82e-26;                               %Naԭ�Ӿ�������(kg)
fD2 = 0.6408;                               %Naԭ��D2�ߵ�����ǿ��
%����������
% % �¶�����
% T = [80 100 120 140 160 180 200 220 240 260 280 300 320];
% % ��������
% V = [-120 -100 -80 -60 -40 -20 0 20 40 60 80 100 120];
T = 80:20:320;
V = -200:20:200;

%��Ч����ɢ��������ѭ��
Data_eff_0 = zeros(size(T,2),size(V,2));
Data_eff_R = zeros(size(T,2),size(V,2));
Data_eff_L = zeros(size(T,2),size(V,2));
% �������ѭ��
Dnu = -300:1:300;
PxT = zeros(21,601);
PxV = zeros(21,601);
for aa = 1:size(Dnu,2)
    %����������
    % �¶�����
%     T = [80 100 120 140 160 180 200 220 240 260 280 300 320];
    T = 80:20:320;
    % ��������
%     V = [-120 -100 -80 -60 -40 -20 0 20 40 60 80 100 120];
    V = -200:20:200;
    for j = 1:length(V)    
        %���߲�����ֵ��
        % Na��6��ԾǨ���ߵ�����Ƶ��ƫ��ֵ(Hz)
        nu01 = (nu_L - 0.7328e9 + 0.6302e9) * ((c+V(j))/c);
        nu02 = (nu_L - 0.6962e9 + 0.6302e9) * ((c+V(j))/c);
        nu03 = (nu_L - 0.6302e9 + 0.6302e9) * ((c+V(j))/c);
        nu04 = (nu_L + 1.0333e9 + 0.6302e9) * ((c+V(j))/c);
        nu05 = (nu_L + 1.0552e9 + 0.6302e9) * ((c+V(j))/c);
        nu06 = (nu_L + 1.0919e9 + 0.6302e9) * ((c+V(j))/c);
        % Na�ĳ���ϸԾǨǿ������
        f1 = 1/32;
        f2 = 5/32;
        f3 = 14/32;
        f4 = 2/32;
        f5 = 5/32;
        f6 = 5/32;
        %��һ������
        nu = nu_L-3e9:1e6:nu_L+3e9;
        g_L = (1./(sqrt(2.*pi).*sigma_L)) .* exp(-(nu-nu_L).^2 ./ (2.*sigma_L.^2));
        % �¶�ѭ��
        for i = 1:length(T)
            %�����������RMS���
            sigma_D_1 = nu01 * sqrt(k_B*T(i) / (M*c^2));
            sigma_D_2 = nu02 * sqrt(k_B*T(i) / (M*c^2));
            sigma_D_3 = nu03 * sqrt(k_B*T(i) / (M*c^2));
            sigma_D_4 = nu04 * sqrt(k_B*T(i) / (M*c^2));
            sigma_D_5 = nu05 * sqrt(k_B*T(i) / (M*c^2));
            sigma_D_6 = nu06 * sqrt(k_B*T(i) / (M*c^2));
            %��ֵ�������ս���
            sigma_0_1 = (1 / (sqrt(2*pi)*sigma_D_1)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
            sigma_0_2 = (1 / (sqrt(2*pi)*sigma_D_2)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
            sigma_0_3 = (1 / (sqrt(2*pi)*sigma_D_3)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
            sigma_0_4 = (1 / (sqrt(2*pi)*sigma_D_4)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
            sigma_0_5 = (1 / (sqrt(2*pi)*sigma_D_5)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
            sigma_0_6 = (1 / (sqrt(2*pi)*sigma_D_6)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
            %����ԾǨ�ߵĵ������ս���
            sigma_abs_1 = sigma_0_1 * exp(-(nu-nu01).^2 / (2*sigma_D_1^2));
            sigma_abs_2 = sigma_0_2 * exp(-(nu-nu02).^2 / (2*sigma_D_2^2));
            sigma_abs_3 = sigma_0_3 * exp(-(nu-nu03).^2 / (2*sigma_D_3^2));
            sigma_abs_4 = sigma_0_4 * exp(-(nu-nu04).^2 / (2*sigma_D_4^2));
            sigma_abs_5 = sigma_0_5 * exp(-(nu-nu05).^2 / (2*sigma_D_5^2));
            sigma_abs_6 = sigma_0_6 * exp(-(nu-nu06).^2 / (2*sigma_D_6^2));
            %�������ս�����Գ���ϸԾǨǿ������
            sigma_abs_1_fi = sigma_abs_1*f1;
            sigma_abs_2_fi = sigma_abs_2*f2;
            sigma_abs_3_fi = sigma_abs_3*f3;
            sigma_abs_4_fi = sigma_abs_4*f4;
            sigma_abs_5_fi = sigma_abs_5*f5;
            sigma_abs_6_fi = sigma_abs_6*f6;
            %���ս������
            sigma_abs = sigma_abs_1_fi + sigma_abs_2_fi + sigma_abs_3_fi + sigma_abs_4_fi + sigma_abs_5_fi + sigma_abs_6_fi;
            %��Ч����ɢ�����
            sigma_eff_fx = conv(sigma_abs,g_L);
            sigma_eff_fx = sigma_eff_fx * 1e6;
            Data_eff_0(i,j) = sigma_eff_fx(6001+Dnu(aa))/(4*pi);
            Data_eff_R(i,j) = sigma_eff_fx(6001+Dnu(aa)+585)/(4*pi);
            Data_eff_L(i,j) = sigma_eff_fx(6001+Dnu(aa)-585)/(4*pi);
        end
    end
    % ������Ϸ��ٺ��¶�����
    T = T';
    T = repmat(T,1,size(V,2));
    V = repmat(V,size(T,1),1);
    % ���±��ʾ���
    R_T = (Data_eff_R + Data_eff_L) ./ (2 * Data_eff_0);
    R_V = (Data_eff_R - Data_eff_L) ./ (Data_eff_0);
    % ���¶Ƚ��ж�ά�������
    ft = fittype( 'poly55' );
    [fit_T, gof_T] = fit( [R_V(:), R_T(:)], T(:), ft );
    PxT(:,aa)=coeffvalues(fit_T)';
    f_T = fit_T;
    % �Է��ٽ��ж�ά�������
    ft = fittype( 'poly55' );
    [fit_V, gof_V] = fit( [R_V(:), R_T(:)], V(:), ft );
    PxV(:,aa)=coeffvalues(fit_V)';
    f_V = fit_V;
    
    READ = '�����ǵ� ' + string(aa) + ' ��ѭ��'
end


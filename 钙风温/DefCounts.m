function [N00, NRR, NLL] = DefCounts(DnuInput, TT, VV, CountInput)

% 光子数假设

% 通用常数
c = 2.99792458e8;                           %光速(m/s)
k_B = 1.3806505e-23;                        %玻尔兹曼常数(J/K)
e = 1.60217662e-19;                         %电子电量(C)
m_e = 9.10938215e-31;                       %电子质量(kg)
epsilon_0 = 8.854187817e-12;                %真空电容率(F/m)

% 激光常数
lambda_L = 393.3663e-9;                     %激光波长(m)
nu_L = c / lambda_L;                        %激光中心频率(Hz)
sigma_FWHM = 180e6;                         %激光线宽
sigma_L = sigma_FWHM / (2*sqrt(2*log(2)));  %高斯线型的RMS宽度

% 偏移量
% DnuR = DnuInput;
% DnuL = -DnuInput;

% 钙离子常数
M = 6.665e-26;                              %钙离子绝对质量(kg)
fD2 = 0.69;                                 %钙离子振子强度

% 风温数据
% T = TT;
% V = VV;

T = 900;
V = 0;

% 有效后向散射截面
nu_0_1 = nu_L * ((c+V)/c);
f1 = 1;
nu = nu_L-3e10:1e6:nu_L+3e10;
g_L = (1./(sqrt(2.*pi).*sigma_L)) .* exp(-(nu-nu_L).^2 ./ (2.*sigma_L.^2));
sigma_D_1 = nu_0_1 * sqrt(k_B*T / (M*c^2));
sigma_0_1 = (1 / (sqrt(2*pi)*sigma_D_1)) * (e^2 / (4*epsilon_0*m_e*c)) * fD2;
sigma_abs_1 = sigma_0_1 * exp(-(nu-nu_0_1).^2 / (2*sigma_D_1^2));
sigma_abs_1_fi = sigma_abs_1*f1;
sigma_abs = sigma_abs_1_fi;
sigma_eff_fx = conv(sigma_abs,g_L);
sigma_eff_fx = sigma_eff_fx * 1e6;
Data_eff_0 = sigma_eff_fx(60001)/(4*pi);
Data_eff_R = sigma_eff_fx(60001+DnuR)/(4*pi);
Data_eff_L = sigma_eff_fx(60001+DnuL)/(4*pi);

% 假定光子数

N00 = CountInput;
NRR = CountInput.*(Data_eff_R./Data_eff_0);
NLL = CountInput.*(Data_eff_L./Data_eff_0);
end



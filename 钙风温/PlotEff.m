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
% 钙离子常数
M = 6.665e-26;                              %钙离子绝对质量(kg)
fD2 = 0.69;                                 %钙离子振子强度

% % 设置帧的宽度和高度
% figurePosition = get(gcf, 'Position');
% filename1 = 'ScatterFWHM.mp4';
% writerObj1 = VideoWriter(filename1, 'MPEG-4');
% writerObj1.FrameRate = 10; % 设置帧速率为每秒2帧
% open(writerObj1);

% Tcl = 180:5:1080;
Tcl = 900;
for i = 1:size(Tcl,2)
% 风温数据
T = Tcl(i);
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
max(sigma_eff_fx)



% 散射截面画图
Eff = sigma_eff_fx(30000:90000)';
Nu = nu-nu_L';
plot(Nu,Eff,'-b','linewidth',2);
xlabel('Frequency Offset(Hz)');
ylabel('Cross-section (m^{2})');
xlim([-3e9 3e9]);
ylim([0 1.5e-15])
grid on;
hold on;
HMindex = find(Eff > 0.5*max(Eff), 1);
HM = Eff(HMindex);
FW1 = Nu(HMindex);
FW2 = -Nu(HMindex);
plot([FW1,FW2],[HM,HM],'-r','linewidth',2)
legstr = string(T)+' K, '+string(max(sigma_eff_fx))+' m/s';
legstr2 = '0.5*FWHM: '+string(FW2/1e6)+' MHz';
legend(legstr,legstr2);
title('钙离子散射截面随温度风速变化');

% % 捕捉figure1内容并写入VideoWriter对象1
% frame1 = getframe;
% writeVideo(writerObj1, frame1);
% cla;
% 
% strDD = string(i);
% disp(['现在是第',strDD]);

end
% 
% % 关闭VideoWriter对象1
% close(writerObj1);



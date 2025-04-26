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
% �����ӳ���
M = 6.665e-26;                              %�����Ӿ�������(kg)
fD2 = 0.69;                                 %����������ǿ��

% % ����֡�Ŀ�Ⱥ͸߶�
% figurePosition = get(gcf, 'Position');
% filename1 = 'ScatterFWHM.mp4';
% writerObj1 = VideoWriter(filename1, 'MPEG-4');
% writerObj1.FrameRate = 10; % ����֡����Ϊÿ��2֡
% open(writerObj1);

% Tcl = 180:5:1080;
Tcl = 900;
for i = 1:size(Tcl,2)
% ��������
T = Tcl(i);
V = 0;

% ��Ч����ɢ�����
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



% ɢ����滭ͼ
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
title('������ɢ��������¶ȷ��ٱ仯');

% % ��׽figure1���ݲ�д��VideoWriter����1
% frame1 = getframe;
% writeVideo(writerObj1, frame1);
% cla;
% 
% strDD = string(i);
% disp(['�����ǵ�',strDD]);

end
% 
% % �ر�VideoWriter����1
% close(writerObj1);



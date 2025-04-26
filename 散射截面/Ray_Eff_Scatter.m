% 30 km ������������ɢ�����ϵ��
ray_k = 4.83e-21;

% �����Ӳ���
lambda = [341.4764 , 371.9937 , 393.3663 , 422.6727 , 589.1583 , 769.898];

% ����������� 30 km ������������ɢ�����
ray_sct_cs = ray_k .* (lambda.^-4);

% �Ѿ�����õĸ�������Ч����ɢ�����
eff_sct_cs = [3.770096e-18;3.011306e-18;1.090969e-16;2.960050e-16;7.350548e-17;3.828095e-17];
format shortE
disp(['Ni   ',num2str(ray_sct_cs(1)),',',num2str(eff_sct_cs(1))]);
disp(['Fe   ',num2str(ray_sct_cs(2)),',',num2str(eff_sct_cs(2))]);
disp(['Ca^+ ',num2str(ray_sct_cs(3)),',',num2str(eff_sct_cs(3))]);
disp(['Ca   ',num2str(ray_sct_cs(4)),',',num2str(eff_sct_cs(4))]);
disp(['Na   ',num2str(ray_sct_cs(5)),',',num2str(eff_sct_cs(5))]);
disp(['K    ',num2str(ray_sct_cs(6)),',',num2str(eff_sct_cs(6))]);

% ���ļ�
fid = fopen('RayEffScatter.txt', 'w');
% ִ��disp��������������ݱ��浽�ļ���
fprintf(fid,'# Ray&Eff Scattering cross sections of 6 Channel\n');
fprintf(fid, '%-8s %+11s %+11s\n', 'Channel', 'RayScatter', 'EffScatter');
fprintf(fid,'Ni       %+11s %+11s\n',num2str(ray_sct_cs(1)),num2str(eff_sct_cs(1)));
fprintf(fid,'Fe       %+11s %+11s\n',num2str(ray_sct_cs(2)),num2str(eff_sct_cs(2)));
fprintf(fid,'Cap     %+11s %+11s\n',num2str(ray_sct_cs(3)),num2str(eff_sct_cs(3)));
fprintf(fid,'Ca       %+11s %+11s\n',num2str(ray_sct_cs(4)),num2str(eff_sct_cs(4)));
fprintf(fid,'Na       %+11s %+11s\n',num2str(ray_sct_cs(5)),num2str(eff_sct_cs(5)));
fprintf(fid,'K         %+11s %+11s\n',num2str(ray_sct_cs(6)),num2str(eff_sct_cs(6)));
% �ر��ļ�
fclose(fid);

ScatterTable = readtable('RayEffScatter.txt');
Scatter = table2array(ScatterTable(:,2:3));

%% ����ʦ�����㷨
R = 8.314;                  % ����ѧ���� J/(K*mol)
NA = 6.02214076e23;         % ����٤���޳���
Z_0 = 30000;                % �ο��߶�,m
% lambda = 393.3663e-9;       % ����,m
lambda = 300e-9:1e-9:800e-9;% ����,m
% �������ģ���¶�
[TemRay,DenRay] = atmosnrlmsise00(Z_0,40.3,112.7,2023,281,17);
T_Z_0 = TemRay(1,2);
P_Z_0 = (sum(DenRay(1,:))/NA)*R*T_Z_0*0.01;
RayScat = ((2.983e-32)*(P_Z_0/T_Z_0)./(lambda.^(4.0117)))./sum(DenRay(1,:));
plot(lambda,RayScat,'linewidth',2)



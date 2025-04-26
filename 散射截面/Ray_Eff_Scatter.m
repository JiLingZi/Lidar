% 30 km 处的瑞利后向散射截面系数
ray_k = 4.83e-21;

% 各粒子波长
lambda = [341.4764 , 371.9937 , 393.3663 , 422.6727 , 589.1583 , 769.898];

% 计算各粒子在 30 km 处的瑞利后向散射截面
ray_sct_cs = ray_k .* (lambda.^-4);

% 已经计算好的各粒子有效后向散射截面
eff_sct_cs = [3.770096e-18;3.011306e-18;1.090969e-16;2.960050e-16;7.350548e-17;3.828095e-17];
format shortE
disp(['Ni   ',num2str(ray_sct_cs(1)),',',num2str(eff_sct_cs(1))]);
disp(['Fe   ',num2str(ray_sct_cs(2)),',',num2str(eff_sct_cs(2))]);
disp(['Ca^+ ',num2str(ray_sct_cs(3)),',',num2str(eff_sct_cs(3))]);
disp(['Ca   ',num2str(ray_sct_cs(4)),',',num2str(eff_sct_cs(4))]);
disp(['Na   ',num2str(ray_sct_cs(5)),',',num2str(eff_sct_cs(5))]);
disp(['K    ',num2str(ray_sct_cs(6)),',',num2str(eff_sct_cs(6))]);

% 打开文件
fid = fopen('RayEffScatter.txt', 'w');
% 执行disp函数并将输出内容保存到文件中
fprintf(fid,'# Ray&Eff Scattering cross sections of 6 Channel\n');
fprintf(fid, '%-8s %+11s %+11s\n', 'Channel', 'RayScatter', 'EffScatter');
fprintf(fid,'Ni       %+11s %+11s\n',num2str(ray_sct_cs(1)),num2str(eff_sct_cs(1)));
fprintf(fid,'Fe       %+11s %+11s\n',num2str(ray_sct_cs(2)),num2str(eff_sct_cs(2)));
fprintf(fid,'Cap     %+11s %+11s\n',num2str(ray_sct_cs(3)),num2str(eff_sct_cs(3)));
fprintf(fid,'Ca       %+11s %+11s\n',num2str(ray_sct_cs(4)),num2str(eff_sct_cs(4)));
fprintf(fid,'Na       %+11s %+11s\n',num2str(ray_sct_cs(5)),num2str(eff_sct_cs(5)));
fprintf(fid,'K         %+11s %+11s\n',num2str(ray_sct_cs(6)),num2str(eff_sct_cs(6)));
% 关闭文件
fclose(fid);

ScatterTable = readtable('RayEffScatter.txt');
Scatter = table2array(ScatterTable(:,2:3));

%% 初老师瑞利算法
R = 8.314;                  % 热力学常数 J/(K*mol)
NA = 6.02214076e23;         % 阿伏伽德罗常数
Z_0 = 30000;                % 参考高度,m
% lambda = 393.3663e-9;       % 波长,m
lambda = 300e-9:1e-9:800e-9;% 波长,m
% 导入大气模型温度
[TemRay,DenRay] = atmosnrlmsise00(Z_0,40.3,112.7,2023,281,17);
T_Z_0 = TemRay(1,2);
P_Z_0 = (sum(DenRay(1,:))/NA)*R*T_Z_0*0.01;
RayScat = ((2.983e-32)*(P_Z_0/T_Z_0)./(lambda.^(4.0117)))./sum(DenRay(1,:));
plot(lambda,RayScat,'linewidth',2)



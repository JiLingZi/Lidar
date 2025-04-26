% 考虑大气衰减和钠原子层的吸收，计算温度和风速，
% 改变脉冲数，分别计算

clc;
clear all;

% 给出数据所在路径，共有三个方向，每一个方向都有三个频率
path='E:\temperature_wind_calculate\data\20160803\';
pulse_number='30';

path1=strcat(path,pulse_number,'\Central\*.dat');
path2=strcat(path,pulse_number,'\East\*.dat');
path3=strcat(path,pulse_number,'\North\*.dat');

file_central=dir(path1);
file_east=dir(path2);
file_north=dir(path3);

k1=length(file_central);
k2=length(file_east);
k3=length(file_north);

% 计算垂直方向--------------------------------------------------------------
photo_number_central=zeros(4096,4);
for i=1:k1
    file_name=strcat(path,pulse_number,'\Central\',file_central(i).name);
    data_file_central=importdata(file_name);

    % 将数据合并成一条曲线
    photo_number_central=photo_number_central+data_file_central.data;
end

% 画出一段时间合并后的光子数曲线
dot_origin=4096; % 原始点数4096
bin_origin=90; % 空间分辨率90m
height_origin=(1:dot_origin).*bin_origin./1000;

figure(1);
semilogy(height_origin,photo_number_central(:,2),'k',height_origin,photo_number_central(:,3),'r',height_origin,photo_number_central(:,4),'b');
legend('F0','F+','F-');
picture_name=strcat('E:\temperature_wind_calculate\picture_final\',pulse_number,'_central_profile_origin');
saveas(gcf,picture_name,'tif');

% 画出n点合并后的曲线
dot_number=10;
max_number=dot_origin-mod(dot_origin,dot_number);
% 进行n点合并
photo_merge_central=reshape(sum(reshape(photo_number_central(1:max_number,:),dot_number,[])),[],4);
photo_middle_central=photo_merge_central(:,2);
photo_plus_central=photo_merge_central(:,3);
photo_minus_central=photo_merge_central(:,4);

dot_merge=floor(dot_origin./dot_number);
height_merge=(1:dot_merge).*bin_origin.*dot_number./1000;
height_merge_slant=height_merge.*sqrt(3)./2;

figure(2);
semilogy(height_merge,photo_middle_central,'k',height_merge,photo_plus_central,'r',height_merge,photo_minus_central,'b');
legend('F0','F+','F-');
picture_name=strcat('E:\temperature_wind_calculate\picture_final\',pulse_number,'_central_profile_merge');
saveas(gcf,picture_name,'tif');

% 先扣除噪声（150-200km），再与30km处相比较
h_30=34;
h_200=222;
h_250=278;

% 计算光子数噪声
photo_noise_middle=mean(photo_middle_central(h_200:h_250));
photo_noise_plus=mean(photo_plus_central(h_200:h_250));
photo_noise_minus=mean(photo_minus_central(h_200:h_250));

photo_middle_central=photo_middle_central-photo_noise_middle;
photo_plus_central=photo_plus_central-photo_noise_plus;
photo_minus_central=photo_minus_central-photo_noise_minus;

photo_middle_30_central=photo_middle_central(h_30);
photo_plus_30_central=photo_plus_central(h_30);
photo_minus_30_central=photo_minus_central(h_30);

photo_middle_cali_central=photo_middle_central./photo_middle_30_central;
photo_plus_cali_central=photo_plus_central./photo_plus_30_central;
photo_minus_cali_central=photo_minus_central./photo_minus_30_central;


% 读取拟合系数
filename1=strcat('E:\temperature_wind_calculate\coefficient_fit_temperature_wind.txt');
data_file1=importdata(filename1);
coe_fit_t=data_file1.data(:,1);
coe_fit_w=data_file1.data(:,2);

% 给出拟合函数类型
fit_type=fittype('poly55');
% 给出温度拟合函数
string_1='sfit(fit_type,';
for i=1:20
    string_1=strcat(string_1,num2str(coe_fit_t(i),20),',');
end
string_1=strcat(string_1,num2str(coe_fit_t(21),20),')');
fun_fit_t=eval(string_1);

%给出风速拟合函数
string_2='sfit(fit_type,';
for i=1:20
    string_2=strcat(string_2,num2str(coe_fit_w(i),20),',');
end
string_2=strcat(string_2,num2str(coe_fit_w(21),20),')');
fun_fit_w=eval(string_2);

wavelength_0=589.158e-9;  % D2线的中心波长
lightspeed=299792458;
frequency_0=lightspeed/wavelength_0;  % D2线的中心频率
mass_sodium=3.82e-26;
boltzmann_constant=1.3806488e-23;

% D2线的六种跃迁频率
frequency_5_1=1.0408e9;  
frequency_6_1=1.0566e9;
frequency_7_1=1.0911e9;
frequency_6_2=-0.7150e9;
frequency_7_2=-0.6806e9;
frequency_8_2=-0.6216e9;

alpha_0=1;
alpha_1=5.625/6;
alpha_2=1.175/2;
alpha_3=0.975/6;
alpha_4=0.825/2;
alpha_5=1.005/6;
alpha_6=1.175/2;
alpha_7=1.12;
alpha_8=4.875/6;
alpha_9=0.825/2;

frequency_laser_middle=-0.6541e9;
frequency_laser_plus=frequency_laser_middle+585e6;
frequency_laser_minus=frequency_laser_middle-585e6;

x_middle=frequency_0-0.6541e9;
x_minus=x_middle-585e6;
x_plus=x_middle+585e6;

% 激光与原子相互作用后，原子后向发射光的频率（10个路径，8个频率）
emission_0_middle=2*frequency_5_1-frequency_laser_middle;
emission_1_middle=2*frequency_6_1-frequency_laser_middle;
emission_2_middle=2*frequency_7_1-frequency_laser_middle;
emission_3_middle=frequency_6_1+frequency_6_2-frequency_laser_middle;
emission_4_middle=frequency_7_1+frequency_7_2-frequency_laser_middle;
emission_5_middle=2*frequency_6_2-frequency_laser_middle;
emission_6_middle=2*frequency_7_2-frequency_laser_middle;
emission_7_middle=2*frequency_8_2-frequency_laser_middle;
emission_8_middle=frequency_6_2+frequency_6_1-frequency_laser_middle;
emission_9_middle=frequency_7_2+frequency_7_1-frequency_laser_middle;

emission_0_plus=2*frequency_5_1-frequency_laser_plus;
emission_1_plus=2*frequency_6_1-frequency_laser_plus;
emission_2_plus=2*frequency_7_1-frequency_laser_plus;
emission_3_plus=frequency_6_1+frequency_6_2-frequency_laser_plus;
emission_4_plus=frequency_7_1+frequency_7_2-frequency_laser_plus;
emission_5_plus=2*frequency_6_2-frequency_laser_plus;
emission_6_plus=2*frequency_7_2-frequency_laser_plus;
emission_7_plus=2*frequency_8_2-frequency_laser_plus;
emission_8_plus=frequency_6_2+frequency_6_1-frequency_laser_plus;
emission_9_plus=frequency_7_2+frequency_7_1-frequency_laser_plus;

emission_0_minus=2*frequency_5_1-frequency_laser_minus;
emission_1_minus=2*frequency_6_1-frequency_laser_minus;
emission_2_minus=2*frequency_7_1-frequency_laser_minus;
emission_3_minus=frequency_6_1+frequency_6_2-frequency_laser_minus;
emission_4_minus=frequency_7_1+frequency_7_2-frequency_laser_minus;
emission_5_minus=2*frequency_6_2-frequency_laser_minus;
emission_6_minus=2*frequency_7_2-frequency_laser_minus;
emission_7_minus=2*frequency_8_2-frequency_laser_minus;
emission_8_minus=frequency_6_2+frequency_6_1-frequency_laser_minus;
emission_9_minus=frequency_7_2+frequency_7_1-frequency_laser_minus;

integral_up_replace_middle=0;
integral_up_replace_plus=0;
integral_up_replace_minus=0;

integral_down_replace_middle_0=0;
integral_down_replace_middle_1=0;
integral_down_replace_middle_2=0;
integral_down_replace_middle_3=0;
integral_down_replace_middle_4=0;
integral_down_replace_middle_5=0;
integral_down_replace_middle_6=0;
integral_down_replace_middle_7=0;
integral_down_replace_middle_8=0;
integral_down_replace_middle_9=0;

integral_down_replace_plus_0=0;
integral_down_replace_plus_1=0;
integral_down_replace_plus_2=0;
integral_down_replace_plus_3=0;
integral_down_replace_plus_4=0;
integral_down_replace_plus_5=0;
integral_down_replace_plus_6=0;
integral_down_replace_plus_7=0;
integral_down_replace_plus_8=0;
integral_down_replace_plus_9=0;

integral_down_replace_minus_0=0;
integral_down_replace_minus_1=0;
integral_down_replace_minus_2=0;
integral_down_replace_minus_3=0;
integral_down_replace_minus_4=0;
integral_down_replace_minus_5=0;
integral_down_replace_minus_6=0;
integral_down_replace_minus_7=0;
integral_down_replace_minus_8=0;
integral_down_replace_minus_9=0;

h_bottom=95;
h_top=116;
height_resolution=bin_origin*dot_number;

temperature_central=zeros(h_top,1);
wind_central=zeros(h_top,1);

attenuation_up_middle=ones(h_top,1);
attenuation_down_effective_middle=ones(h_top,1);
attenuation_up_plus=ones(h_top,1);
attenuation_down_effective_plus=ones(h_top,1);
attenuation_up_minus=ones(h_top,1);
attenuation_down_effective_minus=ones(h_top,1);

% 读取大气密度数据
filename2=strcat('E:\temperature_wind_calculate\atmosphere_density.txt');
data_file2=importdata(filename2);
atmosphere_density=data_file2.data(:,2).*1e3;
atmosphere_density=atmosphere_density./(29.*1e-3).*6.022e23;
extinction_cross_section_rayleigh=3.57e-31;
backscattering_cross_section_rayleigh=4.06e-32;
integral_atmosphere_replace=sum(atmosphere_density(h_30:h_bottom-1))*extinction_cross_section_rayleigh*height_resolution;


for i=h_bottom:h_top
    j=i
    % 本层的光子数除以本层一下的向上衰减和向下衰减
    photo_middle_real_central(i)=photo_middle_cali_central(i)./attenuation_up_middle(i-1)./attenuation_down_effective_middle(i-1);
    photo_plus_real_central(i)=photo_plus_cali_central(i)./attenuation_up_plus(i-1)./attenuation_down_effective_plus(i-1);
    photo_minus_real_central(i)=photo_minus_cali_central(i)./attenuation_up_minus(i-1)./attenuation_down_effective_minus(i-1);
    % 利用真实的光子数计算比例系数
    ratio_t_central(i)=(photo_plus_real_central(i)+photo_minus_real_central(i))./(2.*photo_middle_real_central(i));
    ratio_w_central(i)=(photo_plus_real_central(i)-photo_minus_real_central(i))./photo_middle_real_central(i);
    % 计算本层的温度和风速
    temperature_central(i)=feval(fun_fit_t,ratio_w_central(i),ratio_t_central(i));
    wind_central(i)=feval(fun_fit_w,ratio_w_central(i),ratio_t_central(i));
    % 计算本层的温度和风速误差
    
% 计算拟合曲面的偏导数 t=temperature, v=wind
[dT_dRv, dT_dRt]=differentiate(fun_fit_t,ratio_w_central(i),ratio_t_central(i));
[dV_dRv, dV_dRt]=differentiate(fun_fit_w,ratio_w_central(i),ratio_t_central(i));

% 计算三个频率的delta_N^2
delta_photo_middle_square=((photo_middle_central(i)+photo_noise_middle)./photo_middle_central(i)^2).*photo_middle_real_central(i)^2;
xishu_plus=(photo_middle_30_central./attenuation_up_plus(i-1)./attenuation_down_effective_plus(i-1))./(photo_plus_30_central./attenuation_up_middle(i-1)./attenuation_down_effective_middle(i-1));
delta_photo_plus_square=xishu_plus^2.*(photo_plus_central(i)+photo_noise_plus)./photo_middle_central(i)^2.*photo_middle_real_central(i)^2;
xishu_minus=(photo_middle_30_central./attenuation_up_minus(i-1)./attenuation_down_effective_minus(i-1))./(photo_minus_30_central./attenuation_up_middle(i-1)./attenuation_down_effective_middle(i-1));
delta_photo_minus_square=xishu_minus^2.*(photo_minus_central(i)+photo_noise_minus)./photo_middle_central(i)^2.*photo_middle_real_central(i)^2;

% 计算温度误差的平方
delta_t_square(i)=((0.5.*dT_dRt+dT_dRv)^2.*delta_photo_plus_square+(0.5.*dT_dRt-dT_dRv)^2.*delta_photo_minus_square+  ...
                   (ratio_t_central(i).*dT_dRt+ratio_w_central(i).*dT_dRv)^2.*delta_photo_middle_square   ...
                  )./photo_middle_real_central(i)^2;
% 计算温度误差
delta_temperature(i)=sqrt(delta_t_square(i));

% 计算风速误差的平方
delta_w_square(i)=((0.5.*dV_dRt+dV_dRv)^2*delta_photo_plus_square+(0.5.*dV_dRt-dV_dRv)^2.*delta_photo_minus_square+  ...
                   (ratio_t_central(i).*dV_dRt+ratio_w_central(i).*dV_dRv)^2.*delta_photo_middle_square   ...
                   )./photo_middle_real_central(i)^2;
% 计算风速误差
delta_wind(i)=sqrt(delta_w_square(i));
    
    % 计算本层的钠原子有效后向散射截面-----------------------------------------------------
    frequency_x=frequency_0-3.5e9:1e5:frequency_0+3.5e9;  % [-3GH,+3.5GHz]区间内的频率
    wavelength_x=lightspeed./frequency_x;  % [-3GH,+3GHz]区间内的波长
    ke_sai=mass_sodium.*wavelength_x.^2./(2.*boltzmann_constant);  % 替换常量
    % 计算六种高斯曲线
gaussian_5_1=sqrt(ke_sai./(pi.*temperature_central(i))).*exp(-ke_sai.*(frequency_x-(wind_central(i)./wavelength_x)-frequency_5_1-frequency_0).^2./temperature_central(i));
gaussian_6_1=sqrt(ke_sai./(pi.*temperature_central(i))).*exp(-ke_sai.*(frequency_x-(wind_central(i)./wavelength_x)-frequency_6_1-frequency_0).^2./temperature_central(i));
gaussian_7_1=sqrt(ke_sai./(pi.*temperature_central(i))).*exp(-ke_sai.*(frequency_x-(wind_central(i)./wavelength_x)-frequency_7_1-frequency_0).^2./temperature_central(i));
gaussian_6_2=sqrt(ke_sai./(pi.*temperature_central(i))).*exp(-ke_sai.*(frequency_x-(wind_central(i)./wavelength_x)-frequency_6_2-frequency_0).^2./temperature_central(i));
gaussian_7_2=sqrt(ke_sai./(pi.*temperature_central(i))).*exp(-ke_sai.*(frequency_x-(wind_central(i)./wavelength_x)-frequency_7_2-frequency_0).^2./temperature_central(i));
gaussian_8_2=sqrt(ke_sai./(pi.*temperature_central(i))).*exp(-ke_sai.*(frequency_x-(wind_central(i)./wavelength_x)-frequency_8_2-frequency_0).^2./temperature_central(i));

A_0=6.15e7;
C=0.08502656;
N_1=3./(3+5.*exp(-C./temperature_central(i)));
N_2=5.*exp(-C./temperature_central(i))./(3+5.*exp(-C./temperature_central(i)));

delta_nu=A_0/(2*pi);
lorentzian=delta_nu./(2.*pi.*((frequency_x-frequency_0).^2+(delta_nu./2).^2));

% 将高斯线型和洛伦兹线型进行卷积，计算voigt线型
voigt_5_1=conv(gaussian_5_1,lorentzian,'same').*1e5;
voigt_6_1=conv(gaussian_6_1,lorentzian,'same').*1e5;
voigt_7_1=conv(gaussian_7_1,lorentzian,'same').*1e5;
voigt_6_2=conv(gaussian_6_2,lorentzian,'same').*1e5;
voigt_7_2=conv(gaussian_7_2,lorentzian,'same').*1e5;
voigt_8_2=conv(gaussian_8_2,lorentzian,'same').*1e5;
% 假设激光线型是高斯线型，半高全宽为80MHz，计算激光谱线
sigma_laser=80e6/(2*sqrt(2*log(2)));
lineshape_laser=(1./(sqrt(2.*pi).*sigma_laser)).*exp(-(frequency_x-frequency_0).^2./(2.*sigma_laser.^2));
% voigt线型与激光线型进行卷积
voigt_5_1_effective=conv(voigt_5_1,lineshape_laser,'same').*1e5;
voigt_6_1_effective=conv(voigt_6_1,lineshape_laser,'same').*1e5;
voigt_7_1_effective=conv(voigt_7_1,lineshape_laser,'same').*1e5;
voigt_6_2_effective=conv(voigt_6_2,lineshape_laser,'same').*1e5;
voigt_7_2_effective=conv(voigt_7_2,lineshape_laser,'same').*1e5;
voigt_8_2_effective=conv(voigt_8_2,lineshape_laser,'same').*1e5;
% 计算后向散射截面曲线(考虑激光线宽的情况)
sigma_back_effective=(A_0.*wavelength_x.^2./(32.*pi.^2)).*((N_1./6).*(2.*voigt_5_1_effective+5.5.*voigt_6_1_effective+5.*voigt_7_1_effective)+(N_2./10).*(0.98.*voigt_6_2_effective+5.*voigt_7_2_effective+15.68.*voigt_8_2_effective));
sigma_back_middle(i)=sigma_back_effective(frequency_x==x_middle);
% 计算吸收截面
sigma_absorption=(A_0.*wavelength_x.^2./(8.*pi)).*((N_1./6).*(2.*voigt_5_1_effective+5.*voigt_6_1_effective+5.*voigt_7_1_effective)+(N_2./10).*(voigt_6_2_effective+5.*voigt_7_2_effective+14.*voigt_8_2_effective));
sigma_absorption_middle(i)=sigma_absorption(frequency_x==x_middle);
% 计算大气的瑞利衰减
integral_atmosphere_replace=integral_atmosphere_replace+atmosphere_density(i).*extinction_cross_section_rayleigh*height_resolution;
integral_atmosphere(i)=integral_atmosphere_replace;
% 计算钠原子的密度
sodium_density(i)=((bin_origin.*dot_number.*i).^2./(bin_origin.*dot_number.*h_30).^2).*   ...
                  (backscattering_cross_section_rayleigh./sigma_back_middle(i)).*   ...
                  photo_middle_real_central(i).*exp(2.*integral_atmosphere(i)).*atmosphere_density(h_30);
% 计算向上的衰减
integral_up_replace_middle=integral_up_replace_middle+sodium_density(i)*sigma_absorption_middle(i)*height_resolution;
integral_up_middle(i)=integral_up_replace_middle;
attenuation_up_middle(i)=exp(-integral_up_middle(i));
% 计算向下的衰减
% 计算branching ratio of florescence emission
branching_ratio_0=((A_0.*wavelength_x.^2./(32.*pi.^2)).*(N_1./6).*alpha_0.*2.*voigt_5_1_effective)./sigma_back_effective;
branching_ratio_1=((A_0.*wavelength_x.^2./(32.*pi.^2)).*(N_1./6).*alpha_1.*5.*voigt_6_1_effective)./sigma_back_effective;
branching_ratio_2=((A_0.*wavelength_x.^2./(32.*pi.^2)).*(N_1./6).*alpha_2.*5.*voigt_7_1_effective)./sigma_back_effective;
branching_ratio_3=((A_0.*wavelength_x.^2./(32.*pi.^2)).*(N_1./6).*alpha_3.*5.*voigt_6_1_effective)./sigma_back_effective;
branching_ratio_4=((A_0.*wavelength_x.^2./(32.*pi.^2)).*(N_1./6).*alpha_4.*5.*voigt_7_1_effective)./sigma_back_effective;
branching_ratio_5=((A_0.*wavelength_x.^2./(32.*pi.^2)).*(N_2./10).*alpha_5.*voigt_6_2_effective)./sigma_back_effective;
branching_ratio_6=((A_0.*wavelength_x.^2./(32.*pi.^2)).*(N_2./10).*alpha_6.*5.*voigt_7_2_effective)./sigma_back_effective;
branching_ratio_7=((A_0.*wavelength_x.^2./(32.*pi.^2)).*(N_2./10).*alpha_7.*14.*voigt_8_2_effective)./sigma_back_effective;
branching_ratio_8=((A_0.*wavelength_x.^2./(32.*pi.^2)).*(N_2./10).*alpha_8.*voigt_6_2_effective)./sigma_back_effective;
branching_ratio_9=((A_0.*wavelength_x.^2./(32.*pi.^2)).*(N_2./10).*alpha_9.*5.*voigt_7_2_effective)./sigma_back_effective;
branching_ratio_0_middle_central(i)=branching_ratio_0(frequency_x==x_middle);
branching_ratio_1_middle_central(i)=branching_ratio_1(frequency_x==x_middle);
branching_ratio_2_middle_central(i)=branching_ratio_2(frequency_x==x_middle);
branching_ratio_3_middle_central(i)=branching_ratio_3(frequency_x==x_middle);
branching_ratio_4_middle_central(i)=branching_ratio_4(frequency_x==x_middle);
branching_ratio_5_middle_central(i)=branching_ratio_5(frequency_x==x_middle);
branching_ratio_6_middle_central(i)=branching_ratio_6(frequency_x==x_middle);
branching_ratio_7_middle_central(i)=branching_ratio_7(frequency_x==x_middle);
branching_ratio_8_middle_central(i)=branching_ratio_8(frequency_x==x_middle);
branching_ratio_9_middle_central(i)=branching_ratio_9(frequency_x==x_middle);
% 计算10个发射频率所对应的吸收截面
sigma_absorption_middle_0(i)=sigma_absorption(frequency_x==frequency_0+emission_0_middle);
sigma_absorption_middle_1(i)=sigma_absorption(frequency_x==frequency_0+emission_1_middle);
sigma_absorption_middle_2(i)=sigma_absorption(frequency_x==frequency_0+emission_2_middle);
sigma_absorption_middle_3(i)=sigma_absorption(frequency_x==frequency_0+emission_3_middle);
sigma_absorption_middle_4(i)=sigma_absorption(frequency_x==frequency_0+emission_4_middle);
sigma_absorption_middle_5(i)=sigma_absorption(frequency_x==frequency_0+emission_5_middle);
sigma_absorption_middle_6(i)=sigma_absorption(frequency_x==frequency_0+emission_6_middle);
sigma_absorption_middle_7(i)=sigma_absorption(frequency_x==frequency_0+emission_7_middle);
sigma_absorption_middle_8(i)=sigma_absorption(frequency_x==frequency_0+emission_8_middle);
sigma_absorption_middle_9(i)=sigma_absorption(frequency_x==frequency_0+emission_9_middle);

integral_down_replace_middle_0=integral_down_replace_middle_0+sodium_density(i)*sigma_absorption_middle_0(i)*height_resolution;
integral_down_middle_0(i)=integral_down_replace_middle_0;
attenuation_down_middle_0(i)=exp(-integral_down_middle_0(i));
integral_down_replace_middle_1=integral_down_replace_middle_1+sodium_density(i)*sigma_absorption_middle_1(i)*height_resolution;
integral_down_middle_1(i)=integral_down_replace_middle_1;
attenuation_down_middle_1(i)=exp(-integral_down_middle_1(i));
integral_down_replace_middle_2=integral_down_replace_middle_2+sodium_density(i)*sigma_absorption_middle_2(i)*height_resolution;
integral_down_middle_2(i)=integral_down_replace_middle_2;
attenuation_down_middle_2(i)=exp(-integral_down_middle_2(i));
integral_down_replace_middle_3=integral_down_replace_middle_3+sodium_density(i)*sigma_absorption_middle_3(i)*height_resolution;
integral_down_middle_3(i)=integral_down_replace_middle_3;
attenuation_down_middle_3(i)=exp(-integral_down_middle_3(i));
integral_down_replace_middle_4=integral_down_replace_middle_4+sodium_density(i)*sigma_absorption_middle_4(i)*height_resolution;
integral_down_middle_4(i)=integral_down_replace_middle_4;
attenuation_down_middle_4(i)=exp(-integral_down_middle_4(i));
integral_down_replace_middle_5=integral_down_replace_middle_5+sodium_density(i)*sigma_absorption_middle_5(i)*height_resolution;
integral_down_middle_5(i)=integral_down_replace_middle_5;
attenuation_down_middle_5(i)=exp(-integral_down_middle_5(i));
integral_down_replace_middle_6=integral_down_replace_middle_6+sodium_density(i)*sigma_absorption_middle_6(i)*height_resolution;
integral_down_middle_6(i)=integral_down_replace_middle_6;
attenuation_down_middle_6(i)=exp(-integral_down_middle_6(i));
integral_down_replace_middle_7=integral_down_replace_middle_7+sodium_density(i)*sigma_absorption_middle_7(i)*height_resolution;
integral_down_middle_7(i)=integral_down_replace_middle_7;
attenuation_down_middle_7(i)=exp(-integral_down_middle_7(i));
integral_down_replace_middle_8=integral_down_replace_middle_8+sodium_density(i)*sigma_absorption_middle_8(i)*height_resolution;
integral_down_middle_8(i)=integral_down_replace_middle_8;
attenuation_down_middle_8(i)=exp(-integral_down_middle_8(i));
integral_down_replace_middle_9=integral_down_replace_middle_9+sodium_density(i)*sigma_absorption_middle_9(i)*height_resolution;
integral_down_middle_9(i)=integral_down_replace_middle_9;
attenuation_down_middle_9(i)=exp(-integral_down_middle_9(i));
% 计算向下有效衰减
attenuation_down_effective_middle(i)=branching_ratio_0_middle_central(i).*attenuation_down_middle_0(i)+   ...
                                     branching_ratio_1_middle_central(i).*attenuation_down_middle_1(i)+   ...
                                     branching_ratio_2_middle_central(i).*attenuation_down_middle_2(i)+   ...
                                     branching_ratio_3_middle_central(i).*attenuation_down_middle_3(i)+   ...
                                     branching_ratio_4_middle_central(i).*attenuation_down_middle_4(i)+   ...
                                     branching_ratio_5_middle_central(i).*attenuation_down_middle_5(i)+   ...
                                     branching_ratio_6_middle_central(i).*attenuation_down_middle_6(i)+   ...
                                     branching_ratio_7_middle_central(i).*attenuation_down_middle_7(i)+   ...
                                     branching_ratio_8_middle_central(i).*attenuation_down_middle_8(i)+   ...
                                     branching_ratio_9_middle_central(i).*attenuation_down_middle_9(i);

% 计算plus的衰减------------------------------------------------------------
sigma_absorption_plus(i)=sigma_absorption(frequency_x==x_plus);
% 计算plus的向上衰减
integral_up_replace_plus=integral_up_replace_plus+sodium_density(i)*sigma_absorption_plus(i)*height_resolution;
integral_up_plus(i)=integral_up_replace_plus;
attenuation_up_plus(i)=exp(-integral_up_plus(i));
% 计算plus的向下衰减
% 计算分支比
branching_ratio_0_plus_central(i)=branching_ratio_0(frequency_x==x_plus);
branching_ratio_1_plus_central(i)=branching_ratio_1(frequency_x==x_plus);
branching_ratio_2_plus_central(i)=branching_ratio_2(frequency_x==x_plus);
branching_ratio_3_plus_central(i)=branching_ratio_3(frequency_x==x_plus);
branching_ratio_4_plus_central(i)=branching_ratio_4(frequency_x==x_plus);
branching_ratio_5_plus_central(i)=branching_ratio_5(frequency_x==x_plus);
branching_ratio_6_plus_central(i)=branching_ratio_6(frequency_x==x_plus);
branching_ratio_7_plus_central(i)=branching_ratio_7(frequency_x==x_plus);
branching_ratio_8_plus_central(i)=branching_ratio_8(frequency_x==x_plus);
branching_ratio_9_plus_central(i)=branching_ratio_9(frequency_x==x_plus);
% 计算10个发射频率所对应的吸收截面，发射频率取决于激光的入射频率和能级状态
sigma_absorption_plus_0(i)=sigma_absorption(frequency_x==frequency_0+emission_0_plus);
sigma_absorption_plus_1(i)=sigma_absorption(frequency_x==frequency_0+emission_1_plus);
sigma_absorption_plus_2(i)=sigma_absorption(frequency_x==frequency_0+emission_2_plus);
sigma_absorption_plus_3(i)=sigma_absorption(frequency_x==frequency_0+emission_3_plus);
sigma_absorption_plus_4(i)=sigma_absorption(frequency_x==frequency_0+emission_4_plus);
sigma_absorption_plus_5(i)=sigma_absorption(frequency_x==frequency_0+emission_5_plus);
sigma_absorption_plus_6(i)=sigma_absorption(frequency_x==frequency_0+emission_6_plus);
sigma_absorption_plus_7(i)=sigma_absorption(frequency_x==frequency_0+emission_7_plus);
sigma_absorption_plus_8(i)=sigma_absorption(frequency_x==frequency_0+emission_8_plus);
sigma_absorption_plus_9(i)=sigma_absorption(frequency_x==frequency_0+emission_9_plus);

integral_down_replace_plus_0=integral_down_replace_plus_0+sodium_density(i)*sigma_absorption_plus_0(i)*height_resolution;
integral_down_plus_0(i)=integral_down_replace_plus_0;
attenuation_down_plus_0(i)=exp(-integral_down_plus_0(i));
integral_down_replace_plus_1=integral_down_replace_plus_1+sodium_density(i)*sigma_absorption_plus_1(i)*height_resolution;
integral_down_plus_1(i)=integral_down_replace_plus_1;
attenuation_down_plus_1(i)=exp(-integral_down_plus_1(i));
integral_down_replace_plus_2=integral_down_replace_plus_2+sodium_density(i)*sigma_absorption_plus_2(i)*height_resolution;
integral_down_plus_2(i)=integral_down_replace_plus_2;
attenuation_down_plus_2(i)=exp(-integral_down_plus_2(i));
integral_down_replace_plus_3=integral_down_replace_plus_3+sodium_density(i)*sigma_absorption_plus_3(i)*height_resolution;
integral_down_plus_3(i)=integral_down_replace_plus_3;
attenuation_down_plus_3(i)=exp(-integral_down_plus_3(i));
integral_down_replace_plus_4=integral_down_replace_plus_4+sodium_density(i)*sigma_absorption_plus_4(i)*height_resolution;
integral_down_plus_4(i)=integral_down_replace_plus_4;
attenuation_down_plus_4(i)=exp(-integral_down_plus_4(i));
integral_down_replace_plus_5=integral_down_replace_plus_5+sodium_density(i)*sigma_absorption_plus_5(i)*height_resolution;
integral_down_plus_5(i)=integral_down_replace_plus_5;
attenuation_down_plus_5(i)=exp(-integral_down_plus_5(i));
integral_down_replace_plus_6=integral_down_replace_plus_6+sodium_density(i)*sigma_absorption_plus_6(i)*height_resolution;
integral_down_plus_6(i)=integral_down_replace_plus_6;
attenuation_down_plus_6(i)=exp(-integral_down_plus_6(i));
integral_down_replace_plus_7=integral_down_replace_plus_7+sodium_density(i)*sigma_absorption_plus_7(i)*height_resolution;
integral_down_plus_7(i)=integral_down_replace_plus_7;
attenuation_down_plus_7(i)=exp(-integral_down_plus_7(i));
integral_down_replace_plus_8=integral_down_replace_plus_8+sodium_density(i)*sigma_absorption_plus_8(i)*height_resolution;
integral_down_plus_8(i)=integral_down_replace_plus_8;
attenuation_down_plus_8(i)=exp(-integral_down_plus_8(i));
integral_down_replace_plus_9=integral_down_replace_plus_9+sodium_density(i)*sigma_absorption_plus_9(i)*height_resolution;
integral_down_plus_9(i)=integral_down_replace_plus_9;
attenuation_down_plus_9(i)=exp(-integral_down_plus_9(i));
% 计算plus的有效向下衰减
attenuation_down_effective_plus(i)=branching_ratio_0_plus_central(i).*attenuation_down_plus_0(i)+   ...
                                     branching_ratio_1_plus_central(i).*attenuation_down_plus_1(i)+   ...
                                     branching_ratio_2_plus_central(i).*attenuation_down_plus_2(i)+   ...
                                     branching_ratio_3_plus_central(i).*attenuation_down_plus_3(i)+   ...
                                     branching_ratio_4_plus_central(i).*attenuation_down_plus_4(i)+   ...
                                     branching_ratio_5_plus_central(i).*attenuation_down_plus_5(i)+   ...
                                     branching_ratio_6_plus_central(i).*attenuation_down_plus_6(i)+   ...
                                     branching_ratio_7_plus_central(i).*attenuation_down_plus_7(i)+   ...
                                     branching_ratio_8_plus_central(i).*attenuation_down_plus_8(i)+   ...
                                     branching_ratio_9_plus_central(i).*attenuation_down_plus_9(i);


% 计算minus的衰减------------------------------------------------------------
sigma_absorption_minus(i)=sigma_absorption(frequency_x==x_minus);
% 计算minus的向上衰减
integral_up_replace_minus=integral_up_replace_minus+sodium_density(i)*sigma_absorption_minus(i)*height_resolution;
integral_up_minus(i)=integral_up_replace_minus;
attenuation_up_minus(i)=exp(-integral_up_minus(i));
% 计算minus的向下衰减
% 计算分支比
branching_ratio_0_minus_central(i)=branching_ratio_0(frequency_x==x_minus);
branching_ratio_1_minus_central(i)=branching_ratio_1(frequency_x==x_minus);
branching_ratio_2_minus_central(i)=branching_ratio_2(frequency_x==x_minus);
branching_ratio_3_minus_central(i)=branching_ratio_3(frequency_x==x_minus);
branching_ratio_4_minus_central(i)=branching_ratio_4(frequency_x==x_minus);
branching_ratio_5_minus_central(i)=branching_ratio_5(frequency_x==x_minus);
branching_ratio_6_minus_central(i)=branching_ratio_6(frequency_x==x_minus);
branching_ratio_7_minus_central(i)=branching_ratio_7(frequency_x==x_minus);
branching_ratio_8_minus_central(i)=branching_ratio_8(frequency_x==x_minus);
branching_ratio_9_minus_central(i)=branching_ratio_9(frequency_x==x_minus);
% 计算10个发射频率所对应的吸收截面，发射频率取决于激光的入射频率和能级状态
sigma_absorption_minus_0(i)=sigma_absorption(frequency_x==frequency_0+emission_0_minus);
sigma_absorption_minus_1(i)=sigma_absorption(frequency_x==frequency_0+emission_1_minus);
sigma_absorption_minus_2(i)=sigma_absorption(frequency_x==frequency_0+emission_2_minus);
sigma_absorption_minus_3(i)=sigma_absorption(frequency_x==frequency_0+emission_3_minus);
sigma_absorption_minus_4(i)=sigma_absorption(frequency_x==frequency_0+emission_4_minus);
sigma_absorption_minus_5(i)=sigma_absorption(frequency_x==frequency_0+emission_5_minus);
sigma_absorption_minus_6(i)=sigma_absorption(frequency_x==frequency_0+emission_6_minus);
sigma_absorption_minus_7(i)=sigma_absorption(frequency_x==frequency_0+emission_7_minus);
sigma_absorption_minus_8(i)=sigma_absorption(frequency_x==frequency_0+emission_8_minus);
sigma_absorption_minus_9(i)=sigma_absorption(frequency_x==frequency_0+emission_9_minus);

integral_down_replace_minus_0=integral_down_replace_minus_0+sodium_density(i)*sigma_absorption_minus_0(i)*height_resolution;
integral_down_minus_0(i)=integral_down_replace_minus_0;
attenuation_down_minus_0(i)=exp(-integral_down_minus_0(i));
integral_down_replace_minus_1=integral_down_replace_minus_1+sodium_density(i)*sigma_absorption_minus_1(i)*height_resolution;
integral_down_minus_1(i)=integral_down_replace_minus_1;
attenuation_down_minus_1(i)=exp(-integral_down_minus_1(i));
integral_down_replace_minus_2=integral_down_replace_minus_2+sodium_density(i)*sigma_absorption_minus_2(i)*height_resolution;
integral_down_minus_2(i)=integral_down_replace_minus_2;
attenuation_down_minus_2(i)=exp(-integral_down_minus_2(i));
integral_down_replace_minus_3=integral_down_replace_minus_3+sodium_density(i)*sigma_absorption_minus_3(i)*height_resolution;
integral_down_minus_3(i)=integral_down_replace_minus_3;
attenuation_down_minus_3(i)=exp(-integral_down_minus_3(i));
integral_down_replace_minus_4=integral_down_replace_minus_4+sodium_density(i)*sigma_absorption_minus_4(i)*height_resolution;
integral_down_minus_4(i)=integral_down_replace_minus_4;
attenuation_down_minus_4(i)=exp(-integral_down_minus_4(i));
integral_down_replace_minus_5=integral_down_replace_minus_5+sodium_density(i)*sigma_absorption_minus_5(i)*height_resolution;
integral_down_minus_5(i)=integral_down_replace_minus_5;
attenuation_down_minus_5(i)=exp(-integral_down_minus_5(i));
integral_down_replace_minus_6=integral_down_replace_minus_6+sodium_density(i)*sigma_absorption_minus_6(i)*height_resolution;
integral_down_minus_6(i)=integral_down_replace_minus_6;
attenuation_down_minus_6(i)=exp(-integral_down_minus_6(i));
integral_down_replace_minus_7=integral_down_replace_minus_7+sodium_density(i)*sigma_absorption_minus_7(i)*height_resolution;
integral_down_minus_7(i)=integral_down_replace_minus_7;
attenuation_down_minus_7(i)=exp(-integral_down_minus_7(i));
integral_down_replace_minus_8=integral_down_replace_minus_8+sodium_density(i)*sigma_absorption_minus_8(i)*height_resolution;
integral_down_minus_8(i)=integral_down_replace_minus_8;
attenuation_down_minus_8(i)=exp(-integral_down_minus_8(i));
integral_down_replace_minus_9=integral_down_replace_minus_9+sodium_density(i)*sigma_absorption_minus_9(i)*height_resolution;
integral_down_minus_9(i)=integral_down_replace_minus_9;
attenuation_down_minus_9(i)=exp(-integral_down_minus_9(i));
% 计算minus的有效向下衰减
attenuation_down_effective_minus(i)=branching_ratio_0_minus_central(i).*attenuation_down_minus_0(i)+   ...
                                     branching_ratio_1_minus_central(i).*attenuation_down_minus_1(i)+   ...
                                     branching_ratio_2_minus_central(i).*attenuation_down_minus_2(i)+   ...
                                     branching_ratio_3_minus_central(i).*attenuation_down_minus_3(i)+   ...
                                     branching_ratio_4_minus_central(i).*attenuation_down_minus_4(i)+   ...
                                     branching_ratio_5_minus_central(i).*attenuation_down_minus_5(i)+   ...
                                     branching_ratio_6_minus_central(i).*attenuation_down_minus_6(i)+   ...
                                     branching_ratio_7_minus_central(i).*attenuation_down_minus_7(i)+   ...
                                     branching_ratio_8_minus_central(i).*attenuation_down_minus_8(i)+   ...
                                     branching_ratio_9_minus_central(i).*attenuation_down_minus_9(i);
                                 
                                 
end
% x=(90:116).*0.9;
% plot(x,photo_minus_real_central(90:116),'r',x,photo_minus_cali_central(90:116),'b');

% 写出温度与风速文件
filename2=strcat('E:\temperature_wind_calculate\picture_final\',pulse_number,'_central_temperature_wind.txt');
fid=fopen(filename2,'w');
for i=h_bottom:h_top;
sodium_density_cm3(i)=sodium_density(i)/1e6;
number=i-h_bottom+1;
fprintf(fid,'%10.1f %10.1f %10.3f %10.1f %10.3f %10.1f %10d\r\n',height_merge(i),temperature_central(i),delta_temperature(i),wind_central(i),delta_wind(i),sodium_density_cm3(i),number);
end

% plot(frequency_x,sigma_absorption,'r',frequency_x,sigma_back_effective,'b');


% % 画出垂直温度
% figure(3);
% h_bottom_central=96;
% h_top_central=113;
% plot(temperature_central(h_bottom_central:h_top_central),height_merge(h_bottom_central:h_top_central));
% picture_name=strcat('E:\temperature_wind_calculate\picture_result\',pulse_number,'_central_temperature');
% saveas(gcf,picture_name,'tif');
% 
% % 画出垂直风场
% figure(4);
% plot(wind_central(h_bottom_central:h_top_central),height_merge(h_bottom_central:h_top_central));
% picture_name=strcat('E:\temperature_wind_calculate\picture_result\',pulse_number,'_central_wind');
% saveas(gcf,picture_name,'tif');
% 




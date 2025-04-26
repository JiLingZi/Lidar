%% 读Excel
clc;

XlsName ='D:\KYBF\TrailPPT\Datas\粗高度时间分布.xlsx';
Table = readtable(XlsName,'sheet','CapAlt');
AltCeil = table2array(Table(:,4));
AltCeil = AltCeil(~isnan(AltCeil));
DenSum = table2array(Table(:,5));
DenSum = DenSum(~isnan(DenSum));
AlCnt = table2array(Table(:,6));
AlCnt = AlCnt(~isnan(AlCnt));

%% 各粒子密度-Na
DateList = {'20231022';'20231103';'20240114';'20240117';'20240118';...
            '20240119';'20240120';'20240128';'20240129';'20240519';...
            '20240529';'20240530'};

for jd = 1:size(DateList,1)
% 读取光子数矩阵
ReadDate = char(DateList(jd));
RunDate = ReadDate;

Year = str2double(char(RunDate(1:4)));
Month = str2double(char(RunDate(5:6)));
Day = str2double(char(RunDate(7:8)));
DayNum = Month*30-30+Day;

NajNum = 15; NaiNum = 10;
KjNum = 15; KiNum = 10;
FejNum = 15; FeiNum = 10;
NijNum = 1; NiiNum = 1;
CajNum = 10; CaiNum = 10;
CapjNum = 15; CapiNum = 10;
RayjNum = 1; RayiNum = 1;

% fldstart = ['F:\RawData\ZWDATA\MOHEnew\Na\',ReadDate,'\'];
% folder = [fldstart,'Na\'];
% files = dir(fullfile(folder, '*.dat'));
% num_files = length(files)-1;% 为了统一文件数，减一
% 
% V_F_0 = zeros(8192, num_files);
% V_F_R = zeros(8192, num_files);
% V_F_L = zeros(8192, num_files);
% N_F_0 = zeros(8192, num_files);
% N_F_R = zeros(8192, num_files);
% N_F_L = zeros(8192, num_files);
% E_F_0 = zeros(8192, num_files);
% E_F_R = zeros(8192, num_files);
% E_F_L = zeros(8192, num_files);
% for j = 1:num_files
%     filename = fullfile(folder, files(j).name);
%     Na_table = readtable(filename);
%     Na_data = table2array(Na_table);
%     V_F_0(:,j) = Na_data(1:8192,2);
%     V_F_R(:,j) = Na_data(1:8192,3);
%     V_F_L(:,j) = Na_data(1:8192,4);
%     N_F_0(:,j) = Na_data(1:8192,5);
%     N_F_R(:,j) = Na_data(1:8192,6);
%     N_F_L(:,j) = Na_data(1:8192,7);
%     E_F_0(:,j) = Na_data(1:8192,8);
%     E_F_R(:,j) = Na_data(1:8192,9);
%     E_F_L(:,j) = Na_data(1:8192,10);
% end
% 
% % 读取时间序列
% TimeValues = cell(1,num_files);
% for j = 1:num_files
%     filename = fullfile(folder, files(j).name);
%     fid = fopen(filename);
%     TXT = textscan(fid, '%s', 'Delimiter', '\n');
%     RowTime = TXT{1}{4};
%     TimeData = strsplit(RowTime);
%     fclose(fid);
%     TimeValues(:,j) = cellstr(TimeData{4});
% end

    fldstart = ['F:\RawData\ZWDATA\MOHEnew\Cap\',ReadDate,'\'];
    folder = [fldstart,'CH1\'];
    files = dir(fullfile(folder, '*.txt'));
    num_files = length(files)-1;

    FileCap = num_files;

    DenPh1 = zeros(4096, num_files);
    folder = [fldstart,'CH1\'];
    files = dir(fullfile(folder, '*.txt'));
    for j = 1:num_files
        filename = fullfile(folder,files(j).name);
        fid = fopen(filename);
        PhData1 = textscan(fid, '%f %f','HeaderLines',20);
        fclose(fid);
        DenPh1(:,j) = PhData1{1,2}(1:4096,:);
    end

    DenPh2 = zeros(4096, num_files);
    folder = [fldstart,'CH2\'];
    files = dir(fullfile(folder, '*.txt'));
    for j = 1:num_files
        filename = fullfile(folder,files(j).name);
        fid = fopen(filename);
        PhData2 = textscan(fid, '%f %f','HeaderLines',20);
        fclose(fid);
        DenPh2(:,j) = PhData2{1,2}(1:4096,:);
    end

    DenPh3 = zeros(4096, num_files);
    folder = [fldstart,'CH3\'];
    files = dir(fullfile(folder, '*.txt'));
    for j = 1:num_files
        filename = fullfile(folder,files(j).name);
        fid = fopen(filename);
        PhData3 = textscan(fid, '%f %f','HeaderLines',20);
        fclose(fid);
        DenPh3(:,j) = PhData3{1,2}(1:4096,:);
    end

    DenPh4 = zeros(4096, num_files);
    folder = [fldstart,'CH4\'];
    files = dir(fullfile(folder, '*.txt'));
    for j = 1:num_files
        filename = fullfile(folder,files(j).name);
        fid = fopen(filename);
        PhData4 = textscan(fid, '%f %f','HeaderLines',20);
        fclose(fid);
        DenPh4(:,j) = PhData4{1,2}(1:4096,:);
    end

    DenPh = DenPh1+DenPh2+DenPh3+DenPh4;
    
    TimeValues = cell(1,num_files);
    for j = 1:num_files
        filename = fullfile(folder, files(j).name);
        fid = fopen(filename);
        TXT = textscan(fid, '%s', 'Delimiter', '\n');
        RowTime = TXT{1}{9};
        TimeData = strsplit(RowTime);
        fclose(fid);
        TimeValues(:,j) = cellstr(TimeData{2:end});
    end
    jNum = NajNum;
    iNum = NaiNum;
    TimeList = cell(1,floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
    end

%     TimeXCap = datetime(TimeList, 'InputFormat', ':''yyyy-MM-dd''T''HH:mm:ss.SSS');
    height_num_origin = PhData1{1,1}(1:4096);
%     height_num_origin = Na_data(1:8192,1);


    PhiSum = zeros(floor(4096/iNum),num_files);
    for i = 1:floor(4096/iNum)
        PhiSum(i,:) = sum(DenPh(iNum*(i-1)+1:iNum*i,:),1);
    end
    PhjSumCap = zeros(floor(4096/iNum),floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        PhjSumCap(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
    end
    AltitudeCap = height_num_origin(1:floor(4096/iNum))*iNum;
    KM30Cap = size(AltitudeCap(AltitudeCap<30),1)+1;
    KM75Cap = size(AltitudeCap(AltitudeCap<75),1)+1;
    KM80Cap = size(AltitudeCap(AltitudeCap<80),1)+1;
    KM85Cap = size(AltitudeCap(AltitudeCap<85),1)+1;
    KM90Cap = size(AltitudeCap(AltitudeCap<90),1)+1;
    KM95Cap = size(AltitudeCap(AltitudeCap<95),1)+1;
    KM100Cap = size(AltitudeCap(AltitudeCap<100),1)+1;
    KM105Cap = size(AltitudeCap(AltitudeCap<105),1)+1;
    KM110Cap = size(AltitudeCap(AltitudeCap<110),1)+1;
    KM115Cap = size(AltitudeCap(AltitudeCap<115),1)+1;
    KM120Cap = size(AltitudeCap(AltitudeCap<120),1)+1;
    KM125Cap = size(AltitudeCap(AltitudeCap<125),1)+1;
    
%     TimeShuCap = TimeXCap';

    Noise = mean(PhjSumCap(KM120Cap:KM125Cap,:),1);

    ScatterTable = readtable('RayEffScatter.txt');
    Scatter = table2array(ScatterTable(:,2:3));
    RayScatter = Scatter(3,1);
    EffScatter = Scatter(3,2);

    [TemRay,DenRay] = atmosnrlmsise00(30000,53.3,122.7,Year,DayNum,17);

    Z = AltitudeCap;                       % 高度矩阵
    ZR = 30;                            % 参考高度
    SigmaRay = RayScatter;              % 瑞利后向散射截面
    SigmaEff = EffScatter;              % 有效后向散射截面
    NZ = PhjSumCap;                        % 光子数矩阵
    NB = Noise;                         % 噪声矩阵
    NZR = PhjSumCap(KM30Cap,:);                % 参考高度处光子数
    NumberRay = sum(DenRay)*1e-6;       % 参考高度处大气模型数密度
    NumberZ = (((Z.^2).*SigmaRay.*(NZ-NB))./((ZR^2)*SigmaEff*(NZR-NB)))*NumberRay;

    DayMean = mean(NumberZ,2);
    DenMeanSum(:,jd) = DayMean;
    TS = '目前第'+string(jd)
end
DenMean = mean(DenMeanSum,2);

%% 画图
figure('name','密度图')
yyaxis left
plot(AltCeil,DenSum,'-b','linewidth',1.5);
ylabel('Density (cm^{-3})');
xlabel('Altitude (km)');
title('Density Distribution of Na Meteor Trails');
% xlim([0 20]);
xlim([75 115]);
grid on;

yyaxis right
plot(AltitudeCap,DenMean,'-r','linewidth',1.5);
ylabel('Density (cm^{-3})');
xlabel('Altitude (km)');
title('Mean Density of Na');
% ylim([0 100]);
xlim([75 115]);
grid on;

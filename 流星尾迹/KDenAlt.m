%% 钠钾铁镍尾迹密度与海拔读取

clc; clear all;

% 通用list

DateList = {'20231018';'20231019';'20231020';'20231022';'20231023';...
            '20231026';'20231027';'20231028';'20231102';'20231103';...
            '20231104';'20231105';'20231106';'20231107';'20231111';...
            '20231112';'20231113';'20231118';'20231119';'20231121';...
            '20231122';'20231123';...
            '20240114';'20240115';'20240117';'20240118';'20240119';...
            '20240125';'20240128';'20240131';...
            '20240229';'20240301';'20240304';'20240305';...
            '20240412';'20240415';'20240418';...
            '20240425';'20240428';...
            '20240512';...
            '20240525';...
            '20240528';'20240530';'20240531';};

%% K 读取

DateList = {'20240510';...
            '20240526';};

for jd = 1:size(DateList,1)
% 读取光子数矩阵
ReadDate = char(DateList(jd));
RunDate = ReadDate;

Year = str2double(char(RunDate(1:4)));
Month = str2double(char(RunDate(5:6)));
Day = str2double(char(RunDate(7:8)));
DayNum = Month*30-30+Day;

CHN = 'K';
XlsName = char("F:\TrailIMG\" + CHN + ".xlsx");
SDate = RunDate(5:8);
SheetName = char('S'+string(SDate));
Table = readtable(XlsName,'sheet',SheetName);
NaVList = table2array(Table(:,4));
NaVList = NaVList(~isnan(NaVList));
NaMList = table2array(Table(:,5));
NaMList = NaMList(~isnan(NaMList));
NaZList = table2array(Table(:,6));
NaZList = NaZList(~isnan(NaZList));
KList = table2array(Table(:,7));
KList = KList(~isnan(KList));
FeList = table2array(Table(:,8));
FeList = FeList(~isnan(FeList));
NiList = table2array(Table(:,9));
NiList = NiList(~isnan(NiList));
CaList = table2array(Table(:,10));
CaList = CaList(~isnan(CaList));
CapList = table2array(Table(:,11));
CapList = CapList(~isnan(CapList));
RayList = table2array(Table(:,12));
RayList = RayList(~isnan(RayList));

NajNum = 1; NaiNum = 1;
KjNum = 1; KiNum = 1;
FejNum = 1; FeiNum = 1;
NijNum = 1; NiiNum = 1;
CajNum = 1; CaiNum = 1;
CapjNum = 1; CapiNum = 1;
RayjNum = 1; RayiNum = 1;

KListNum = size(KList,1);

NumLins = KListNum;

if NumLins>0

    fldstart = ['F:\RawData\ZWDATA\MOHEnew\K532\',ReadDate,'\'];
    folder = [fldstart,'CH1\'];
    files = dir(fullfile(folder, '*.txt'));
    num_files = length(files)-1;

    FileK = num_files;

    DenPh = zeros(4096, num_files);
    for j = 1:num_files
        filename = fullfile(folder, files(j).name);
        fid = fopen(filename);
        PhData = textscan(fid, '%f %f','HeaderLines',20);   % 7月22日及之前均为23，后应改为20
        fclose(fid);
        DenPh(:,j) = PhData{1,2}(1:4096);
    end
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
    jNum = KjNum;
    iNum = KiNum;
    TimeList = cell(1,floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
    end

    TimeXK = datetime(TimeList, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');
    height_num_origin = PhData{1,1}(1:4096);

    PhiSum = zeros(floor(4096/iNum),num_files);
    for i = 1:floor(4096/iNum)
        PhiSum(i,:) = sum(DenPh(iNum*(i-1)+1:iNum*i,:),1);
    end
    PhjSumK = zeros(floor(4096/iNum),floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        PhjSumK(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
    end
    AltitudeK = height_num_origin(1:floor(4096/iNum))*iNum;
    KM30K = size(AltitudeK(AltitudeK<30),1)+1;
    KM75K = size(AltitudeK(AltitudeK<75),1)+1;
    KM80K = size(AltitudeK(AltitudeK<80),1)+1;
    KM85K = size(AltitudeK(AltitudeK<85),1)+1;
    KM90K = size(AltitudeK(AltitudeK<90),1)+1;
    KM95K = size(AltitudeK(AltitudeK<95),1)+1;
    KM100K = size(AltitudeK(AltitudeK<100),1)+1;
    KM105K = size(AltitudeK(AltitudeK<105),1)+1;
    KM110K = size(AltitudeK(AltitudeK<110),1)+1;
    KM115K = size(AltitudeK(AltitudeK<115),1)+1;
    KM120K = size(AltitudeK(AltitudeK<120),1)+1;
    KM125K = size(AltitudeK(AltitudeK<125),1)+1;
    
    TimeShuK = TimeXK';

    Noise = mean(PhjSumK(KM120K:KM125K,:),1);

    ScatterTable = readtable('RayEffScatter.txt');
    Scatter = table2array(ScatterTable(:,2:3));
    RayScatter = Scatter(6,1);
    EffScatter = Scatter(6,2);

    [TemRay,DenRay] = atmosnrlmsise00(30000,53.3,122.7,Year,DayNum,17);

    Z = AltitudeK;                       % 高度矩阵
    ZR = 30;                            % 参考高度
    SigmaRay = RayScatter;              % 瑞利后向散射截面
    SigmaEff = EffScatter;              % 有效后向散射截面
    NZ = PhjSumK;                        % 光子数矩阵
    NB = Noise;                         % 噪声矩阵
    NZR = PhjSumK(KM30K,:);                % 参考高度处光子数
    NumberRay = sum(DenRay)*1e-6;       % 参考高度处大气模型数密度
    NumberZ = (((Z.^2).*SigmaRay.*(NZ-NB))./((ZR^2)*SigmaEff*(NZR-NB)))*NumberRay;

    Numjs = size(KList,1);
    
    for jpf = 1:Numjs
        Prof = KList(jpf);

        fK = figure('name','K');
        
        [Max75,Aax75] = max(NumberZ(KM75K:KM80K,Prof));
        Aax75 = round(AltitudeK(Aax75+KM75K),2);
        AVG75 = (sum(NumberZ(KM75K:KM80K,Prof))-Max75)./size(NumberZ(KM75K:KM80K,Prof),1);
        Max75 = round(Max75-AVG75);
        [Max80,Aax80] = max(NumberZ(KM80K:KM85K,Prof));
        Aax80 = round(AltitudeK(Aax80+KM80K),2);
        AVG80 = (sum(NumberZ(KM80K:KM85K,Prof))-Max80)./size(NumberZ(KM80K:KM85K,Prof),1);
        Max80 = round(Max80-AVG80);
        [Max85,Aax85] = max(NumberZ(KM85K:KM90K,Prof));
        Aax85 = round(AltitudeK(Aax85+KM85K),2);
        AVG85 = (sum(NumberZ(KM85K:KM90K,Prof))-Max85)./size(NumberZ(KM85K:KM90K,Prof),1);
        Max85 = round(Max85-AVG85);
        [Max90,Aax90] = max(NumberZ(KM90K:KM95K,Prof));
        Aax90 = round(AltitudeK(Aax90+KM90K),2);
        AVG90 = (sum(NumberZ(KM90K:KM95K,Prof))-Max90)./size(NumberZ(KM90K:KM95K,Prof),1);
        Max90 = round(Max90-AVG90);
        [Max95,Aax95] = max(NumberZ(KM95K:KM100K,Prof));
        Aax95 = round(AltitudeK(Aax95+KM95K),2);
        AVG95 = (sum(NumberZ(KM95K:KM100K,Prof))-Max95)./size(NumberZ(KM95K:KM100K,Prof),1);
        Max95 = round(Max95-AVG95);
        [Max100,Aax100] = max(NumberZ(KM100K:KM105K,Prof));
        Aax100 = round(AltitudeK(Aax100+KM100K),2);
        AVG100 = (sum(NumberZ(KM100K:KM105K,Prof))-Max100)./size(NumberZ(KM100K:KM105K,Prof),1);
        Max100 = round(Max100-AVG100);
        [Max105,Aax105] = max(NumberZ(KM105K:KM110K,Prof));
        Aax105 = round(AltitudeK(Aax105+KM105K),2);
        AVG105 = (sum(NumberZ(KM105K:KM110K,Prof))-Max105)./size(NumberZ(KM105K:KM110K,Prof),1);
        Max105 = round(Max105-AVG105);
        [Max110,Aax110] = max(NumberZ(KM110K:KM115K,Prof));
        Aax110 = round(AltitudeK(Aax110+KM110K),2);
        AVG110 = (sum(NumberZ(KM110K:KM115K,Prof))-Max110)./size(NumberZ(KM110K:KM115K,Prof),1);
        Max110 = round(Max110-AVG110);
        
        Ns = NumberZ(:,Prof);
        plot(Z,Ns,'-','linewidth',1.5)
        grid on;
        ystr = string(Prof)+" | "+string(TimeShuK(Prof));
        xstr = string(Max75)+" | "+string(Max80)+" | "+string(Max85)+" | "+string(Max90)+" | "+string(Max95)+" | "+string(Max100)+" | "+string(Max105)+" | "+string(Max110);
        tstr = string(Aax75)+" | "+string(Aax80)+" | "+string(Aax85)+" | "+string(Aax90)+" | "+string(Aax95)+" | "+string(Aax100)+" | "+string(Aax105)+" | "+string(Aax110);
        title(tstr)
        xlim([75 115]);
        xlabel(xstr);
        ylabel(ystr);
        
        % 指定保存路径和文件名
        figPath = "D:\KYBF\TimingTrailIMG\K\"+ReadDate+"\";
        if ~exist(figPath, 'dir')
            mkdir(figPath);
        end
        figV = "K-"+string(Prof);
        print(fK, '-dpng', '-r300', figPath+figV);
        close all;
    end
end
    
    
    


fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);
    
end


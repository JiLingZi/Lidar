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

%% Na 读取

DateList = {'20231102';'20231103';...
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

for jd = 1:size(DateList,1)
% 读取光子数矩阵
ReadDate = char(DateList(jd));
RunDate = ReadDate;

Year = str2double(char(RunDate(1:4)));
Month = str2double(char(RunDate(5:6)));
Day = str2double(char(RunDate(7:8)));
DayNum = Month*30-30+Day;

CHN = 'Na';
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

NaVListNum = size(NaVList,1);
NaMListNum = size(NaMList,1);
NaZListNum = size(NaZList,1);

NumLins = NaVListNum + NaMListNum + NaZListNum;

if NumLins>0

    fldstart = ['F:\RawData\ZWDATA\MOHEnew\Na\',ReadDate,'\'];
    folder = [fldstart,'Na\'];
    files = dir(fullfile(folder, '*.dat'));
    num_files = length(files)-1;
    num_files = 3*num_files;

    FileNa = num_files;

    V_F_0 = zeros(8192, num_files/3);
    V_F_R = zeros(8192, num_files/3);
    V_F_L = zeros(8192, num_files/3);
    N_F_0 = zeros(8192, num_files/3);
    N_F_R = zeros(8192, num_files/3);
    N_F_L = zeros(8192, num_files/3);
    E_F_0 = zeros(8192, num_files/3);
    E_F_R = zeros(8192, num_files/3);
    E_F_L = zeros(8192, num_files/3);
    for j = 1:(num_files/3)
        filename = fullfile(folder, files(j).name);
        Na_table = readtable(filename);
        Na_data = table2array(Na_table);
        V_F_0(:,j) = Na_data(1:8192,2);
        V_F_R(:,j) = Na_data(1:8192,3);
        V_F_L(:,j) = Na_data(1:8192,4);
        N_F_0(:,j) = Na_data(1:8192,5);
        N_F_R(:,j) = Na_data(1:8192,6);
        N_F_L(:,j) = Na_data(1:8192,7);
        E_F_0(:,j) = Na_data(1:8192,8);
        E_F_R(:,j) = Na_data(1:8192,9);
        E_F_L(:,j) = Na_data(1:8192,10);
    end
    DenPhV = zeros(8192, num_files);
    DenPhN = zeros(8192, num_files);
    DenPhE = zeros(8192, num_files);
    for jD = 1:(num_files/3)
        DenPhV(:,3*(jD-1)+1) = V_F_0(:,jD);
        DenPhV(:,3*(jD-1)+2) = V_F_R(:,jD);
        DenPhV(:,3*(jD-1)+3) = V_F_L(:,jD);
        DenPhN(:,3*(jD-1)+1) = N_F_0(:,jD);
        DenPhN(:,3*(jD-1)+2) = N_F_R(:,jD);
        DenPhN(:,3*(jD-1)+3) = N_F_L(:,jD);
        DenPhE(:,3*(jD-1)+1) = E_F_0(:,jD);
        DenPhE(:,3*(jD-1)+2) = E_F_R(:,jD);
        DenPhE(:,3*(jD-1)+3) = E_F_L(:,jD);
    end
    TimeValues = cell(1,num_files);
    for j = 1:(num_files/3)
        filename = fullfile(folder, files(j).name);
        fid = fopen(filename);
        TXT = textscan(fid, '%s', 'Delimiter', '\n');
        RowTime = TXT{1}{4};
        TimeData = strsplit(RowTime);
        fclose(fid);
        TimeValues(:,3*(j-1)+1) = cellstr(TimeData{4});
        TimeValues(:,3*(j-1)+2) = cellstr(TimeData{4});
        TimeValues(:,3*(j-1)+3) = cellstr(TimeData{4});
    end
    jNum = NajNum;
    iNum = NaiNum;
    TimeList = cell(1,floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
    end
    TimeXNa = datetime(TimeList, 'InputFormat', ':yyyy-MM-dd''T''HH:mm:ss.SSS');
    height_num_origin = Na_data(1:8192,1);
    
    % V
    
    if NaVListNum>0
    
    PhiSum = zeros(floor(8192/iNum),num_files);
    for i = 1:floor(8192/iNum)
        PhiSum(i,:) = sum(DenPhV(iNum*(i-1)+1:iNum*i,:),1);
    end
    PhjSumNaV = zeros(floor(8192/iNum),floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        PhjSumNaV(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
    end
    AltitudeNaV = height_num_origin(1:floor(8192/iNum))*iNum;
    KM30NaV = size(AltitudeNaV(AltitudeNaV<30),1)+1;
    KM75NaV = size(AltitudeNaV(AltitudeNaV<75),1)+1;
    KM80NaV = size(AltitudeNaV(AltitudeNaV<80),1)+1;
    KM85NaV = size(AltitudeNaV(AltitudeNaV<85),1)+1;
    KM90NaV = size(AltitudeNaV(AltitudeNaV<90),1)+1;
    KM95NaV = size(AltitudeNaV(AltitudeNaV<95),1)+1;
    KM100NaV = size(AltitudeNaV(AltitudeNaV<100),1)+1;
    KM105NaV = size(AltitudeNaV(AltitudeNaV<105),1)+1;
    KM110NaV = size(AltitudeNaV(AltitudeNaV<110),1)+1;
    KM115NaV = size(AltitudeNaV(AltitudeNaV<115),1)+1;
    KM120NaV = size(AltitudeNaV(AltitudeNaV<120),1)+1;
    KM125NaV = size(AltitudeNaV(AltitudeNaV<125),1)+1;
    
    TimeShuNa = TimeXNa';

    Noise = mean(PhjSumNaV(KM120NaV:KM125NaV,:),1);

    ScatterTable = readtable('RayEffScatter.txt');
    Scatter = table2array(ScatterTable(:,2:3));
    RayScatter = Scatter(5,1);
    EffScatter = Scatter(5,2);
    [TemRay,DenRay] = atmosnrlmsise00(30000,53.3,122.7,Year,DayNum,17);
    Z = AltitudeNaV;                       % 高度矩阵
    ZR = 30;                            % 参考高度
    SigmaRay = RayScatter;              % 瑞利后向散射截面
    SigmaEff = EffScatter;              % 有效后向散射截面
    NZ = PhjSumNaV;                        % 光子数矩阵
    NB = Noise;                         % 噪声矩阵
    NZR = PhjSumNaV(KM30NaV,:);                % 参考高度处光子数
    NumberRay = sum(DenRay)*1e-6;       % 参考高度处大气模型数密度
    NumberZ = (((Z.^2).*SigmaRay.*(NZ-NB))./((ZR^2)*SigmaEff*(NZR-NB)))*NumberRay;

    Numjs = size(NaVList,1);
    
    for jpf = 1:Numjs
        Prof = ceil(NaVList(jpf)./3).*3;

        fNaV = figure('name','NaV','position',[50,100,1450,335]);
        
        subplot(1,3,1)
        
        [Max75,Aax75] = max(NumberZ(KM75NaV:KM80NaV,Prof-2));
        Aax75 = round(AltitudeNaV(Aax75+KM75NaV),2);
        AVG75 = (sum(NumberZ(KM75NaV:KM80NaV,Prof-2))-Max75)./size(NumberZ(KM75NaV:KM80NaV,Prof-2),1);
        Max75 = round(Max75-AVG75);
        [Max80,Aax80] = max(NumberZ(KM80NaV:KM85NaV,Prof-2));
        Aax80 = round(AltitudeNaV(Aax80+KM80NaV),2);
        AVG80 = (sum(NumberZ(KM80NaV:KM85NaV,Prof-2))-Max80)./size(NumberZ(KM80NaV:KM85NaV,Prof-2),1);
        Max80 = round(Max80-AVG80);
        [Max85,Aax85] = max(NumberZ(KM85NaV:KM90NaV,Prof-2));
        Aax85 = round(AltitudeNaV(Aax85+KM85NaV),2);
        AVG85 = (sum(NumberZ(KM85NaV:KM90NaV,Prof-2))-Max85)./size(NumberZ(KM85NaV:KM90NaV,Prof-2),1);
        Max85 = round(Max85-AVG85);
        [Max90,Aax90] = max(NumberZ(KM90NaV:KM95NaV,Prof-2));
        Aax90 = round(AltitudeNaV(Aax90+KM90NaV),2);
        AVG90 = (sum(NumberZ(KM90NaV:KM95NaV,Prof-2))-Max90)./size(NumberZ(KM90NaV:KM95NaV,Prof-2),1);
        Max90 = round(Max90-AVG90);
        [Max95,Aax95] = max(NumberZ(KM95NaV:KM100NaV,Prof-2));
        Aax95 = round(AltitudeNaV(Aax95+KM95NaV),2);
        AVG95 = (sum(NumberZ(KM95NaV:KM100NaV,Prof-2))-Max95)./size(NumberZ(KM95NaV:KM100NaV,Prof-2),1);
        Max95 = round(Max95-AVG95);
        [Max100,Aax100] = max(NumberZ(KM100NaV:KM105NaV,Prof-2));
        Aax100 = round(AltitudeNaV(Aax100+KM100NaV),2);
        AVG100 = (sum(NumberZ(KM100NaV:KM105NaV,Prof-2))-Max100)./size(NumberZ(KM100NaV:KM105NaV,Prof-2),1);
        Max100 = round(Max100-AVG100);
        [Max105,Aax105] = max(NumberZ(KM105NaV:KM110NaV,Prof-2));
        Aax105 = round(AltitudeNaV(Aax105+KM105NaV),2);
        AVG105 = (sum(NumberZ(KM105NaV:KM110NaV,Prof-2))-Max105)./size(NumberZ(KM105NaV:KM110NaV,Prof-2),1);
        Max105 = round(Max105-AVG105);
        [Max110,Aax110] = max(NumberZ(KM110NaV:KM115NaV,Prof-2));
        Aax110 = round(AltitudeNaV(Aax110+KM110NaV),2);
        AVG110 = (sum(NumberZ(KM110NaV:KM115NaV,Prof-2))-Max110)./size(NumberZ(KM110NaV:KM115NaV,Prof-2),1);
        Max110 = round(Max110-AVG110);
        
        Ns = NumberZ(:,Prof-2);
        plot(Z,Ns,'-','linewidth',1.5)
        grid on;
        ystr = string(Prof-2)+" | "+string(TimeShuNa(Prof-2));
        xstr = string(Max75)+" | "+string(Max80)+" | "+string(Max85)+" | "+string(Max90)+" | "+string(Max95)+" | "+string(Max100)+" | "+string(Max105)+" | "+string(Max110);
        tstr = string(Aax75)+" | "+string(Aax80)+" | "+string(Aax85)+" | "+string(Aax90)+" | "+string(Aax95)+" | "+string(Aax100)+" | "+string(Aax105)+" | "+string(Aax110);
        title(tstr)
        xlim([75 115]);
        xlabel(xstr);
        ylabel(ystr);
        
        subplot(1,3,2)
        
        [Max75,Aax75] = max(NumberZ(KM75NaV:KM80NaV,Prof-1));
        Aax75 = round(AltitudeNaV(Aax75+KM75NaV),2);
        AVG75 = (sum(NumberZ(KM75NaV:KM80NaV,Prof-1))-Max75)./size(NumberZ(KM75NaV:KM80NaV,Prof-1),1);
        Max75 = round(Max75-AVG75);
        [Max80,Aax80] = max(NumberZ(KM80NaV:KM85NaV,Prof-1));
        Aax80 = round(AltitudeNaV(Aax80+KM80NaV),2);
        AVG80 = (sum(NumberZ(KM80NaV:KM85NaV,Prof-1))-Max80)./size(NumberZ(KM80NaV:KM85NaV,Prof-1),1);
        Max80 = round(Max80-AVG80);
        [Max85,Aax85] = max(NumberZ(KM85NaV:KM90NaV,Prof-1));
        Aax85 = round(AltitudeNaV(Aax85+KM85NaV),2);
        AVG85 = (sum(NumberZ(KM85NaV:KM90NaV,Prof-1))-Max85)./size(NumberZ(KM85NaV:KM90NaV,Prof-1),1);
        Max85 = round(Max85-AVG85);
        [Max90,Aax90] = max(NumberZ(KM90NaV:KM95NaV,Prof-1));
        Aax90 = round(AltitudeNaV(Aax90+KM90NaV),2);
        AVG90 = (sum(NumberZ(KM90NaV:KM95NaV,Prof-1))-Max90)./size(NumberZ(KM90NaV:KM95NaV,Prof-1),1);
        Max90 = round(Max90-AVG90);
        [Max95,Aax95] = max(NumberZ(KM95NaV:KM100NaV,Prof-1));
        Aax95 = round(AltitudeNaV(Aax95+KM95NaV),2);
        AVG95 = (sum(NumberZ(KM95NaV:KM100NaV,Prof-1))-Max95)./size(NumberZ(KM95NaV:KM100NaV,Prof-1),1);
        Max95 = round(Max95-AVG95);
        [Max100,Aax100] = max(NumberZ(KM100NaV:KM105NaV,Prof-1));
        Aax100 = round(AltitudeNaV(Aax100+KM100NaV),2);
        AVG100 = (sum(NumberZ(KM100NaV:KM105NaV,Prof-1))-Max100)./size(NumberZ(KM100NaV:KM105NaV,Prof-1),1);
        Max100 = round(Max100-AVG100);
        [Max105,Aax105] = max(NumberZ(KM105NaV:KM110NaV,Prof-1));
        Aax105 = round(AltitudeNaV(Aax105+KM105NaV),2);
        AVG105 = (sum(NumberZ(KM105NaV:KM110NaV,Prof-1))-Max105)./size(NumberZ(KM105NaV:KM110NaV,Prof-1),1);
        Max105 = round(Max105-AVG105);
        [Max110,Aax110] = max(NumberZ(KM110NaV:KM115NaV,Prof-1));
        Aax110 = round(AltitudeNaV(Aax110+KM110NaV),2);
        AVG110 = (sum(NumberZ(KM110NaV:KM115NaV,Prof-1))-Max110)./size(NumberZ(KM110NaV:KM115NaV,Prof-1),1);
        Max110 = round(Max110-AVG110);
        
        Ns = NumberZ(:,Prof-1);
        plot(Z,Ns,'-','linewidth',1.5)
        grid on;
        ystr = string(Prof-1)+" | "+string(TimeShuNa(Prof-1));
        xstr = string(Max75)+" | "+string(Max80)+" | "+string(Max85)+" | "+string(Max90)+" | "+string(Max95)+" | "+string(Max100)+" | "+string(Max105)+" | "+string(Max110);
        tstr = string(Aax75)+" | "+string(Aax80)+" | "+string(Aax85)+" | "+string(Aax90)+" | "+string(Aax95)+" | "+string(Aax100)+" | "+string(Aax105)+" | "+string(Aax110);
        title(tstr)
        xlim([75 115]);
        xlabel(xstr);
        ylabel(ystr);
        
        subplot(1,3,3)
        
        [Max75,Aax75] = max(NumberZ(KM75NaV:KM80NaV,Prof));
        Aax75 = round(AltitudeNaV(Aax75+KM75NaV),2);
        AVG75 = (sum(NumberZ(KM75NaV:KM80NaV,Prof))-Max75)./size(NumberZ(KM75NaV:KM80NaV,Prof),1);
        Max75 = round(Max75-AVG75);
        [Max80,Aax80] = max(NumberZ(KM80NaV:KM85NaV,Prof));
        Aax80 = round(AltitudeNaV(Aax80+KM80NaV),2);
        AVG80 = (sum(NumberZ(KM80NaV:KM85NaV,Prof))-Max80)./size(NumberZ(KM80NaV:KM85NaV,Prof),1);
        Max80 = round(Max80-AVG80);
        [Max85,Aax85] = max(NumberZ(KM85NaV:KM90NaV,Prof));
        Aax85 = round(AltitudeNaV(Aax85+KM85NaV),2);
        AVG85 = (sum(NumberZ(KM85NaV:KM90NaV,Prof))-Max85)./size(NumberZ(KM85NaV:KM90NaV,Prof),1);
        Max85 = round(Max85-AVG85);
        [Max90,Aax90] = max(NumberZ(KM90NaV:KM95NaV,Prof));
        Aax90 = round(AltitudeNaV(Aax90+KM90NaV),2);
        AVG90 = (sum(NumberZ(KM90NaV:KM95NaV,Prof))-Max90)./size(NumberZ(KM90NaV:KM95NaV,Prof),1);
        Max90 = round(Max90-AVG90);
        [Max95,Aax95] = max(NumberZ(KM95NaV:KM100NaV,Prof));
        Aax95 = round(AltitudeNaV(Aax95+KM95NaV),2);
        AVG95 = (sum(NumberZ(KM95NaV:KM100NaV,Prof))-Max95)./size(NumberZ(KM95NaV:KM100NaV,Prof),1);
        Max95 = round(Max95-AVG95);
        [Max100,Aax100] = max(NumberZ(KM100NaV:KM105NaV,Prof));
        Aax100 = round(AltitudeNaV(Aax100+KM100NaV),2);
        AVG100 = (sum(NumberZ(KM100NaV:KM105NaV,Prof))-Max100)./size(NumberZ(KM100NaV:KM105NaV,Prof),1);
        Max100 = round(Max100-AVG100);
        [Max105,Aax105] = max(NumberZ(KM105NaV:KM110NaV,Prof));
        Aax105 = round(AltitudeNaV(Aax105+KM105NaV),2);
        AVG105 = (sum(NumberZ(KM105NaV:KM110NaV,Prof))-Max105)./size(NumberZ(KM105NaV:KM110NaV,Prof),1);
        Max105 = round(Max105-AVG105);
        [Max110,Aax110] = max(NumberZ(KM110NaV:KM115NaV,Prof));
        Aax110 = round(AltitudeNaV(Aax110+KM110NaV),2);
        AVG110 = (sum(NumberZ(KM110NaV:KM115NaV,Prof))-Max110)./size(NumberZ(KM110NaV:KM115NaV,Prof),1);
        Max110 = round(Max110-AVG110);
        
        Ns = NumberZ(:,Prof);
        plot(Z,Ns,'-','linewidth',1.5)
        grid on;
        ystr = string(Prof)+" | "+string(TimeShuNa(Prof));
        xstr = string(Max75)+" | "+string(Max80)+" | "+string(Max85)+" | "+string(Max90)+" | "+string(Max95)+" | "+string(Max100)+" | "+string(Max105)+" | "+string(Max110);
        tstr = string(Aax75)+" | "+string(Aax80)+" | "+string(Aax85)+" | "+string(Aax90)+" | "+string(Aax95)+" | "+string(Aax100)+" | "+string(Aax105)+" | "+string(Aax110);
        title(tstr)
        xlim([75 115]);
        xlabel(xstr);
        ylabel(ystr);
        
        % 指定保存路径和文件名
        figPath = "D:\KYBF\TimingTrailIMG\Na\"+ReadDate+"\";
        if ~exist(figPath, 'dir')
            mkdir(figPath);
        end
        figV = "Na-V-"+string(Prof);
        print(fNaV, '-dpng', '-r300', figPath+figV);
        close all;
    end
    end
    
    % N
    
    if NaMListNum>0
    
    PhiSum = zeros(floor(8192/iNum),num_files);
    for i = 1:floor(8192/iNum)
        PhiSum(i,:) = sum(DenPhN(iNum*(i-1)+1:iNum*i,:),1);
    end
    PhjSumNaM = zeros(floor(8192/iNum),floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        PhjSumNaM(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
    end
    AltitudeNaM = height_num_origin(1:floor(8192/iNum))*iNum*(sqrt(3)/2);
    KM30NaM = size(AltitudeNaM(AltitudeNaM<30),1)+1;
    KM75NaM = size(AltitudeNaM(AltitudeNaM<75),1)+1;
    KM80NaM = size(AltitudeNaM(AltitudeNaM<80),1)+1;
    KM85NaM = size(AltitudeNaM(AltitudeNaM<85),1)+1;
    KM90NaM = size(AltitudeNaM(AltitudeNaM<90),1)+1;
    KM95NaM = size(AltitudeNaM(AltitudeNaM<95),1)+1;
    KM100NaM = size(AltitudeNaM(AltitudeNaM<100),1)+1;
    KM105NaM = size(AltitudeNaM(AltitudeNaM<105),1)+1;
    KM110NaM = size(AltitudeNaM(AltitudeNaM<110),1)+1;
    KM115NaM = size(AltitudeNaM(AltitudeNaM<115),1)+1;
    KM120NaM = size(AltitudeNaM(AltitudeNaM<120),1)+1;
    KM125NaM = size(AltitudeNaM(AltitudeNaM<125),1)+1;
    
    TimeShuNa = TimeXNa';

    Noise = mean(PhjSumNaM(KM120NaM:KM125NaM,:),1);

    ScatterTable = readtable('RayEffScatter.txt');
    Scatter = table2array(ScatterTable(:,2:3));
    RayScatter = Scatter(5,1);
    EffScatter = Scatter(5,2);
    [TemRay,DenRay] = atmosnrlmsise00(30000,53.3,122.7,Year,DayNum,17);
    Z = AltitudeNaM;                       % 高度矩阵
    ZR = 30;                            % 参考高度
    SigmaRay = RayScatter;              % 瑞利后向散射截面
    SigmaEff = EffScatter;              % 有效后向散射截面
    NZ = PhjSumNaM;                        % 光子数矩阵
    NB = Noise;                         % 噪声矩阵
    NZR = PhjSumNaM(KM30NaM,:);                % 参考高度处光子数
    NumberRay = sum(DenRay)*1e-6;       % 参考高度处大气模型数密度
    NumberZ = (((Z.^2).*SigmaRay.*(NZ-NB))./((ZR^2)*SigmaEff*(NZR-NB)))*NumberRay;

    Numjs = size(NaMList,1);
    
    for jpf = 1:Numjs
        Prof = ceil(NaMList(jpf)./3).*3;

        fNaM = figure('name','NaM','position',[50,100,1450,335]);
        
        subplot(1,3,1)
        
        [Max75,Aax75] = max(NumberZ(KM75NaM:KM80NaM,Prof-2));
        Aax75 = round(AltitudeNaM(Aax75+KM75NaM),2);
        AMG75 = (sum(NumberZ(KM75NaM:KM80NaM,Prof-2))-Max75)./size(NumberZ(KM75NaM:KM80NaM,Prof-2),1);
        Max75 = round(Max75-AMG75);
        [Max80,Aax80] = max(NumberZ(KM80NaM:KM85NaM,Prof-2));
        Aax80 = round(AltitudeNaM(Aax80+KM80NaM),2);
        AMG80 = (sum(NumberZ(KM80NaM:KM85NaM,Prof-2))-Max80)./size(NumberZ(KM80NaM:KM85NaM,Prof-2),1);
        Max80 = round(Max80-AMG80);
        [Max85,Aax85] = max(NumberZ(KM85NaM:KM90NaM,Prof-2));
        Aax85 = round(AltitudeNaM(Aax85+KM85NaM),2);
        AMG85 = (sum(NumberZ(KM85NaM:KM90NaM,Prof-2))-Max85)./size(NumberZ(KM85NaM:KM90NaM,Prof-2),1);
        Max85 = round(Max85-AMG85);
        [Max90,Aax90] = max(NumberZ(KM90NaM:KM95NaM,Prof-2));
        Aax90 = round(AltitudeNaM(Aax90+KM90NaM),2);
        AMG90 = (sum(NumberZ(KM90NaM:KM95NaM,Prof-2))-Max90)./size(NumberZ(KM90NaM:KM95NaM,Prof-2),1);
        Max90 = round(Max90-AMG90);
        [Max95,Aax95] = max(NumberZ(KM95NaM:KM100NaM,Prof-2));
        Aax95 = round(AltitudeNaM(Aax95+KM95NaM),2);
        AMG95 = (sum(NumberZ(KM95NaM:KM100NaM,Prof-2))-Max95)./size(NumberZ(KM95NaM:KM100NaM,Prof-2),1);
        Max95 = round(Max95-AMG95);
        [Max100,Aax100] = max(NumberZ(KM100NaM:KM105NaM,Prof-2));
        Aax100 = round(AltitudeNaM(Aax100+KM100NaM),2);
        AMG100 = (sum(NumberZ(KM100NaM:KM105NaM,Prof-2))-Max100)./size(NumberZ(KM100NaM:KM105NaM,Prof-2),1);
        Max100 = round(Max100-AMG100);
        [Max105,Aax105] = max(NumberZ(KM105NaM:KM110NaM,Prof-2));
        Aax105 = round(AltitudeNaM(Aax105+KM105NaM),2);
        AMG105 = (sum(NumberZ(KM105NaM:KM110NaM,Prof-2))-Max105)./size(NumberZ(KM105NaM:KM110NaM,Prof-2),1);
        Max105 = round(Max105-AMG105);
        [Max110,Aax110] = max(NumberZ(KM110NaM:KM115NaM,Prof-2));
        Aax110 = round(AltitudeNaM(Aax110+KM110NaM),2);
        AMG110 = (sum(NumberZ(KM110NaM:KM115NaM,Prof-2))-Max110)./size(NumberZ(KM110NaM:KM115NaM,Prof-2),1);
        Max110 = round(Max110-AMG110);
        
        Ns = NumberZ(:,Prof-2);
        plot(Z,Ns,'-','linewidth',1.5)
        grid on;
        ystr = string(Prof-2)+" | "+string(TimeShuNa(Prof-2));
        xstr = string(Max75)+" | "+string(Max80)+" | "+string(Max85)+" | "+string(Max90)+" | "+string(Max95)+" | "+string(Max100)+" | "+string(Max105)+" | "+string(Max110);
        tstr = string(Aax75)+" | "+string(Aax80)+" | "+string(Aax85)+" | "+string(Aax90)+" | "+string(Aax95)+" | "+string(Aax100)+" | "+string(Aax105)+" | "+string(Aax110);
        title(tstr)
        xlim([75 115]);
        xlabel(xstr);
        ylabel(ystr);
        
        subplot(1,3,2)
        
        [Max75,Aax75] = max(NumberZ(KM75NaM:KM80NaM,Prof-1));
        Aax75 = round(AltitudeNaM(Aax75+KM75NaM),2);
        AMG75 = (sum(NumberZ(KM75NaM:KM80NaM,Prof-1))-Max75)./size(NumberZ(KM75NaM:KM80NaM,Prof-1),1);
        Max75 = round(Max75-AMG75);
        [Max80,Aax80] = max(NumberZ(KM80NaM:KM85NaM,Prof-1));
        Aax80 = round(AltitudeNaM(Aax80+KM80NaM),2);
        AMG80 = (sum(NumberZ(KM80NaM:KM85NaM,Prof-1))-Max80)./size(NumberZ(KM80NaM:KM85NaM,Prof-1),1);
        Max80 = round(Max80-AMG80);
        [Max85,Aax85] = max(NumberZ(KM85NaM:KM90NaM,Prof-1));
        Aax85 = round(AltitudeNaM(Aax85+KM85NaM),2);
        AMG85 = (sum(NumberZ(KM85NaM:KM90NaM,Prof-1))-Max85)./size(NumberZ(KM85NaM:KM90NaM,Prof-1),1);
        Max85 = round(Max85-AMG85);
        [Max90,Aax90] = max(NumberZ(KM90NaM:KM95NaM,Prof-1));
        Aax90 = round(AltitudeNaM(Aax90+KM90NaM),2);
        AMG90 = (sum(NumberZ(KM90NaM:KM95NaM,Prof-1))-Max90)./size(NumberZ(KM90NaM:KM95NaM,Prof-1),1);
        Max90 = round(Max90-AMG90);
        [Max95,Aax95] = max(NumberZ(KM95NaM:KM100NaM,Prof-1));
        Aax95 = round(AltitudeNaM(Aax95+KM95NaM),2);
        AMG95 = (sum(NumberZ(KM95NaM:KM100NaM,Prof-1))-Max95)./size(NumberZ(KM95NaM:KM100NaM,Prof-1),1);
        Max95 = round(Max95-AMG95);
        [Max100,Aax100] = max(NumberZ(KM100NaM:KM105NaM,Prof-1));
        Aax100 = round(AltitudeNaM(Aax100+KM100NaM),2);
        AMG100 = (sum(NumberZ(KM100NaM:KM105NaM,Prof-1))-Max100)./size(NumberZ(KM100NaM:KM105NaM,Prof-1),1);
        Max100 = round(Max100-AMG100);
        [Max105,Aax105] = max(NumberZ(KM105NaM:KM110NaM,Prof-1));
        Aax105 = round(AltitudeNaM(Aax105+KM105NaM),2);
        AMG105 = (sum(NumberZ(KM105NaM:KM110NaM,Prof-1))-Max105)./size(NumberZ(KM105NaM:KM110NaM,Prof-1),1);
        Max105 = round(Max105-AMG105);
        [Max110,Aax110] = max(NumberZ(KM110NaM:KM115NaM,Prof-1));
        Aax110 = round(AltitudeNaM(Aax110+KM110NaM),2);
        AMG110 = (sum(NumberZ(KM110NaM:KM115NaM,Prof-1))-Max110)./size(NumberZ(KM110NaM:KM115NaM,Prof-1),1);
        Max110 = round(Max110-AMG110);
        
        Ns = NumberZ(:,Prof-1);
        plot(Z,Ns,'-','linewidth',1.5)
        grid on;
        ystr = string(Prof-1)+" | "+string(TimeShuNa(Prof-1));
        xstr = string(Max75)+" | "+string(Max80)+" | "+string(Max85)+" | "+string(Max90)+" | "+string(Max95)+" | "+string(Max100)+" | "+string(Max105)+" | "+string(Max110);
        tstr = string(Aax75)+" | "+string(Aax80)+" | "+string(Aax85)+" | "+string(Aax90)+" | "+string(Aax95)+" | "+string(Aax100)+" | "+string(Aax105)+" | "+string(Aax110);
        title(tstr)
        xlim([75 115]);
        xlabel(xstr);
        ylabel(ystr);
        
        subplot(1,3,3)
        
        [Max75,Aax75] = max(NumberZ(KM75NaM:KM80NaM,Prof));
        Aax75 = round(AltitudeNaM(Aax75+KM75NaM),2);
        AMG75 = (sum(NumberZ(KM75NaM:KM80NaM,Prof))-Max75)./size(NumberZ(KM75NaM:KM80NaM,Prof),1);
        Max75 = round(Max75-AMG75);
        [Max80,Aax80] = max(NumberZ(KM80NaM:KM85NaM,Prof));
        Aax80 = round(AltitudeNaM(Aax80+KM80NaM),2);
        AMG80 = (sum(NumberZ(KM80NaM:KM85NaM,Prof))-Max80)./size(NumberZ(KM80NaM:KM85NaM,Prof),1);
        Max80 = round(Max80-AMG80);
        [Max85,Aax85] = max(NumberZ(KM85NaM:KM90NaM,Prof));
        Aax85 = round(AltitudeNaM(Aax85+KM85NaM),2);
        AMG85 = (sum(NumberZ(KM85NaM:KM90NaM,Prof))-Max85)./size(NumberZ(KM85NaM:KM90NaM,Prof),1);
        Max85 = round(Max85-AMG85);
        [Max90,Aax90] = max(NumberZ(KM90NaM:KM95NaM,Prof));
        Aax90 = round(AltitudeNaM(Aax90+KM90NaM),2);
        AMG90 = (sum(NumberZ(KM90NaM:KM95NaM,Prof))-Max90)./size(NumberZ(KM90NaM:KM95NaM,Prof),1);
        Max90 = round(Max90-AMG90);
        [Max95,Aax95] = max(NumberZ(KM95NaM:KM100NaM,Prof));
        Aax95 = round(AltitudeNaM(Aax95+KM95NaM),2);
        AMG95 = (sum(NumberZ(KM95NaM:KM100NaM,Prof))-Max95)./size(NumberZ(KM95NaM:KM100NaM,Prof),1);
        Max95 = round(Max95-AMG95);
        [Max100,Aax100] = max(NumberZ(KM100NaM:KM105NaM,Prof));
        Aax100 = round(AltitudeNaM(Aax100+KM100NaM),2);
        AMG100 = (sum(NumberZ(KM100NaM:KM105NaM,Prof))-Max100)./size(NumberZ(KM100NaM:KM105NaM,Prof),1);
        Max100 = round(Max100-AMG100);
        [Max105,Aax105] = max(NumberZ(KM105NaM:KM110NaM,Prof));
        Aax105 = round(AltitudeNaM(Aax105+KM105NaM),2);
        AMG105 = (sum(NumberZ(KM105NaM:KM110NaM,Prof))-Max105)./size(NumberZ(KM105NaM:KM110NaM,Prof),1);
        Max105 = round(Max105-AMG105);
        [Max110,Aax110] = max(NumberZ(KM110NaM:KM115NaM,Prof));
        Aax110 = round(AltitudeNaM(Aax110+KM110NaM),2);
        AMG110 = (sum(NumberZ(KM110NaM:KM115NaM,Prof))-Max110)./size(NumberZ(KM110NaM:KM115NaM,Prof),1);
        Max110 = round(Max110-AMG110);
        
        Ns = NumberZ(:,Prof);
        plot(Z,Ns,'-','linewidth',1.5)
        grid on;
        ystr = string(Prof)+" | "+string(TimeShuNa(Prof));
        xstr = string(Max75)+" | "+string(Max80)+" | "+string(Max85)+" | "+string(Max90)+" | "+string(Max95)+" | "+string(Max100)+" | "+string(Max105)+" | "+string(Max110);
        tstr = string(Aax75)+" | "+string(Aax80)+" | "+string(Aax85)+" | "+string(Aax90)+" | "+string(Aax95)+" | "+string(Aax100)+" | "+string(Aax105)+" | "+string(Aax110);
        title(tstr)
        xlim([75 115]);
        xlabel(xstr);
        ylabel(ystr);
        
        % 指定保存路径和文件名
        figPath = "D:\KYBF\TimingTrailIMG\Na\"+ReadDate+"\";
        if ~exist(figPath, 'dir')
            mkdir(figPath);
        end
        figV = "Na-M-"+string(Prof);
        print(fNaM, '-dpng', '-r300', figPath+figV);
        close all;
    end
    end
    
    % E
    
    if NaZListNum>0
    
    PhiSum = zeros(floor(8192/iNum),num_files);
    for i = 1:floor(8192/iNum)
        PhiSum(i,:) = sum(DenPhE(iNum*(i-1)+1:iNum*i,:),1);
    end
    PhjSumNaZ = zeros(floor(8192/iNum),floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        PhjSumNaZ(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
    end
    AltitudeNaZ = height_num_origin(1:floor(8192/iNum))*iNum*(sqrt(3)/2);
    KM30NaZ = size(AltitudeNaZ(AltitudeNaZ<30),1)+1;
    KM75NaZ = size(AltitudeNaZ(AltitudeNaZ<75),1)+1;
    KM80NaZ = size(AltitudeNaZ(AltitudeNaZ<80),1)+1;
    KM85NaZ = size(AltitudeNaZ(AltitudeNaZ<85),1)+1;
    KM90NaZ = size(AltitudeNaZ(AltitudeNaZ<90),1)+1;
    KM95NaZ = size(AltitudeNaZ(AltitudeNaZ<95),1)+1;
    KM100NaZ = size(AltitudeNaZ(AltitudeNaZ<100),1)+1;
    KM105NaZ = size(AltitudeNaZ(AltitudeNaZ<105),1)+1;
    KM110NaZ = size(AltitudeNaZ(AltitudeNaZ<110),1)+1;
    KM115NaZ = size(AltitudeNaZ(AltitudeNaZ<115),1)+1;
    KM120NaZ = size(AltitudeNaZ(AltitudeNaZ<120),1)+1;
    KM125NaZ = size(AltitudeNaZ(AltitudeNaZ<125),1)+1;

    TimeShuNa = TimeXNa';

    Noise = mean(PhjSumNaZ(KM120NaZ:KM125NaZ,:),1);

    ScatterTable = readtable('RayEffScatter.txt');
    Scatter = table2array(ScatterTable(:,2:3));
    RayScatter = Scatter(5,1);
    EffScatter = Scatter(5,2);
    [TemRay,DenRay] = atmosnrlmsise00(30000,53.3,122.7,Year,DayNum,17);
    Z = AltitudeNaZ;                       % 高度矩阵
    ZR = 30;                            % 参考高度
    SigmaRay = RayScatter;              % 瑞利后向散射截面
    SigmaEff = EffScatter;              % 有效后向散射截面
    NZ = PhjSumNaZ;                        % 光子数矩阵
    NB = Noise;                         % 噪声矩阵
    NZR = PhjSumNaZ(KM30NaZ,:);                % 参考高度处光子数
    NumberRay = sum(DenRay)*1e-6;       % 参考高度处大气模型数密度
    NumberZ = (((Z.^2).*SigmaRay.*(NZ-NB))./((ZR^2)*SigmaEff*(NZR-NB)))*NumberRay;

    Numjs = size(NaZList,1);
    
    for jpf = 1:Numjs
        Prof = ceil(NaZList(jpf)./3).*3;

        fNaZ = figure('name','NaZ','position',[50,100,1450,335]);
        
        subplot(1,3,1)
        
        [Max75,Aax75] = max(NumberZ(KM75NaZ:KM80NaZ,Prof-2));
        Aax75 = round(AltitudeNaZ(Aax75+KM75NaZ),2);
        AZG75 = (sum(NumberZ(KM75NaZ:KM80NaZ,Prof-2))-Max75)./size(NumberZ(KM75NaZ:KM80NaZ,Prof-2),1);
        Max75 = round(Max75-AZG75);
        [Max80,Aax80] = max(NumberZ(KM80NaZ:KM85NaZ,Prof-2));
        Aax80 = round(AltitudeNaZ(Aax80+KM80NaZ),2);
        AZG80 = (sum(NumberZ(KM80NaZ:KM85NaZ,Prof-2))-Max80)./size(NumberZ(KM80NaZ:KM85NaZ,Prof-2),1);
        Max80 = round(Max80-AZG80);
        [Max85,Aax85] = max(NumberZ(KM85NaZ:KM90NaZ,Prof-2));
        Aax85 = round(AltitudeNaZ(Aax85+KM85NaZ),2);
        AZG85 = (sum(NumberZ(KM85NaZ:KM90NaZ,Prof-2))-Max85)./size(NumberZ(KM85NaZ:KM90NaZ,Prof-2),1);
        Max85 = round(Max85-AZG85);
        [Max90,Aax90] = max(NumberZ(KM90NaZ:KM95NaZ,Prof-2));
        Aax90 = round(AltitudeNaZ(Aax90+KM90NaZ),2);
        AZG90 = (sum(NumberZ(KM90NaZ:KM95NaZ,Prof-2))-Max90)./size(NumberZ(KM90NaZ:KM95NaZ,Prof-2),1);
        Max90 = round(Max90-AZG90);
        [Max95,Aax95] = max(NumberZ(KM95NaZ:KM100NaZ,Prof-2));
        Aax95 = round(AltitudeNaZ(Aax95+KM95NaZ),2);
        AZG95 = (sum(NumberZ(KM95NaZ:KM100NaZ,Prof-2))-Max95)./size(NumberZ(KM95NaZ:KM100NaZ,Prof-2),1);
        Max95 = round(Max95-AZG95);
        [Max100,Aax100] = max(NumberZ(KM100NaZ:KM105NaZ,Prof-2));
        Aax100 = round(AltitudeNaZ(Aax100+KM100NaZ),2);
        AZG100 = (sum(NumberZ(KM100NaZ:KM105NaZ,Prof-2))-Max100)./size(NumberZ(KM100NaZ:KM105NaZ,Prof-2),1);
        Max100 = round(Max100-AZG100);
        [Max105,Aax105] = max(NumberZ(KM105NaZ:KM110NaZ,Prof-2));
        Aax105 = round(AltitudeNaZ(Aax105+KM105NaZ),2);
        AZG105 = (sum(NumberZ(KM105NaZ:KM110NaZ,Prof-2))-Max105)./size(NumberZ(KM105NaZ:KM110NaZ,Prof-2),1);
        Max105 = round(Max105-AZG105);
        [Max110,Aax110] = max(NumberZ(KM110NaZ:KM115NaZ,Prof-2));
        Aax110 = round(AltitudeNaZ(Aax110+KM110NaZ),2);
        AZG110 = (sum(NumberZ(KM110NaZ:KM115NaZ,Prof-2))-Max110)./size(NumberZ(KM110NaZ:KM115NaZ,Prof-2),1);
        Max110 = round(Max110-AZG110);
        
        Ns = NumberZ(:,Prof-2);
        plot(Z,Ns,'-','linewidth',1.5)
        grid on;
        ystr = string(Prof-2)+" | "+string(TimeShuNa(Prof-2));
        xstr = string(Max75)+" | "+string(Max80)+" | "+string(Max85)+" | "+string(Max90)+" | "+string(Max95)+" | "+string(Max100)+" | "+string(Max105)+" | "+string(Max110);
        tstr = string(Aax75)+" | "+string(Aax80)+" | "+string(Aax85)+" | "+string(Aax90)+" | "+string(Aax95)+" | "+string(Aax100)+" | "+string(Aax105)+" | "+string(Aax110);
        title(tstr)
        xlim([75 115]);
        xlabel(xstr);
        ylabel(ystr);
        
        subplot(1,3,2)
        
        [Max75,Aax75] = max(NumberZ(KM75NaZ:KM80NaZ,Prof-1));
        Aax75 = round(AltitudeNaZ(Aax75+KM75NaZ),2);
        AZG75 = (sum(NumberZ(KM75NaZ:KM80NaZ,Prof-1))-Max75)./size(NumberZ(KM75NaZ:KM80NaZ,Prof-1),1);
        Max75 = round(Max75-AZG75);
        [Max80,Aax80] = max(NumberZ(KM80NaZ:KM85NaZ,Prof-1));
        Aax80 = round(AltitudeNaZ(Aax80+KM80NaZ),2);
        AZG80 = (sum(NumberZ(KM80NaZ:KM85NaZ,Prof-1))-Max80)./size(NumberZ(KM80NaZ:KM85NaZ,Prof-1),1);
        Max80 = round(Max80-AZG80);
        [Max85,Aax85] = max(NumberZ(KM85NaZ:KM90NaZ,Prof-1));
        Aax85 = round(AltitudeNaZ(Aax85+KM85NaZ),2);
        AZG85 = (sum(NumberZ(KM85NaZ:KM90NaZ,Prof-1))-Max85)./size(NumberZ(KM85NaZ:KM90NaZ,Prof-1),1);
        Max85 = round(Max85-AZG85);
        [Max90,Aax90] = max(NumberZ(KM90NaZ:KM95NaZ,Prof-1));
        Aax90 = round(AltitudeNaZ(Aax90+KM90NaZ),2);
        AZG90 = (sum(NumberZ(KM90NaZ:KM95NaZ,Prof-1))-Max90)./size(NumberZ(KM90NaZ:KM95NaZ,Prof-1),1);
        Max90 = round(Max90-AZG90);
        [Max95,Aax95] = max(NumberZ(KM95NaZ:KM100NaZ,Prof-1));
        Aax95 = round(AltitudeNaZ(Aax95+KM95NaZ),2);
        AZG95 = (sum(NumberZ(KM95NaZ:KM100NaZ,Prof-1))-Max95)./size(NumberZ(KM95NaZ:KM100NaZ,Prof-1),1);
        Max95 = round(Max95-AZG95);
        [Max100,Aax100] = max(NumberZ(KM100NaZ:KM105NaZ,Prof-1));
        Aax100 = round(AltitudeNaZ(Aax100+KM100NaZ),2);
        AZG100 = (sum(NumberZ(KM100NaZ:KM105NaZ,Prof-1))-Max100)./size(NumberZ(KM100NaZ:KM105NaZ,Prof-1),1);
        Max100 = round(Max100-AZG100);
        [Max105,Aax105] = max(NumberZ(KM105NaZ:KM110NaZ,Prof-1));
        Aax105 = round(AltitudeNaZ(Aax105+KM105NaZ),2);
        AZG105 = (sum(NumberZ(KM105NaZ:KM110NaZ,Prof-1))-Max105)./size(NumberZ(KM105NaZ:KM110NaZ,Prof-1),1);
        Max105 = round(Max105-AZG105);
        [Max110,Aax110] = max(NumberZ(KM110NaZ:KM115NaZ,Prof-1));
        Aax110 = round(AltitudeNaZ(Aax110+KM110NaZ),2);
        AZG110 = (sum(NumberZ(KM110NaZ:KM115NaZ,Prof-1))-Max110)./size(NumberZ(KM110NaZ:KM115NaZ,Prof-1),1);
        Max110 = round(Max110-AZG110);
        
        Ns = NumberZ(:,Prof-1);
        plot(Z,Ns,'-','linewidth',1.5)
        grid on;
        ystr = string(Prof-1)+" | "+string(TimeShuNa(Prof-1));
        xstr = string(Max75)+" | "+string(Max80)+" | "+string(Max85)+" | "+string(Max90)+" | "+string(Max95)+" | "+string(Max100)+" | "+string(Max105)+" | "+string(Max110);
        tstr = string(Aax75)+" | "+string(Aax80)+" | "+string(Aax85)+" | "+string(Aax90)+" | "+string(Aax95)+" | "+string(Aax100)+" | "+string(Aax105)+" | "+string(Aax110);
        title(tstr)
        xlim([75 115]);
        xlabel(xstr);
        ylabel(ystr);
        
        subplot(1,3,3)
        
        [Max75,Aax75] = max(NumberZ(KM75NaZ:KM80NaZ,Prof));
        Aax75 = round(AltitudeNaZ(Aax75+KM75NaZ),2);
        AZG75 = (sum(NumberZ(KM75NaZ:KM80NaZ,Prof))-Max75)./size(NumberZ(KM75NaZ:KM80NaZ,Prof),1);
        Max75 = round(Max75-AZG75);
        [Max80,Aax80] = max(NumberZ(KM80NaZ:KM85NaZ,Prof));
        Aax80 = round(AltitudeNaZ(Aax80+KM80NaZ),2);
        AZG80 = (sum(NumberZ(KM80NaZ:KM85NaZ,Prof))-Max80)./size(NumberZ(KM80NaZ:KM85NaZ,Prof),1);
        Max80 = round(Max80-AZG80);
        [Max85,Aax85] = max(NumberZ(KM85NaZ:KM90NaZ,Prof));
        Aax85 = round(AltitudeNaZ(Aax85+KM85NaZ),2);
        AZG85 = (sum(NumberZ(KM85NaZ:KM90NaZ,Prof))-Max85)./size(NumberZ(KM85NaZ:KM90NaZ,Prof),1);
        Max85 = round(Max85-AZG85);
        [Max90,Aax90] = max(NumberZ(KM90NaZ:KM95NaZ,Prof));
        Aax90 = round(AltitudeNaZ(Aax90+KM90NaZ),2);
        AZG90 = (sum(NumberZ(KM90NaZ:KM95NaZ,Prof))-Max90)./size(NumberZ(KM90NaZ:KM95NaZ,Prof),1);
        Max90 = round(Max90-AZG90);
        [Max95,Aax95] = max(NumberZ(KM95NaZ:KM100NaZ,Prof));
        Aax95 = round(AltitudeNaZ(Aax95+KM95NaZ),2);
        AZG95 = (sum(NumberZ(KM95NaZ:KM100NaZ,Prof))-Max95)./size(NumberZ(KM95NaZ:KM100NaZ,Prof),1);
        Max95 = round(Max95-AZG95);
        [Max100,Aax100] = max(NumberZ(KM100NaZ:KM105NaZ,Prof));
        Aax100 = round(AltitudeNaZ(Aax100+KM100NaZ),2);
        AZG100 = (sum(NumberZ(KM100NaZ:KM105NaZ,Prof))-Max100)./size(NumberZ(KM100NaZ:KM105NaZ,Prof),1);
        Max100 = round(Max100-AZG100);
        [Max105,Aax105] = max(NumberZ(KM105NaZ:KM110NaZ,Prof));
        Aax105 = round(AltitudeNaZ(Aax105+KM105NaZ),2);
        AZG105 = (sum(NumberZ(KM105NaZ:KM110NaZ,Prof))-Max105)./size(NumberZ(KM105NaZ:KM110NaZ,Prof),1);
        Max105 = round(Max105-AZG105);
        [Max110,Aax110] = max(NumberZ(KM110NaZ:KM115NaZ,Prof));
        Aax110 = round(AltitudeNaZ(Aax110+KM110NaZ),2);
        AZG110 = (sum(NumberZ(KM110NaZ:KM115NaZ,Prof))-Max110)./size(NumberZ(KM110NaZ:KM115NaZ,Prof),1);
        Max110 = round(Max110-AZG110);
        
        Ns = NumberZ(:,Prof);
        plot(Z,Ns,'-','linewidth',1.5)
        grid on;
        ystr = string(Prof)+" | "+string(TimeShuNa(Prof));
        xstr = string(Max75)+" | "+string(Max80)+" | "+string(Max85)+" | "+string(Max90)+" | "+string(Max95)+" | "+string(Max100)+" | "+string(Max105)+" | "+string(Max110);
        tstr = string(Aax75)+" | "+string(Aax80)+" | "+string(Aax85)+" | "+string(Aax90)+" | "+string(Aax95)+" | "+string(Aax100)+" | "+string(Aax105)+" | "+string(Aax110);
        title(tstr)
        xlim([75 115]);
        xlabel(xstr);
        ylabel(ystr);
        
        % 指定保存路径和文件名
        figPath = "D:\KYBF\TimingTrailIMG\Na\"+ReadDate+"\";
        if ~exist(figPath, 'dir')
            mkdir(figPath);
        end
        figV = "Na-Z-"+string(Prof);
        print(fNaZ, '-dpng', '-r300', figPath+figV);
        close all;
    end
    end

end

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);
    
end


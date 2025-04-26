%% �Ƽ�����β���ܶ��뺣�ζ�ȡ

clc; clear all;

% ͨ��list

DateList = {'20231018';'20231019';'20231020';'20231022';'20231023';...
            '20231026';'20231027';'20231028';'20231102';'20231103';...
            '20231104';'20231105';'20231106';'20231107';'20231111';...
            '20231112';'20231113';'20231118';'20231119';'20231121';...
            '20231122';'20231123';...
            '20240114';'20240115';'20240117';'20240118';'20240119';...
            '20240120';'20240121';'20240124';'20240125';'20240128';...
            '20240129';'20240131';...
            '20240311';'20240510';'20240511';...
            '20240512';'20240514';'20240515';'20240518';'20240519';...
            '20240524';'20240525';'20240526';...
            '20240528';'20240529';'20240530';'20240531';'20240601';...
            '20240629';'20240707';};

%% K ��ȡ

DateList = {'20231018';'20231019';'20231020';'20231022';'20231023';...
            '20231026';'20231027';'20231028';'20231102';'20231103';...
            '20231104';'20231105';'20231106';'20231107';'20231111';...
            '20231112';'20231113';'20231118';'20231119';'20231121';...
            '20231122';'20231123';...
            '20240114';'20240115';'20240117';'20240118';'20240119';...
            '20240120';'20240121';'20240124';'20240125';'20240128';...
            '20240129';'20240131';...
            '20240311';'20240510';'20240511';...
            '20240512';'20240514';'20240515';'20240518';'20240519';...
            '20240524';'20240525';'20240526';...
            '20240528';'20240529';'20240530';'20240531';'20240601';...
            '20240629';'20240707';};

for jd = 1:size(DateList,1)
% ��ȡ����������
ReadDate = char(DateList(jd));
RunDate = ReadDate;

Year = str2double(char(RunDate(1:4)));
Month = str2double(char(RunDate(5:6)));
Day = str2double(char(RunDate(7:8)));
DayNum = Month*30-30+Day;

CHN = 'Cap';
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

CapListNum = size(CapList,1);

NumLins = CapListNum;

if NumLins>0

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
    jNum = CapjNum;
    iNum = CapiNum;
    TimeList = cell(1,floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
    end

    TimeXCap = datetime(TimeList, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');
    height_num_origin = PhData1{1,1}(1:4096);

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
    
    TimeShuCap = TimeXCap';

    Noise = mean(PhjSumCap(KM120Cap:KM125Cap,:),1);

    ScatterTable = readtable('RayEffScatter.txt');
    Scatter = table2array(ScatterTable(:,2:3));
    RayScatter = Scatter(3,1);
    EffScatter = Scatter(3,2);

    [TemRay,DenRay] = atmosnrlmsise00(30000,53.3,122.7,Year,DayNum,17);

    Z = AltitudeCap;                       % �߶Ⱦ���
    ZR = 30;                            % �ο��߶�
    SigmaRay = RayScatter;              % ��������ɢ�����
    SigmaEff = EffScatter;              % ��Ч����ɢ�����
    NZ = PhjSumCap;                        % ����������
    NB = Noise;                         % ��������
    NZR = PhjSumCap(KM30Cap,:);                % �ο��߶ȴ�������
    NumberRay = sum(DenRay)*1e-6;       % �ο��߶ȴ�����ģ�����ܶ�
    NumberZ = (((Z.^2).*SigmaRay.*(NZ-NB))./((ZR^2)*SigmaEff*(NZR-NB)))*NumberRay;

    Numjs = size(CapList,1);
    
    for jpf = 1:Numjs
        Prof = CapList(jpf);

        fCap = figure('name','Cap');
        
        [Max75,Aax75] = max(NumberZ(KM75Cap:KM80Cap,Prof));
        Aax75 = round(AltitudeCap(Aax75+KM75Cap),2);
        AVG75 = (sum(NumberZ(KM75Cap:KM80Cap,Prof))-Max75)./size(NumberZ(KM75Cap:KM80Cap,Prof),1);
        Max75 = round(Max75-AVG75);
        [Max80,Aax80] = max(NumberZ(KM80Cap:KM85Cap,Prof));
        Aax80 = round(AltitudeCap(Aax80+KM80Cap),2);
        AVG80 = (sum(NumberZ(KM80Cap:KM85Cap,Prof))-Max80)./size(NumberZ(KM80Cap:KM85Cap,Prof),1);
        Max80 = round(Max80-AVG80);
        [Max85,Aax85] = max(NumberZ(KM85Cap:KM90Cap,Prof));
        Aax85 = round(AltitudeCap(Aax85+KM85Cap),2);
        AVG85 = (sum(NumberZ(KM85Cap:KM90Cap,Prof))-Max85)./size(NumberZ(KM85Cap:KM90Cap,Prof),1);
        Max85 = round(Max85-AVG85);
        [Max90,Aax90] = max(NumberZ(KM90Cap:KM95Cap,Prof));
        Aax90 = round(AltitudeCap(Aax90+KM90Cap),2);
        AVG90 = (sum(NumberZ(KM90Cap:KM95Cap,Prof))-Max90)./size(NumberZ(KM90Cap:KM95Cap,Prof),1);
        Max90 = round(Max90-AVG90);
        [Max95,Aax95] = max(NumberZ(KM95Cap:KM100Cap,Prof));
        Aax95 = round(AltitudeCap(Aax95+KM95Cap),2);
        AVG95 = (sum(NumberZ(KM95Cap:KM100Cap,Prof))-Max95)./size(NumberZ(KM95Cap:KM100Cap,Prof),1);
        Max95 = round(Max95-AVG95);
        [Max100,Aax100] = max(NumberZ(KM100Cap:KM105Cap,Prof));
        Aax100 = round(AltitudeCap(Aax100+KM100Cap),2);
        AVG100 = (sum(NumberZ(KM100Cap:KM105Cap,Prof))-Max100)./size(NumberZ(KM100Cap:KM105Cap,Prof),1);
        Max100 = round(Max100-AVG100);
        [Max105,Aax105] = max(NumberZ(KM105Cap:KM110Cap,Prof));
        Aax105 = round(AltitudeCap(Aax105+KM105Cap),2);
        AVG105 = (sum(NumberZ(KM105Cap:KM110Cap,Prof))-Max105)./size(NumberZ(KM105Cap:KM110Cap,Prof),1);
        Max105 = round(Max105-AVG105);
        [Max110,Aax110] = max(NumberZ(KM110Cap:KM115Cap,Prof));
        Aax110 = round(AltitudeCap(Aax110+KM110Cap),2);
        AVG110 = (sum(NumberZ(KM110Cap:KM115Cap,Prof))-Max110)./size(NumberZ(KM110Cap:KM115Cap,Prof),1);
        Max110 = round(Max110-AVG110);
        
        Ns = NumberZ(:,Prof);
        plot(Z,Ns,'-','linewidth',1.5)
        grid on;
        ystr = string(Prof)+" | "+string(TimeShuCap(Prof));
        xstr = string(Max75)+" | "+string(Max80)+" | "+string(Max85)+" | "+string(Max90)+" | "+string(Max95)+" | "+string(Max100)+" | "+string(Max105)+" | "+string(Max110);
        tstr = string(Aax75)+" | "+string(Aax80)+" | "+string(Aax85)+" | "+string(Aax90)+" | "+string(Aax95)+" | "+string(Aax100)+" | "+string(Aax105)+" | "+string(Aax110);
        title(tstr)
        xlim([75 115]);
        xlabel(xstr);
        ylabel(ystr);
        
        % ָ������·�����ļ���
        figPath = "D:\KYBF\TimingTrailIMG\Cap\"+ReadDate+"\";
        if ~exist(figPath, 'dir')
            mkdir(figPath);
        end
        figV = "Cap-"+string(Prof);
        print(fCap, '-dpng', '-r300', figPath+figV);
        close all;
    end
end
end


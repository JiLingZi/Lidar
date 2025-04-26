%% �Ƽ�����β���ܶ��뺣�ζ�ȡ

clc; clear all;

% ͨ��list

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

%% Ni ��ȡ

DateList = {'20231018';'20231019';'20231020';'20231022';'20231023';...
            '20231026';'20231027';'20231028';'20231102';'20231103';...
            '20231104';'20231105';'20231106';'20231107';'20231111';...
            '20231112';'20231113';'20231118';'20231119';'20231121';...
            '20231122';'20231123';};

for jd = 1:size(DateList,1)
% ��ȡ����������
ReadDate = char(DateList(jd));
RunDate = ReadDate;

Year = str2double(char(RunDate(1:4)));
Month = str2double(char(RunDate(5:6)));
Day = str2double(char(RunDate(7:8)));
DayNum = Month*30-30+Day;

CHN = 'Ni';
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

NiListNum = size(NiList,1);

NumLins = NiListNum;

if NumLins>0

    fldstart = ['F:\RawData\ZWDATA\MOHEnew\Ni\',ReadDate,'\'];
    folder = [fldstart,'CH1\'];
    files = dir(fullfile(folder, '*.txt'));
    num_files = length(files)-1;

    FileNi = num_files;

    DenPh = zeros(4096, num_files);
    for j = 1:num_files
        filename = fullfile(folder, files(j).name);
        fid = fopen(filename);
        PhData = textscan(fid, '%f %f','HeaderLines',20);   % 7��22�ռ�֮ǰ��Ϊ23����Ӧ��Ϊ20
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
    jNum = NijNum;
    iNum = NiiNum;
    TimeList = cell(1,floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
    end

    TimeXNi = datetime(TimeList, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');
    height_num_origin = PhData{1,1}(1:4096);

    PhiSum = zeros(floor(4096/iNum),num_files);
    for i = 1:floor(4096/iNum)
        PhiSum(i,:) = sum(DenPh(iNum*(i-1)+1:iNum*i,:),1);
    end
    PhjSumNi = zeros(floor(4096/iNum),floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        PhjSumNi(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
    end
    AltitudeNi = height_num_origin(1:floor(4096/iNum))*iNum;
    KM30Ni = size(AltitudeNi(AltitudeNi<30),1)+1;
    KM75Ni = size(AltitudeNi(AltitudeNi<75),1)+1;
    KM80Ni = size(AltitudeNi(AltitudeNi<80),1)+1;
    KM85Ni = size(AltitudeNi(AltitudeNi<85),1)+1;
    KM90Ni = size(AltitudeNi(AltitudeNi<90),1)+1;
    KM95Ni = size(AltitudeNi(AltitudeNi<95),1)+1;
    KM100Ni = size(AltitudeNi(AltitudeNi<100),1)+1;
    KM105Ni = size(AltitudeNi(AltitudeNi<105),1)+1;
    KM110Ni = size(AltitudeNi(AltitudeNi<110),1)+1;
    KM115Ni = size(AltitudeNi(AltitudeNi<115),1)+1;
    KM120Ni = size(AltitudeNi(AltitudeNi<120),1)+1;
    KM125Ni = size(AltitudeNi(AltitudeNi<125),1)+1;
    
    TimeShuNi = TimeXNi';

    Noise = mean(PhjSumNi(KM120Ni:KM125Ni,:),1);

    ScatterTable = readtable('RayEffScatter.txt');
    Scatter = table2array(ScatterTable(:,2:3));
    RayScatter = Scatter(6,1);
    EffScatter = Scatter(6,2);

    [TemRay,DenRay] = atmosnrlmsise00(30000,53.3,122.7,Year,DayNum,17);

    Z = AltitudeNi;                       % �߶Ⱦ���
    ZR = 30;                            % �ο��߶�
    SigmaRay = RayScatter;              % ��������ɢ�����
    SigmaEff = EffScatter;              % ��Ч����ɢ�����
    NZ = PhjSumNi;                        % ����������
    NB = Noise;                         % ��������
    NZR = PhjSumNi(KM30Ni,:);                % �ο��߶ȴ�������
    NumberRay = sum(DenRay)*1e-6;       % �ο��߶ȴ�����ģ�����ܶ�
    NumberZ = (((Z.^2).*SigmaRay.*(NZ-NB))./((ZR^2)*SigmaEff*(NZR-NB)))*NumberRay;

    Numjs = size(NiList,1);
    
    for jpf = 1:Numjs
        Prof = NiList(jpf);

        fNi = figure('name','Ni');
        
        [Max75,Aax75] = max(NumberZ(KM75Ni:KM80Ni,Prof));
        Aax75 = round(AltitudeNi(Aax75+KM75Ni),2);
        AVG75 = (sum(NumberZ(KM75Ni:KM80Ni,Prof))-Max75)./size(NumberZ(KM75Ni:KM80Ni,Prof),1);
        Max75 = round(Max75-AVG75);
        [Max80,Aax80] = max(NumberZ(KM80Ni:KM85Ni,Prof));
        Aax80 = round(AltitudeNi(Aax80+KM80Ni),2);
        AVG80 = (sum(NumberZ(KM80Ni:KM85Ni,Prof))-Max80)./size(NumberZ(KM80Ni:KM85Ni,Prof),1);
        Max80 = round(Max80-AVG80);
        [Max85,Aax85] = max(NumberZ(KM85Ni:KM90Ni,Prof));
        Aax85 = round(AltitudeNi(Aax85+KM85Ni),2);
        AVG85 = (sum(NumberZ(KM85Ni:KM90Ni,Prof))-Max85)./size(NumberZ(KM85Ni:KM90Ni,Prof),1);
        Max85 = round(Max85-AVG85);
        [Max90,Aax90] = max(NumberZ(KM90Ni:KM95Ni,Prof));
        Aax90 = round(AltitudeNi(Aax90+KM90Ni),2);
        AVG90 = (sum(NumberZ(KM90Ni:KM95Ni,Prof))-Max90)./size(NumberZ(KM90Ni:KM95Ni,Prof),1);
        Max90 = round(Max90-AVG90);
        [Max95,Aax95] = max(NumberZ(KM95Ni:KM100Ni,Prof));
        Aax95 = round(AltitudeNi(Aax95+KM95Ni),2);
        AVG95 = (sum(NumberZ(KM95Ni:KM100Ni,Prof))-Max95)./size(NumberZ(KM95Ni:KM100Ni,Prof),1);
        Max95 = round(Max95-AVG95);
        [Max100,Aax100] = max(NumberZ(KM100Ni:KM105Ni,Prof));
        Aax100 = round(AltitudeNi(Aax100+KM100Ni),2);
        AVG100 = (sum(NumberZ(KM100Ni:KM105Ni,Prof))-Max100)./size(NumberZ(KM100Ni:KM105Ni,Prof),1);
        Max100 = round(Max100-AVG100);
        [Max105,Aax105] = max(NumberZ(KM105Ni:KM110Ni,Prof));
        Aax105 = round(AltitudeNi(Aax105+KM105Ni),2);
        AVG105 = (sum(NumberZ(KM105Ni:KM110Ni,Prof))-Max105)./size(NumberZ(KM105Ni:KM110Ni,Prof),1);
        Max105 = round(Max105-AVG105);
        [Max110,Aax110] = max(NumberZ(KM110Ni:KM115Ni,Prof));
        Aax110 = round(AltitudeNi(Aax110+KM110Ni),2);
        AVG110 = (sum(NumberZ(KM110Ni:KM115Ni,Prof))-Max110)./size(NumberZ(KM110Ni:KM115Ni,Prof),1);
        Max110 = round(Max110-AVG110);
        
        Ns = NumberZ(:,Prof);
        plot(Z,Ns,'-','linewidth',1.5)
        grid on;
        ystr = string(Prof)+" | "+string(TimeShuNi(Prof));
        xstr = string(Max75)+" | "+string(Max80)+" | "+string(Max85)+" | "+string(Max90)+" | "+string(Max95)+" | "+string(Max100)+" | "+string(Max105)+" | "+string(Max110);
        tstr = string(Aax75)+" | "+string(Aax80)+" | "+string(Aax85)+" | "+string(Aax90)+" | "+string(Aax95)+" | "+string(Aax100)+" | "+string(Aax105)+" | "+string(Aax110);
        title(tstr)
        xlim([75 115]);
        xlabel(xstr);
        ylabel(ystr);
        
        % ָ������·�����ļ���
        figPath = "D:\KYBF\TimingTrailIMG\Ni\"+ReadDate+"\";
        if ~exist(figPath, 'dir')
            mkdir(figPath);
        end
        figV = "Ni-"+string(Prof);
        print(fNi, '-dpng', '-r300', figPath+figV);
        close all;
    end
end
    
    
    


fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);
    
end


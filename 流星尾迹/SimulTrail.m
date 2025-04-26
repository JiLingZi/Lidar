% Date--1
% Time--2
% Rank--3
% NaV---4
% NaM---5
% NaZ---6
% K-----7
% Fe----8
% Ni----9
% Ca----10
% Cap---11
% Ray---12

clc; clear all;

%%

RunDate = '20231019';
CHN = 'Fe';
XlsName = char("F:\TrailIMG\" + CHN + ".xlsx");
SDate = RunDate(5:8);
SheetName = char('S'+string(SDate));
Table = readtable(XlsName,'sheet',SheetName);
NaVList = table2array(Table(:,4));
NaMList = table2array(Table(:,5));
NaZList = table2array(Table(:,6));
KList = table2array(Table(:,7));
FeList = table2array(Table(:,8));
NiList = table2array(Table(:,9));
CaList = table2array(Table(:,10));
CapList = table2array(Table(:,11));
RayList = table2array(Table(:,12));

NajNum = 1; NaiNum = 1;NaDex = 1/1;
KjNum = 1; KiNum = 1;KDex = 1/2;
FejNum = 1; FeiNum = 1;FeDex = 1/2;
NijNum = 1; NiiNum = 1;NiDex = 1/1;
CajNum = 1; CaiNum = 1;CaDex = 1/1;
CapjNum = 1; CapiNum = 1;CapDex = 1/1;
RayjNum = 1; RayiNum = 1;RayDex = KDex;

List = FeList;

PF = List(~isnan(List));
Ttk = 2/1440;
TtkUp = 3/1440;
TtkDp = 3/1440;

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = 'Excel ok '

ReadDate = RunDate;
%% Na
ReadDate = RunDate;
fldstart = ['F:\RawData\ZWDATA\MOHEnew\Na\',ReadDate,'\'];
folder = [fldstart,'Na\'];
files = dir(fullfile(folder, '*.dat'));
num_files = length(files)-1;
num_files = 3*num_files;

FileNa = num_files;

if FileNa>=0
    V_F_0 = zeros(4096, num_files/3);
    V_F_R = zeros(4096, num_files/3);
    V_F_L = zeros(4096, num_files/3);
    N_F_0 = zeros(4096, num_files/3);
    N_F_R = zeros(4096, num_files/3);
    N_F_L = zeros(4096, num_files/3);
    E_F_0 = zeros(4096, num_files/3);
    E_F_R = zeros(4096, num_files/3);
    E_F_L = zeros(4096, num_files/3);
    for j = 1:(num_files/3)
        filename = fullfile(folder, files(j).name);
        Na_table = readtable(filename);
        Na_data = table2array(Na_table);
        V_F_0(:,j) = Na_data(1:4096,2);
        V_F_R(:,j) = Na_data(1:4096,3);
        V_F_L(:,j) = Na_data(1:4096,4);
        N_F_0(:,j) = Na_data(1:4096,5);
        N_F_R(:,j) = Na_data(1:4096,6);
        N_F_L(:,j) = Na_data(1:4096,7);
        E_F_0(:,j) = Na_data(1:4096,8);
        E_F_R(:,j) = Na_data(1:4096,9);
        E_F_L(:,j) = Na_data(1:4096,10);
    end
    DenPhV = zeros(4096, num_files);
    DenPhN = zeros(4096, num_files);
    DenPhE = zeros(4096, num_files);
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
    height_num_origin = Na_data(1:4096,1);
    
    % V
    PhiSum = zeros(floor(4096/iNum),num_files);
    for i = 1:floor(4096/iNum)
        PhiSum(i,:) = sum(DenPhV(iNum*(i-1)+1:iNum*i,:),1);
    end
    PhjSumNaV = zeros(floor(4096/iNum),floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        PhjSumNaV(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
    end
    AltitudeNaV = height_num_origin(1:floor(4096/iNum))*iNum;
    KM30NaV = size(AltitudeNaV(AltitudeNaV<30),1)+1;
    KM75NaV = size(AltitudeNaV(AltitudeNaV<75),1)+1;
    KM80NaV = size(AltitudeNaV(AltitudeNaV<80),1)+1;
    KM110NaV = size(AltitudeNaV(AltitudeNaV<110),1)+1;
    KM115NaV = size(AltitudeNaV(AltitudeNaV<115),1)+1;
    KM120NaV = size(AltitudeNaV(AltitudeNaV<120),1)+1;
    KM125NaV = size(AltitudeNaV(AltitudeNaV<125),1)+1;
    
    % N
    PhiSum = zeros(floor(4096/iNum),num_files);
    for i = 1:floor(4096/iNum)
        PhiSum(i,:) = sum(DenPhN(iNum*(i-1)+1:iNum*i,:),1);
    end
    PhjSumNaM = zeros(floor(4096/iNum),floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        PhjSumNaM(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
    end
    AltitudeNaM = height_num_origin(1:floor(4096/iNum))*iNum*(sqrt(3)/2);
    KM30NaM = size(AltitudeNaM(AltitudeNaM<30),1)+1;
    KM75NaM = size(AltitudeNaM(AltitudeNaM<75),1)+1;
    KM80NaM = size(AltitudeNaM(AltitudeNaM<80),1)+1;
    KM110NaM = size(AltitudeNaM(AltitudeNaM<110),1)+1;
    KM115NaM = size(AltitudeNaM(AltitudeNaM<115),1)+1;
    KM120NaM = size(AltitudeNaM(AltitudeNaM<120),1)+1;
    KM125NaM = size(AltitudeNaM(AltitudeNaM<125),1)+1;
    
    % E
    PhiSum = zeros(floor(4096/iNum),num_files);
    for i = 1:floor(4096/iNum)
        PhiSum(i,:) = sum(DenPhE(iNum*(i-1)+1:iNum*i,:),1);
    end
    PhjSumNaZ = zeros(floor(4096/iNum),floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        PhjSumNaZ(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
    end
    AltitudeNaZ = height_num_origin(1:floor(4096/iNum))*iNum*(sqrt(3)/2);
    KM30NaZ = size(AltitudeNaZ(AltitudeNaZ<30),1)+1;
    KM75NaZ = size(AltitudeNaZ(AltitudeNaZ<75),1)+1;
    KM80NaZ = size(AltitudeNaZ(AltitudeNaZ<80),1)+1;
    KM110NaZ = size(AltitudeNaZ(AltitudeNaZ<110),1)+1;
    KM115NaZ = size(AltitudeNaZ(AltitudeNaZ<115),1)+1;
    KM120NaZ = size(AltitudeNaZ(AltitudeNaZ<120),1)+1;
    KM125NaZ = size(AltitudeNaZ(AltitudeNaZ<125),1)+1;
end

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = 'Na ok '

%% K
fldstart = ['F:\RawData\ZWDATA\MOHEnew\K532\',ReadDate,'\'];
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1;

FileK = num_files;

if FileK>=0

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
    KM110K = size(AltitudeK(AltitudeK<110),1)+1;
    KM115K = size(AltitudeK(AltitudeK<115),1)+1;
    KM120K = size(AltitudeK(AltitudeK<120),1)+1;
    KM125K = size(AltitudeK(AltitudeK<125),1)+1;

end

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = 'K ok '

%% Fe
fldstart = ['F:\RawData\ZWDATA\MOHEnew\Fe\',ReadDate,'\'];
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1;

FileFe = num_files;

if FileFe>=0

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
    jNum = FejNum;
    iNum = FeiNum;
    TimeList = cell(1,floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
    end

    TimeXFe = datetime(TimeList, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');
    height_num_origin = PhData{1,1}(1:4096);

    PhiSum = zeros(floor(4096/iNum),num_files);
    for i = 1:floor(4096/iNum)
        PhiSum(i,:) = sum(DenPh(iNum*(i-1)+1:iNum*i,:),1);
    end
    PhjSumFe = zeros(floor(4096/iNum),floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        PhjSumFe(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
    end
    AltitudeFe = height_num_origin(1:floor(4096/iNum))*iNum;
    KM30Fe = size(AltitudeFe(AltitudeFe<30),1)+1;
    KM75Fe = size(AltitudeFe(AltitudeFe<75),1)+1;
    KM80Fe = size(AltitudeFe(AltitudeFe<80),1)+1;
    KM110Fe = size(AltitudeFe(AltitudeFe<110),1)+1;
    KM115Fe = size(AltitudeFe(AltitudeFe<115),1)+1;
    KM120Fe = size(AltitudeFe(AltitudeFe<120),1)+1;
    KM125Fe = size(AltitudeFe(AltitudeFe<125),1)+1;

end

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = 'Fe ok '


%% Ni
fldstart = ['F:\RawData\ZWDATA\MOHEnew\Ni\',ReadDate,'\'];
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1;

FileNi = num_files;

if FileNi>=0

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
    KM110Ni = size(AltitudeNi(AltitudeNi<110),1)+1;
    KM115Ni = size(AltitudeNi(AltitudeNi<115),1)+1;
    KM120Ni = size(AltitudeNi(AltitudeNi<120),1)+1;
    KM125Ni = size(AltitudeNi(AltitudeNi<125),1)+1;

end

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = 'Ni ok '



%% Ca
fldstart = ['F:\RawData\ZWDATA\MOHEnew\Ca\',ReadDate,'\'];
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1;

FileCa = num_files;

if FileCa>=0

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
    jNum = CajNum;
    iNum = CaiNum;
    TimeList = cell(1,floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
    end

    TimeXCa = datetime(TimeList, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');
    height_num_origin = PhData1{1,1}(1:4096);

    PhiSum = zeros(floor(4096/iNum),num_files);
    for i = 1:floor(4096/iNum)
        PhiSum(i,:) = sum(DenPh(iNum*(i-1)+1:iNum*i,:),1);
    end
    PhjSumCa = zeros(floor(4096/iNum),floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        PhjSumCa(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
    end
    AltitudeCa = height_num_origin(1:floor(4096/iNum))*iNum;
    KM30Ca = size(AltitudeCa(AltitudeCa<30),1)+1;
    KM75Ca = size(AltitudeCa(AltitudeCa<75),1)+1;
    KM80Ca = size(AltitudeCa(AltitudeCa<80),1)+1;
    KM110Ca = size(AltitudeCa(AltitudeCa<110),1)+1;
    KM115Ca = size(AltitudeCa(AltitudeCa<115),1)+1;
    KM120Ca = size(AltitudeCa(AltitudeCa<120),1)+1;
    KM125Ca = size(AltitudeCa(AltitudeCa<125),1)+1;

end

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = 'Ca ok '


%% Cap
fldstart = ['F:\RawData\ZWDATA\MOHEnew\Cap\',ReadDate,'\'];
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1;

FileCap = num_files;

if FileCap>=0

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
    KM110Cap = size(AltitudeCap(AltitudeCap<110),1)+1;
    KM115Cap = size(AltitudeCap(AltitudeCap<115),1)+1;
    KM120Cap = size(AltitudeCap(AltitudeCap<120),1)+1;
    KM125Cap = size(AltitudeCap(AltitudeCap<125),1)+1;

end

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = 'Cap ok '


%% Ray
fldstart = ['F:\RawData\ZWDATA\MOHEnew\K532\',ReadDate,'\'];
folder = [fldstart,'CH2\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1;

FileRay = num_files;

if FileRay>=0

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
    jNum = RayjNum;
    iNum = RayiNum;
    TimeList = cell(1,floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        TimeList(:,j) = TimeValues(:,(j-1)*jNum+1);
    end

    TimeXRay = datetime(TimeList, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');
    height_num_origin = PhData{1,1}(1:4096);

    PhiSum = zeros(floor(4096/iNum),num_files);
    for i = 1:floor(4096/iNum)
        PhiSum(i,:) = sum(DenPh(iNum*(i-1)+1:iNum*i,:),1);
    end
    PhjSumRay = zeros(floor(4096/iNum),floor(num_files/jNum));
    for j = 1:floor(num_files/jNum)
        PhjSumRay(:,j) = sum(PhiSum(:,jNum*(j-1)+1:jNum*j),2);
    end
    AltitudeRay = height_num_origin(1:floor(4096/iNum))*iNum;
    KM30Ray = size(AltitudeRay(AltitudeRay<30),1)+1;
    KM75Ray = size(AltitudeRay(AltitudeRay<75),1)+1;
    KM80Ray = size(AltitudeRay(AltitudeRay<80),1)+1;
    KM110Ray = size(AltitudeRay(AltitudeRay<110),1)+1;
    KM115Ray = size(AltitudeRay(AltitudeRay<115),1)+1;
    KM120Ray = size(AltitudeRay(AltitudeRay<120),1)+1;
    KM125Ray = size(AltitudeRay(AltitudeRay<125),1)+1;

end

fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*800*t);sound(y, fs);

TS = 'Ray ok '



%% Plot Time [10 50 2400 1200]

% FileNa = -1;
% FileNi = -1;
% FileCap = -1;
% FileRay = -1;

TimeSum = TimeXFe;
YL = 80;
YR = 90;

PF = 684

for jt = PF'

start_time = datenum(TimeSum(:,jt))-TtkUp;
end_time = datenum(TimeSum(:,jt))+TtkDp;

fC = figure('name','全队列','position',[10 50 1500 800]);

% NaV

subplot(3,3,1);

if FileNa >= 0

    PhX = PhjSumNaV(KM75NaV:KM115NaV,:);
    PhMax = max(PhX, [], 1);
    PhMin = PhjSumNaV(KM75NaV:KM115NaV,:).*((1/1440)./PhMax).*NaDex;
    NaV = PhMin;
    TimeD = datenum(TimeXNa);
    for jn = 1:size(PhMin,2)
        NaV(:,jn) = datenum(TimeXNa(:,jn)) + PhMin(:,jn);
    end
    plot(NaV,AltitudeNaV(KM75NaV:KM115NaV,:),'linewidth',1.5)
    grid on;
    ylim([YL YR])
    xlim([start_time,end_time]);
    x_ticks = start_time:Ttk:end_time;
    set(gca, 'XTick', x_ticks);
    datetick('x','HH:MM','keeplimits');
    NaVStr  = "NaV " + RunDate;
    title(NaVStr)

    colorX = ['g','b','r'];
    ColorD = repmat(colorX,1,FileNa/3);
    h = findobj(gca, 'Type', 'line');
    for i = 1:numel(h)
        set(h(i), 'Color', ColorD(i));
    end

end

% NaN

subplot(3,3,4);

if FileNa >= 0

    PhX = PhjSumNaM(KM75NaM:KM115NaM,:);
    PhMax = max(PhX, [], 1);
    PhMin = PhjSumNaM(KM75NaM:KM115NaM,:).*((1/1440)./PhMax).*NaDex;
    NaM = PhMin;
    TimeD = datenum(TimeXNa);
    for jn = 1:size(PhMin,2)
        NaM(:,jn) = datenum(TimeXNa(:,jn)) + PhMin(:,jn);
    end
    plot(NaM,AltitudeNaM(KM75NaM:KM115NaM,:),'linewidth',1.5)
    grid on;
    ylim([YL YR])
    xlim([start_time,end_time]);
    x_ticks = start_time:Ttk:end_time;
    set(gca, 'XTick', x_ticks);
    datetick('x','HH:MM','keeplimits');
    NaMStr  = "NaM " + RunDate;
    title(NaMStr)

    colorX = ['g','b','r'];
    ColorD = repmat(colorX,1,FileNa/3);
    h = findobj(gca, 'Type', 'line');
    for i = 1:numel(h)
        set(h(i), 'Color', ColorD(i));
    end

end

% NaE

subplot(3,3,7);

if FileNa >= 0

    PhX = PhjSumNaZ(KM75NaZ:KM115NaZ,:);
    PhMax = max(PhX, [], 1);
    PhMin = PhjSumNaZ(KM75NaZ:KM115NaZ,:).*((1/1440)./PhMax).*NaDex;
    NaZ = PhMin;
    TimeD = datenum(TimeXNa);
    for jn = 1:size(PhMin,2)
        NaZ(:,jn) = datenum(TimeXNa(:,jn)) + PhMin(:,jn);
    end
    plot(NaZ,AltitudeNaZ(KM75NaZ:KM115NaZ,:),'linewidth',1.5)
    grid on;
    ylim([YL YR])
    xlim([start_time,end_time]);
    x_ticks = start_time:Ttk:end_time;
    set(gca, 'XTick', x_ticks);
    datetick('x','HH:MM','keeplimits');
    NaZStr  = "NaZ " + RunDate;
    title(NaZStr)

    colorX = ['g','b','r'];
    ColorD = repmat(colorX,1,FileNa/3);
    h = findobj(gca, 'Type', 'line');
    for i = 1:numel(h)
        set(h(i), 'Color', ColorD(i));
    end

end

% K

subplot(3,3,2);

if FileK >= 0

    PhX = PhjSumK(KM75K:KM115K,:);
    PhMax = max(PhX, [], 1);
    PhMin = PhjSumK(KM75K:KM115K,:).*((1/1440)./PhMax).*KDex;
    K = PhMin;
    TimeD = datenum(TimeXK);
    for jn = 1:size(PhMin,2)
        K(:,jn) = datenum(TimeXK(:,jn)) + PhMin(:,jn);
    end
    plot(K,AltitudeK(KM75K:KM115K,:),'linewidth',1.5)
    grid on;
    ylim([YL YR])
    xlim([start_time,end_time]);
    x_ticks = start_time:Ttk:end_time;
    set(gca, 'XTick', x_ticks);
    datetick('x','HH:MM','keeplimits');
    KStr  = "K " + RunDate;
    title(KStr)
    
end

subplot(3,3,5);

if FileFe >= 0

    PhX = PhjSumFe(KM75Fe:KM115Fe,:);
    PhMax = max(PhX, [], 1);
    PhMin = PhjSumFe(KM75Fe:KM115Fe,:).*((1/1440)./PhMax).*FeDex;
    Fe = PhMin;
    TimeD = datenum(TimeXFe);
    for jn = 1:size(PhMin,2)
        Fe(:,jn) = datenum(TimeXFe(:,jn)) + PhMin(:,jn);
    end
    plot(Fe,AltitudeFe(KM75Fe:KM115Fe,:),'linewidth',1.5)
    grid on;
    ylim([YL YR])
    xlim([start_time,end_time]);
    x_ticks = start_time:Ttk:end_time;
    set(gca, 'XTick', x_ticks);
    datetick('x','HH:MM','keeplimits');
    FeStr  = "Fe " + RunDate;
    title(FeStr)
    
end

subplot(3,3,8);

if FileNi >= 0

    PhX = PhjSumNi(KM75Ni:KM115Ni,:);
    PhMax = max(PhX, [], 1);
    PhMin = PhjSumNi(KM75Ni:KM115Ni,:).*((1/1440)./PhMax).*NiDex;
    Ni = PhMin;
    TimeD = datenum(TimeXNi);
    for jn = 1:size(PhMin,2)
        Ni(:,jn) = datenum(TimeXNi(:,jn)) + PhMin(:,jn);
    end
    plot(Ni,AltitudeNi(KM75Ni:KM115Ni,:),'linewidth',1.5)
    grid on;
    ylim([YL YR])
    xlim([start_time,end_time]);
    x_ticks = start_time:Ttk:end_time;
    set(gca, 'XTick', x_ticks);
    datetick('x','HH:MM','keeplimits');
    NiStr  = "Ni " + RunDate;
    title(NiStr)
    
end

subplot(3,3,3);

if FileCa >= 0

    PhX = PhjSumCa(KM75Ca:KM115Ca,:);
    PhMax = max(PhX, [], 1);
    PhMin = PhjSumCa(KM75Ca:KM115Ca,:).*((1/1440)./PhMax).*CaDex;
    Ca = PhMin;
    TimeD = datenum(TimeXCa);
    for jn = 1:size(PhMin,2)
        Ca(:,jn) = datenum(TimeXCa(:,jn)) + PhMin(:,jn);
    end
    plot(Ca,AltitudeCa(KM75Ca:KM115Ca,:),'linewidth',1.5)
    grid on;
    ylim([YL YR])
    xlim([start_time,end_time]);
    x_ticks = start_time:Ttk:end_time;
    set(gca, 'XTick', x_ticks);
    datetick('x','HH:MM','keeplimits');
    CaStr  = "Ca " + RunDate;
    title(CaStr)
    
end

subplot(3,3,6);

if FileCap >= 0

    PhX = PhjSumCap(KM75Cap:KM115Cap,:);
    PhMax = max(PhX, [], 1);
    PhMin = PhjSumCap(KM75Cap:KM115Cap,:).*((1/1440)./PhMax).*CapDex;
    Cap = PhMin;
    TimeD = datenum(TimeXCap);
    for jn = 1:size(PhMin,2)
        Cap(:,jn) = datenum(TimeXCap(:,jn)) + PhMin(:,jn);
    end
    plot(Cap,AltitudeCap(KM75Cap:KM115Cap,:),'linewidth',1.5)
    grid on;
    ylim([YL YR])
    xlim([start_time,end_time]);
    x_ticks = start_time:Ttk:end_time;
    set(gca, 'XTick', x_ticks);
    datetick('x','HH:MM','keeplimits');
    CapStr  = "Cap " + RunDate;
    title(CapStr)
    
end

subplot(3,3,9);

if FileRay >= 0

    PhX = PhjSumRay(KM75Ray:KM115Ray,:);
    PhMax = max(PhX, [], 1);
    PhMin = PhjSumRay(KM75Ray:KM115Ray,:).*((1/1440)./PhMax).*RayDex;
    Ray = PhMin;
    TimeD = datenum(TimeXRay);
    for jn = 1:size(PhMin,2)
        Ray(:,jn) = datenum(TimeXRay(:,jn)) + PhMin(:,jn);
    end
    plot(Ray,AltitudeRay(KM75Ray:KM115Ray,:),'linewidth',1.5)
    grid on;
    ylim([YL YR])
    xlim([start_time,end_time]);
    x_ticks = start_time:Ttk:end_time;
    set(gca, 'XTick', x_ticks);
    datetick('x','HH:MM','keeplimits');
    RayStr  = "Ray " + RunDate;
    title(RayStr)
    
end

%     % 指定保存路径和文件名
%     figPath = "F:\TrailIMG\SumTrailFind\"+ReadDate+"\";
%     if ~exist(figPath, 'dir')
%         mkdir(figPath);
%     end
%     figV = datestr(TimeSum(:,jt),'yyyymmdd-HH-MM-SS')+"-Fe-"+string(jt);
%     print(fC, '-dpng', '-r300', figPath+figV);
%     close all;

end

TS = 'Figures ok '


fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*1000*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.3;y = sin(2*pi*1200*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.1;y = sin(2*pi*900*t);sound(y, fs);
pause(0.3);
fs = 10000;t = 0:1/fs:0.5;y = sin(2*pi*800*t);sound(y, fs);



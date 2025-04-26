% CaIon密度
% 读取光子数矩阵,钙原子、钙离子有4个通道
ReadDate = '20231106';
fldstart = ['F:\RawData\ZWDATA\MOHEnew\Cap\',ReadDate,'\'];
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1;% 为了统一文件数，减一
% 减一是因为那个破机器动不动少一个文件，导致不晓得哪个通道文件数不同
% 格式不要再改了，探测高度上限也不要再改了，否则实现自动化反演困难重重

bin = 26432;
% bin = 16384;
% bin = 16362;
eveHang = 23;
tt = 9;
ttstr = 'yyyy-MM-dd''T''HH:mm:ss.SSS';

DenPh1 = zeros(bin, num_files);
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));

Lines = cell(26452,1);
for j = 1:num_files
    filename = fullfile(folder,files(j).name);
    fid = fopen(filename,'r');
    for i = 1:26452
        Lines{i,:} = fgetl(fid);
    end
    fclose(fid);
    datacell = Lines(21:end,:);
    % 初始化一个空的数组来存储数据
    data = zeros(size(datacell,1), 2);

    % 对于 cell 数组中的每个元素
    for i = 1:size(datacell,1)
        % 使用 str2num 函数将字符串转换为数字数组
        data(i, :) = str2num(datacell{i});
    end
    DenPh1(:,j) = data(:,2);
end
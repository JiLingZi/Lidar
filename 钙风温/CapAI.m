% CaIon�ܶ�
% ��ȡ����������,��ԭ�ӡ���������4��ͨ��
ReadDate = '20231106';
fldstart = ['F:\RawData\ZWDATA\MOHEnew\Cap\',ReadDate,'\'];
folder = [fldstart,'CH1\'];
files = dir(fullfile(folder, '*.txt'));
num_files = length(files)-1;% Ϊ��ͳһ�ļ�������һ
% ��һ����Ϊ�Ǹ��ƻ�����������һ���ļ������²������ĸ�ͨ���ļ�����ͬ
% ��ʽ��Ҫ�ٸ��ˣ�̽��߶�����Ҳ��Ҫ�ٸ��ˣ�����ʵ���Զ���������������

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
    % ��ʼ��һ���յ��������洢����
    data = zeros(size(datacell,1), 2);

    % ���� cell �����е�ÿ��Ԫ��
    for i = 1:size(datacell,1)
        % ʹ�� str2num �������ַ���ת��Ϊ��������
        data(i, :) = str2num(datacell{i});
    end
    DenPh1(:,j) = data(:,2);
end
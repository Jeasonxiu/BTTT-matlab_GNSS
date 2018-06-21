function [yrs, dvs, dvm, res, YMD] = rmoutl(file_gds)
%
% ���߿� �ۼ�
%

%% �Լ��ۼ� �� ����� �׽�Ʈ
% file_gds = 'SUWN.rad';

%% ���� ���� c1(yrscale), c2(deviation), c3(sigma), c4(SITE), c5(LAT/LON/RAD), c6(YMD)
fid_gds = fopen(file_gds,'r');
cell_llh = textscan(fid_gds, '%f %f %f %s %s %s'); %: c1(), c2(), c3(),
fclose(fid_gds);

%% cell�� vector�� ��ȯ
yrs = cell_llh{1,1};
dvs = cell_llh{1,2};
sig = cell_llh{1,3};
ymd = cell_llh{1,6};    % ��� ����
%% c3(sig)�� ���� �� �̻��̸� ����
ind_sig = find(sig > 1.0); %: ������ֵ���
yrs(ind_sig) = [];
dvs(ind_sig) = [];
sig(ind_sig) = [];


%% sig �� ���� �� �̻��� ��¥�� �����Ұ�� �ʱ� �� ����
if length(ind_sig)>0
    YMD=ymd(ind_sig);       % �̻��� �߻���
else
    YMD = {};       % �̻��� �߻� ���� �� ����
end
%% ����ȸ�ͺм� ��� 3*sigma �̻��̸� ���� - �ݺ� ���
fprintf('\n %s --- \n', file_gds);
fprintf('sig 1.0 �̻�: %3d ��\n', length(ind_sig))
kIter = 0;

P = polyfit(yrs,dvs,1);
dvm = P(1)*yrs + P(2);

res = dvs - dvm;
std_res = std(res);
ind_out = find(abs(res) > 3*std_res); %: ������ֵ���
%% 0�� 3*sigma �̻��� �� ����
YMD=[YMD;ymd(ind_out)];
fprintf('        0 ��: %3d ��\n', length(ind_out))

while ~isempty(ind_out)
    kIter = kIter + 1;
    yrs(ind_out) = [];
    dvs(ind_out) = [];
    sig(ind_out) = [];
    
    P = polyfit(yrs,dvs,1);
    dvm = P(1)*yrs + P(2);
    res = dvs - dvm;
    std_res = std(res);
    ind_out = find(abs(res) > 3*std_res); %: ������ֵ���
    fprintf('       %2d ��: %3d ��\n', kIter, length(ind_out))
    YMD = [YMD;ymd(ind_out)];      % �̻��� �߻� ���� �� ����
end
YMD = lower(unique(YMD));          % �� �ҹ��� ��ȯ

% subplot(2,1,1)
% plot(yrs, dvs, 'or', yrs, dvm, '-b')
% subplot(2,1,2)
% plot(yrs, res, 'x:')




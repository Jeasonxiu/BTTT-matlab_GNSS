function [yrs, dvs, dvm, res, YMD] = rmoutl(file_gds)
%
% 나중에 작성
%

%% 함수작성 전 입출력 테스트
% file_gds = 'SUWN.rad';

%% 파일 지정 c1(yrscale), c2(deviation), c3(sigma), c4(SITE), c5(LAT/LON/RAD), c6(YMD)
fid_gds = fopen(file_gds,'r');
cell_llh = textscan(fid_gds, '%f %f %f %s %s %s'); %: c1(), c2(), c3(),
fclose(fid_gds);

%% cell을 vector로 변환
yrs = cell_llh{1,1};
dvs = cell_llh{1,2};
sig = cell_llh{1,3};
ymd = cell_llh{1,6};    % 계산 일자
%% c3(sig)가 일정 값 이상이면 제거
ind_sig = find(sig > 1.0); %: 출력해주도록
yrs(ind_sig) = [];
dvs(ind_sig) = [];
sig(ind_sig) = [];


%% sig 가 일정 값 이상인 날짜가 존재할경우 초기 셀 생성
if length(ind_sig)>0
    YMD=ymd(ind_sig);       % 이상점 발생일
else
    YMD = {};       % 이상점 발생 일자 빈셀 생성
end
%% 선형회귀분석 결과 3*sigma 이상이면 제거 - 반복 계산
fprintf('\n %s --- \n', file_gds);
fprintf('sig 1.0 이상: %3d 개\n', length(ind_sig))
kIter = 0;

P = polyfit(yrs,dvs,1);
dvm = P(1)*yrs + P(2);

res = dvs - dvm;
std_res = std(res);
ind_out = find(abs(res) > 3*std_res); %: 출력해주도록
%% 0차 3*sigma 이상점 셀 생성
YMD=[YMD;ymd(ind_out)];
fprintf('        0 차: %3d 개\n', length(ind_out))

while ~isempty(ind_out)
    kIter = kIter + 1;
    yrs(ind_out) = [];
    dvs(ind_out) = [];
    sig(ind_out) = [];
    
    P = polyfit(yrs,dvs,1);
    dvm = P(1)*yrs + P(2);
    res = dvs - dvm;
    std_res = std(res);
    ind_out = find(abs(res) > 3*std_res); %: 출력해주도록
    fprintf('       %2d 차: %3d 개\n', kIter, length(ind_out))
    YMD = [YMD;ymd(ind_out)];      % 이상점 발생 일자 셀 생성
end
YMD = lower(unique(YMD));          % 월 소문자 변환

% subplot(2,1,1)
% plot(yrs, dvs, 'or', yrs, dvm, '-b')
% subplot(2,1,2)
% plot(yrs, res, 'x:')




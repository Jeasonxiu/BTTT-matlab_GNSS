function [arrQM, FinalPRNs, FinalTTs] = ReadQM(qm_file)
%
%function [arrQM, FinalPRNs, FinalTTs] = ReadQM(qm_file)
%
% DO: Get the unique set of PRNs and TT(Time Tag)s in the given QM file
%
% <input>   qm_file: Quick Measurements file
%           
% <output>  arrQM: QM array c1(gs), c2(PRN), c3(Obs. Type), c4(Measurement)
%           FinalPRNs: PRN array from QM file
%           FinalTTs:  Time-Tag array from QM files
%   
% Copyright: Coded by Kwan-Dong Park, January 12, 2014
%               - Modified the original version dated 2006: use "unique"
%               instead of 'for' loops

%% 함수 만들기 전 입출력 테스트
% qm_file = 'QM07195_is23';
%%
fid_qm = fopen(qm_file,'r');
arrQM = fscanf(fid_qm,'%f %d %d %f',[4 Inf]);
arrQM = arrQM';

TTs = arrQM(:,1);
PRNs = arrQM(:,2);

FinalTTs = unique(TTs);
FinalPRNs = unique(PRNs);
%% 깔끔한 뒷 정리 - 파일 핸들 닫기
fclose(fid_qm);

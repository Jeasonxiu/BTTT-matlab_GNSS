function [arrQM, FinalPRNs, FinalTTs] = ReadQM(qm_file)
%
%function [arrQM, FinalPRNs, FinalTTs] = ReadQM(qm_file)
%
% DO: Get the unique set of PRNs and TT(Time Tag)s in the given QM file
%
% Copyright: Coded by Kwan-Dong Park, July 23, 2006
%
%* Variables
%  TTs:  Array, Time Tags from the original QM file
%  PRNs: Array, PRN number from the original QM file
%  MaxPRNs: Maximum number of PRNS, currently set to 38
%  MaxTTs:  Maximum number of Time Tages, currently set to 86400
%                1 day = 86400 seconds
%  FinalPRNs: Array, Sorted PRNs from the QM file
%  FinalTTs:  Array, Sorted TTs  from the QM file

fid_qm = fopen(qm_file,'r');
arrQM = fscanf(fid_qm,'%f %d %d %f',[4 Inf]);
arrQM = arrQM';

TTs = arrQM(:,1);
PRNs = arrQM(:,2);

qmLen = length(arrQM);

MaxPRNs = 38;
UniqPRNs = zeros(MaxPRNs,1);
MaxTTs = 86400;
UniqTTs = zeros(MaxTTs,1);

NoUniqTT = 1;   %* Number of unique TT
NoUniqPRN = 1;  %* Number of unique PRN
CurrentUniqTT = TTs(1); UniqTTs(1) = TTs(1);
CurrentUniqPRN = PRNs(1); UniqPRNs(1) = PRNs(1);

for i = 2:qmLen
    if TTs(i) ~= CurrentUniqTT
        NoUniqTT = NoUniqTT + 1;
        CurrentUniqTT = TTs(i);
        UniqTTs(NoUniqTT) = TTs(i);
    end
    found = find(UniqPRNs == PRNs(i));
    if length(found) == 0
        NoUniqPRN = NoUniqPRN + 1;
        UniqPRNs(NoUniqPRN) = PRNs(i);
    end
end

FinalPRNs = sort(UniqPRNs(1:NoUniqPRN));
FinalTTs = UniqTTs(1:NoUniqTT);
%--
fclose(fid_qm);

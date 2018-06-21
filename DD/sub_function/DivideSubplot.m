function [nPage, subp] = DivideSubplot(nFigs_page, subp)
%
%function [nPage, subp] = DivideSubplot(nFigs_page, subp)
%
% DO: Decide the page number and location of subplot
%
% <input>   nFigs_page: Number of (subplot) figures per page
%           subp: The number or order of subplots
%
% <output>  nPage: The number or order of 'whole' figure page
%           subp: The number or order of subplot in the 'nPage'
%
% Copyright: 13-12-03, Kwan-Dong Park @LDEO
%
% �߰� ����: subplot�� �� �������� ���� �� �׸����� �� �� �ش� subplot�� �� ��° ������, �� ��° ��ġ��
% �������� ������ �ִ� �Լ���
%
%% �Լ� �ۼ� �� ����� �׽�Ʈ
% nFigs_page = 6;
% for subp = 1:14
%     nPage = floor(subp/nFigs_page) + 1;
%     subp = mod(subp,nFigs_page);
%     if subp == 0
%         nPage = nPage - 1;
%         subp = subp + nFigs_page;
%     end
%     fprintf('%4d %4d\n', nPage, subp);
% end
%%
nPage = floor(subp/nFigs_page) + 1;
subp = mod(subp, nFigs_page);
if subp == 0
    nPage = nPage - 1;
    subp = subp + nFigs_page;
end
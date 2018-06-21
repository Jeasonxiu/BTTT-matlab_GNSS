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
% 추가 설명: subplot을 한 페이지에 여러 개 그리려고 할 때 해당 subplot이 몇 번째 페이지, 몇 번째 위치에
% 들어가는지를 결정해 주는 함수임
%
%% 함수 작성 전 입출력 테스트
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
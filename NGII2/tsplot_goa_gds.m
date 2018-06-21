function [] = tsplot_goa_gds(site)
%function [] = tsplot_goa_gds(site1)
%
% figure(1): Plot time series - lat.site1, lon.site1, and rad.site1
% figure(2): Plot de-trended time series
%

site_cap = upper(site);
file_gds1 = strcat(site_cap,'.lat');
file_gds2 = strcat(site_cap,'.lon');
file_gds3 = strcat(site_cap,'.rad');

%% rmoutl �Լ� ��� ���� �߰��߽��ϴ�. ymd : yearmonthday
%% rmoutl_stat �Լ��� 'mm' ���� ��ºκ��� �����Ͽ����ϴ�.
%% North - latitude
[yrs, dvs, dvm, res, ymd1] = rmoutl(file_gds1); tS = yrs(1); tE = yrs(end);
rmoutl_stat(yrs, dvs, res);
figure(1); subplot(3,1,1)
plot(yrs,dvs,'.r:', yrs,dvm, '-k')
ylabel('\Delta N (cm)'); xlim([tS tE]); title(site_cap); 
figure(2); subplot(3,1,1);
plot(yrs,res,'or');
ylabel('\Delta N (cm)'); xlim([tS tE]); title(site_cap); 

%% East - longitude
[yrs, dvs, dvm, res, ymd2] = rmoutl(file_gds2);
rmoutl_stat(yrs, dvs, res);
figure(1); subplot(3,1,2)
plot(yrs,dvs,'.r:', yrs,dvm, '-k')
ylabel('\Delta E (cm)'); xlim([tS tE]);
figure(2); subplot(3,1,2);
plot(yrs,res,'or'); xlim([tS tE]);
ylabel('\Delta E (cm)')

%% Up - radial
[yrs, dvs, dvm, res, ymd3] = rmoutl(file_gds3);
rmoutl_stat(yrs, dvs, res);
figure(1); subplot(3,1,3)
plot(yrs,dvs,'.r:', yrs,dvm, '-k')
ylabel('\Delta V (cm)'); xlabel('Year');  xlim([tS tE]);
figure(2); subplot(3,1,3);
plot(yrs,res,'or'); 
ylabel('\Delta V (cm)'); xlabel('Year'); xlim([tS tE]);

%% �̻��� �߻� ��¥, ����Ʈ ���� ���� ����
%% ���� ���ϸ� 'RMSTAS_'+site��, ex> RMSTAS_chcn
ymd = [ymd1;ymd2;ymd3];
YMD = unique(ymd(:,1));     % �ߺ�������
for i = 1:length(YMD)
    YMD{i,2} = site;          % YMD���� ����Ʈ�� �߰�(2��), �Էµ� ����Ʈ�� �״��
%     YMD{i,2} = site_cap;     % YMD���� ����Ʈ�� �߰�(2��), �Էµ� ����Ʈ�� �빮�ڷ�
end
YMD = YMD';             % ���ϻ����� �����ϱ� ���� �� ��ġ
YMD_file = fopen(strcat('RMSTAS_',site),'w');             % �Էµ� ����Ʈ�� �״�� ���ϻ���
% YMD_file = fopen(strcat('RMSTAS_',site_cap),'w');         % �Է¿� ����Ʈ�� �빮�ڷ� ���ϻ���
fprintf(YMD_file,'%s %s \n',YMD{:});
fclose(YMD_file);

% fprintf(YMD_file,'%s \n',YMD{:});
% fclose(YMD_file)
% figure(1); subplot(3,1,3)
% plot(yr,ts,'.k:',yr,b*yr+a,'-r');
% title3=strcat(site,'.....',num2str(b),'\pm',sigVstr,'(',repeatStr,')');
% title(title3);
% xlabel('Year')
% ylabel('\Delta V (mm)')
% 
% figure(2); subplot(3,1,3)
% plot(yr, ts - ( b * yr + a), '.k-')
% axis([xmin xmax -30 30])
% ylabel('\Delta V (mm)');
% xlabel('Year')

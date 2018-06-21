function Skymap(cut_ele)

% -- skyplot �� �׸��� �Լ�
%
%     Usage:       skyplot(Az, El, cut_ele) 
%
%          eg:       skyplot(120, 30, 10)
% 
% -- ������, ����, �Ӱ������ degree�� �Է��Ѵ�(������ ��ġ �Ǵ� �迭����).
% -- Hye-In, Kim

hold on;

circlea = 0:pi/30:2*pi;
cx  = cos(circlea);
cy  = sin(circlea);

% 30��, 60��, 90���� ��Ÿ���� ���� �Ǽ����� �׸���(������ ��Ÿ���� �κ�).
for i= [30 60 90]
    plot(cx*i,cy*i,'-','color','k','linewidth',1);
end

% 15��, 45��, 75���� ��Ÿ���� ���� �������� �׸���(������ ��Ÿ���� �κ�).
for i=[15 45 75]
    plot(cx*i,cy*i,':','color','k','linewidth',1);
end

% �Ӱ������ �Ķ��� �Ǽ����� ��Ÿ����. 
plot(cx*(cut_ele-90),cy*(cut_ele-90),'-','color','b','linewidth',1);

lenmax = 90;

% ���� 12������� ������(�������� ��Ÿ���� �κ�).
circleTick = (1:6)*pi/6;
cosct = cos(circleTick); 
sinct = sin(circleTick);
cax = [-cosct; cosct];
say = [-sinct; sinct];
plot(lenmax*cax,lenmax*say,'-','color','k','linewidth',1);

%filename = strcat('[ "', filename,'" Invisibility Sky Plot ]');
% title(site, 'position',[-60 95]);

% ������ ��Ÿ���� �ؽ�Ʈ�� �����Ѵ�(30, 60, 90, 120�� ....).
rlen = 1.06*lenmax;
for i = 1:length(circleTick) 
    ticm1 = int2str(i*30);
    ticm2 = int2str(180+i*30);
    
    % ������ 360���� ��Ÿ���� �κ��� ������ ��Ÿ���� 'N'���� ǥ���Ѵ�.
    if ticm2 == '360'
        ticm2 ='N';
    end
    
    text( rlen*sinct(i), rlen*cosct(i),ticm1,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','horizontalalignment','center');     
    text(-rlen*sinct(i),-rlen*cosct(i),ticm2,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','horizontalalignment','center');
end

axis('equal');
axis('off');

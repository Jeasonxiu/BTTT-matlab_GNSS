function fig2pdf
%
%   �����ִ� ��� figureâ�� pdf�� ����
%
% Coded by Joonseong, 18,Jul,2018

hf = get(0,'children');
hf = flipud(hf);
for ix=1:length(hf)
    set(hf(ix),'Units','centimeters');
    pos = get(hf(ix),'Position');
    set(hf(ix),'PaperPositionMode','Auto',...
        'PaperUnits','centimeters','PaperSize',[pos(3), pos(4)]);
    figurename = ['figure' num2str(ix)];
    print(hf(ix), figurename, '-dpdf', '-painters');
end
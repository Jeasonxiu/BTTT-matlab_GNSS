%    Draw Circle
close all

circlea = 0:pi/30:2*pi;
cx = cos(circlea);
cy = sin(circlea);
for j = 1:3
    figure(j)
    plot([0,0],[6,0],'r')
    title('Hotizontal Error','Fontsize',20)
    hold on
    grid on
    plot([0,0],[-6,0],'r-')
    plot([-6,0],[0,0],'r-')
    plot([6,0],[0,0],'r-')
    switch j
        case 1
            axis([-2.4 2.4 -2.4 2.4])
            error1=[0.4 2];
            error2=[0.8 1.2 1.6];
        case 2
            axis([-6 6 -6 6])
            error1=[1 5];
            error2=[2 3 4];
        case 3
            axis([-0.06 0.06 -0.06 0.06])
            error1=[0.01 0.05];
            error2=[0.02 0.03 0.04];
    end
    
    for i= error1
        plot(cx*i, cy*i, '-', 'color', 'k', 'linewidth', 1);
        if error1(1,1) < 0.4
            text(0.002,i+0.005,strcat(num2str(i*100),'cm'),'Fontsize',15)
            xlabel('East (cm)')
            ylabel('North (cm)')
        elseif error1(1,1) < 1
            text(0.1,i+0.2,strcat(num2str(i),'m'),'Fontsize',15)
            xlabel('East (meter)')
            ylabel('North (meter)')
        else
            text(0.2,i+0.4,strcat(num2str(i),'m'),'Fontsize',15)
            xlabel('East (meter)')
            ylabel('North (meter)')
        end
        axis square
    end
    
    for i=error2
        plot(cx*i, cy*i, ':', 'color', 'k', 'linewidth', 1);
        if error1(1,1) < 0.4
            text(0.002,i+0.005,strcat(num2str(i*100),'cm'),'Fontsize',12)
            xlabel('East (cm)')
            ylabel('North (cm)')
        elseif error1(1,1) < 1
            text(0.1,i+0.2,strcat(num2str(i),'m'),'Fontsize',12)
            xlabel('East (meter)')
            ylabel('North (meter)')
        else
            text(0.2,i+0.4,strcat(num2str(i),'m'),'Fontsize',12)
            xlabel('East (meter)')
            ylabel('North (meter)')
        end
        axis square
    end
end

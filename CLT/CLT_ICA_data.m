clear all
close all
load('CLT_ICA_data.mat');

TruePos = [-3012527.5380 4078281.3285 3856555.4327];
gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2);
for j=1:3
    for k=1:4
        if j==1
            temp1 = (['J_1Hz_L',num2str(k)]);
            temp2 = (['J_1Hz_R',num2str(k)]);
            temp3 = (['J_1Hz_',num2str(k),'_average']);
            temp4 = (['FinalTTs=intersect(',temp1,'(:,1),',temp2,'(:,1));']);
            eval(temp4);
        elseif j==2
            temp1 = (['J_5Hz_L',num2str(k)]);
            temp2 = (['J_5Hz_R',num2str(k)]);
            temp3 = (['J_5Hz_',num2str(k),'_average']);
            temp4 = (['FinalTTs=intersect(',temp1,'(:,1),',temp2,'(:,1));']);
            eval(temp4);
        elseif j==3
            temp1 = (['J_10Hz_L',num2str(k)]);
            temp2 = (['J_10Hz_R',num2str(k)]);
            temp3 = (['J_10Hz_',num2str(k),'_average']);
            temp4 = (['FinalTTs=intersect(',temp1,'(:,1),',temp2,'(:,1));']);
            eval(temp4);
        end
        for i=1:length(FinalTTs)
            gs = FinalTTs(i,1);
            temp5 =([temp3,'(i,1)= gs;']);
            temp6 =([temp3,'(i,2:4)=','(',temp1,'(find(',temp1,'(:,1)==gs),2:4)+',...
                temp2,'(find(',temp2,'(:,1)==gs),2:4))/2;']);
            eval(temp5);
            eval(temp6);
            temp7 =(['GGA=',temp3,'(i,2:4);']);
            eval(temp7);
            dXYZ = [GGA(1), GGA(2), GGA(3)] - TruePos;
            dNEV = xyz2topo(dXYZ, TrueLat, TrueLon);
            dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
            dNE = sqrt(dN^2 + dE^2);        %rmsH = myRMS(dNE);
            d2D(i,1) = dNE;
            d3 = sqrt(dN.^2 + dE.^2 + dV.^2); %rms3 = myRMS(d3);
            d3D(i,1) = d3;
            temp8 = ([temp3,'_result(i,:) = [gs, dN, dE, dV, dNE, d3];']);
            eval(temp8)
        end
    end
end

for i=1:12
    
    open airport_dNE.fig
    hold on
    plot(gd(:,2),gd(:,1),'o:')
    axis([-20.243590755988617 14.756409244011394 -3.0657140562647016 31.934285943735308])
    xlabel('\delta E')
    ylabel('\delta N')
    title('170112 인천공항 1Hz, 5Hz, 10Hz')
    axis square
    figure(i)

    switch i
        case 1
            plot(J_1Hz_1_average_result(:,3),J_1Hz_1_average_result(:,2),'bo')
            hold on; grid on;
        case 2
            plot(J_1Hz_2_average_result(:,3),J_1Hz_2_average_result(:,2),'bo')
            hold on; grid on;
        case 3
            plot(J_1Hz_3_average_result(:,3),J_1Hz_3_average_result(:,2),'bo')
            hold on; grid on;
        case 4
            plot(J_1Hz_4_average_result(:,3),J_1Hz_4_average_result(:,2),'bo')
            hold on; grid on;
        case 5
            plot(J_5Hz_1_average_result(:,3),J_5Hz_1_average_result(:,2),'bo')
            hold on; grid on;
        case 6
            plot(J_5Hz_2_average_result(:,3),J_5Hz_2_average_result(:,2),'bo')
            hold on; grid on;
        case 7
            plot(J_5Hz_3_average_result(:,3),J_5Hz_3_average_result(:,2),'bo')
            hold on; grid on;
        case 8
            plot(J_5Hz_4_average_result(:,3),J_5Hz_4_average_result(:,2),'bo')
            hold on; grid on;
        case 9
            plot(J_10Hz_1_average_result(:,3),J_10Hz_1_average_result(:,2),'bo')
            hold on; grid on;
        case 10
            plot(J_10Hz_2_average_result(:,3),J_10Hz_2_average_result(:,2),'bo')
            hold on; grid on;
        case 11
            plot(J_10Hz_3_average_result(:,3),J_10Hz_3_average_result(:,2),'bo')
            hold on; grid on;
        case 12
            plot(J_10Hz_4_average_result(:,3),J_10Hz_4_average_result(:,2),'bo')
            hold on; grid on;
    end
end
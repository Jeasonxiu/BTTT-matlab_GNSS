clear all; close all;
tic;
%% ÁÂ¿ì ¾ÈÅ×³ª À§Ä¡ °á°ú load
file_L = 'test3_L_EstOut.mat';
file_R = 'test3_R_EstOut.mat';
load(file_L); load(file_R);

if file_L(end-5:end-4) == 'el'
    EstOut_L = EstOut_L_el;
    EstOut_R = EstOut_R_el;
    title_txt = 'EKF with el weighting';
else
    title_txt = 'EKF with SNR weighting';
end
%% VRS load
load('PTCO1_170116_adm.txt');
vrs= PTCO1_170116_adm;
vrs(:,1) = vrs(:,1)+18;

%% ±âÁØ ÁÂÇ¥
TruePos = [-3058799.61420451,4083265.35912516,3814946.87192938];
gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2);
%% °øÅë½Ã°£
FinalTTs = intersect(EstOut_L(:, 1), EstOut_R(:, 1));
FinalTTs1 = FinalTTs(1:592,1);
FinalTTs2 = FinalTTs(620:1140,1);
FinalTTs3 = FinalTTs(1158:1636,1);

for j=1:3
    switch j
        case 1
            tic;
            time = FinalTTs1;
            for i=1:length(time)
                gs = time(i);
                %% ÇØÅÁ ½Ã°¢ vrs, ÁÂ, ¿ì ¾ÈÅ×³ª ÁÂÇ¥°ª È¹µæ
                vrs_coordi = vrs(find(vrs(:,1)==gs),5:7);
                left = EstOut_L(find(EstOut_L(:,1) == gs),2:4);
                right = EstOut_R(find(EstOut_R(:,1) == gs),2:4);
                
                %% dNEV °è»ê
                dXYZ_L = left - TruePos;
                dXYZ_R = right - TruePos;
                dXYZ_vrs = vrs_coordi - TruePos;
                dNEV_L = xyz2topo(dXYZ_L, TrueLat, TrueLon);
                dNEV_R = xyz2topo(dXYZ_R, TrueLat, TrueLon);
                dNEV_vrs = xyz2topo(dXYZ_vrs, TrueLat, TrueLon);
                dN_L = dNEV_L(:,1); dE_L = dNEV_L(:,2); dV_L = dNEV_L(:,3);
                dN_R = dNEV_R(:,1); dE_R = dNEV_R(:,2); dV_R = dNEV_R(:,3);
                dN_vrs = dNEV_vrs(:,1); dE_vrs = dNEV_vrs(:,2); dV_vrs = dNEV_vrs(:,3);
                dNE_L = sqrt(dN_L^2 + dE_L^2);
                dNE_R = sqrt(dN_R^2 + dE_R^2);
                dNE_vrs = sqrt(dN_vrs^2 + dE_vrs^2);
                d2D_L(i,1) = dNE_L;
                d2D_R(i,1) = dNE_R;
                d2D_vrs(i,1) = dNE_vrs;
                d3_L = sqrt(dN_L.^2 + dE_L.^2 + dV_L.^2);
                d3_R = sqrt(dN_R.^2 + dE_R.^2 + dV_R.^2);
                d3_vrs = sqrt(dN_vrs.^2 + dE_vrs.^2 + dV_vrs.^2);
                d3D_L(i,1) = d3_L;
                d3D_R(i,1) = d3_R;
                d3D_vrs(i,1) = d3_vrs;
                result_L(i,1:6) = [gs, dN_L, dE_L, dV_L, dNE_L, d3_L];
                result_R(i,1:6) = [gs, dN_R, dE_R, dV_R, dNE_R, d3_R];
                result_vrs(i,1:6) = [gs, dN_vrs, dE_vrs, dV_vrs, dNE_R, d3_R];
                figure(j)
                title(title_txt)
                plot(dE_L,dN_L,'b.','MarkerSize',10)
                hold on; grid on;
                plot(dE_R,dN_R,'g.','MarkerSize',10)
                plot(dE_vrs, dN_vrs,'r.','MarkerSize',20)
                plot((dE_L+dE_R)/2,(dN_L+dN_R)/2,'*k')
                plot([dE_L,dE_R],[dN_L,dN_R],'m-')
                legend('Left','Right','VRS','Position','Location','Best')
                axis equal
%                 drawnow
            end
            toc;
        case 2
            tic;
            time = FinalTTs2;
            for i=1:length(time)
                gs = time(i);
                %% ÇØÅÁ ½Ã°¢ vrs, ÁÂ, ¿ì ¾ÈÅ×³ª ÁÂÇ¥°ª È¹µæ
                vrs_coordi = vrs(find(vrs(:,1)==gs),5:7);
                left = EstOut_L(find(EstOut_L(:,1) == gs),2:4);
                right = EstOut_R(find(EstOut_R(:,1) == gs),2:4);
                
                %% dNEV °è»ê
                dXYZ_L = left - TruePos;
                dXYZ_R = right - TruePos;
                dXYZ_vrs = vrs_coordi - TruePos;
                dNEV_L = xyz2topo(dXYZ_L, TrueLat, TrueLon);
                dNEV_R = xyz2topo(dXYZ_R, TrueLat, TrueLon);
                dNEV_vrs = xyz2topo(dXYZ_vrs, TrueLat, TrueLon);
                dN_L = dNEV_L(:,1); dE_L = dNEV_L(:,2); dV_L = dNEV_L(:,3);
                dN_R = dNEV_R(:,1); dE_R = dNEV_R(:,2); dV_R = dNEV_R(:,3);
                dN_vrs = dNEV_vrs(:,1); dE_vrs = dNEV_vrs(:,2); dV_vrs = dNEV_vrs(:,3);
                dNE_L = sqrt(dN_L^2 + dE_L^2);
                dNE_R = sqrt(dN_R^2 + dE_R^2);
                dNE_vrs = sqrt(dN_vrs^2 + dE_vrs^2);
                d2D_L(i,1) = dNE_L;
                d2D_R(i,1) = dNE_R;
                d2D_vrs(i,1) = dNE_vrs;
                d3_L = sqrt(dN_L.^2 + dE_L.^2 + dV_L.^2);
                d3_R = sqrt(dN_R.^2 + dE_R.^2 + dV_R.^2);
                d3_vrs = sqrt(dN_vrs.^2 + dE_vrs.^2 + dV_vrs.^2);
                d3D_L(i,1) = d3_L;
                d3D_R(i,1) = d3_R;
                d3D_vrs(i,1) = d3_vrs;
                result_L(i,1:6) = [gs, dN_L, dE_L, dV_L, dNE_L, d3_L];
                result_R(i,1:6) = [gs, dN_R, dE_R, dV_R, dNE_R, d3_R];
                result_vrs(i,1:6) = [gs, dN_vrs, dE_vrs, dV_vrs, dNE_R, d3_R];
                figure(j)
                title(title_txt)
                plot(dE_L,dN_L,'b.','MarkerSize',10)
                hold on; grid on;
                plot(dE_R,dN_R,'g.','MarkerSize',10)
                plot(dE_vrs, dN_vrs,'r.','MarkerSize',20)
                plot((dE_L+dE_R)/2,(dN_L+dN_R)/2,'*k')
                plot([dE_L,dE_R],[dN_L,dN_R],'m-')
                legend('Left','Right','VRS','Position','Location','Best')
                axis equal
%                 drawnow
            end
            toc;
        case 3
            tic;
            time = FinalTTs3;
            for i=1:length(time)
                gs = time(i);
                %% ÇØÅÁ ½Ã°¢ vrs, ÁÂ, ¿ì ¾ÈÅ×³ª ÁÂÇ¥°ª È¹µæ
                vrs_coordi = vrs(find(vrs(:,1)==gs),5:7);
                left = EstOut_L(find(EstOut_L(:,1) == gs),2:4);
                right = EstOut_R(find(EstOut_R(:,1) == gs),2:4);
                
                %% dNEV °è»ê
                dXYZ_L = left - TruePos;
                dXYZ_R = right - TruePos;
                dXYZ_vrs = vrs_coordi - TruePos;
                dNEV_L = xyz2topo(dXYZ_L, TrueLat, TrueLon);
                dNEV_R = xyz2topo(dXYZ_R, TrueLat, TrueLon);
                dNEV_vrs = xyz2topo(dXYZ_vrs, TrueLat, TrueLon);
                dN_L = dNEV_L(:,1); dE_L = dNEV_L(:,2); dV_L = dNEV_L(:,3);
                dN_R = dNEV_R(:,1); dE_R = dNEV_R(:,2); dV_R = dNEV_R(:,3);
                dN_vrs = dNEV_vrs(:,1); dE_vrs = dNEV_vrs(:,2); dV_vrs = dNEV_vrs(:,3);
                dNE_L = sqrt(dN_L^2 + dE_L^2);
                dNE_R = sqrt(dN_R^2 + dE_R^2);
                dNE_vrs = sqrt(dN_vrs^2 + dE_vrs^2);
                d2D_L(i,1) = dNE_L;
                d2D_R(i,1) = dNE_R;
                d2D_vrs(i,1) = dNE_vrs;
                d3_L = sqrt(dN_L.^2 + dE_L.^2 + dV_L.^2);
                d3_R = sqrt(dN_R.^2 + dE_R.^2 + dV_R.^2);
                d3_vrs = sqrt(dN_vrs.^2 + dE_vrs.^2 + dV_vrs.^2);
                d3D_L(i,1) = d3_L;
                d3D_R(i,1) = d3_R;
                d3D_vrs(i,1) = d3_vrs;
                result_L(i,1:6) = [gs, dN_L, dE_L, dV_L, dNE_L, d3_L];
                result_R(i,1:6) = [gs, dN_R, dE_R, dV_R, dNE_R, d3_R];
                result_vrs(i,1:6) = [gs, dN_vrs, dE_vrs, dV_vrs, dNE_R, d3_R];
                figure(j)
                title(title_txt)
                plot(dE_L,dN_L,'b.','MarkerSize',10)
                hold on; grid on;
                plot(dE_R,dN_R,'g.','MarkerSize',10)
                plot(dE_vrs, dN_vrs,'r.','MarkerSize',20)
                plot((dE_L+dE_R)/2,(dN_L+dN_R)/2,'*k')
                plot([dE_L,dE_R],[dN_L,dN_R],'m-')
                legend('Left','Right','VRS','Position','Location','Best')
                axis equal
%                 drawnow
            end
            toc;
    end
end
 toc;           
         
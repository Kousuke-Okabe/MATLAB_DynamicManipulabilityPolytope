%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   行列形DMP並進とテンソル形DMP並進の比較
%
%                                                           19.06.11 by OKB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

param = get_parameter;
    Tlim = param.Tlim;
    N = param.Link;

% DMP表示倍率
Gr = 1/500;
Gz = Gr/20;

dz = [0,35,50];

FH = 1;
figure(FH)
clf(FH)
for i = 1:size(dz,2)
    % Motion on ExOS
    x = [ 0.2; 0.2; pi/6 ];
    dx = [ 0; 0; dz(i)];

    % Motoion on Joint space
    q = fIKinematics(x(1:N-1),x(N));
    
     fDraw_DMP(FH,1, Gr,Gz, x,dx)
%     fDraw_DMP(FH,0, Gr,Gz, x,dx)
    Vt = get_TranslateVector(x,dx);
    text(x(1)+Vt(1)*Gr,x(2)+Vt(2)*Gr,Vt(3)*Gz,strcat('$$J_{C,ijk}\dot{x}_j\dot{x}_k$$ with $$\dot{z}$$=',num2str(dx(3))) ,'interpreter','latex', 'fontsize',20)
    
    
    fDraw_DMP_v1(FH,1, Gr,Gz, x,dx)
%     fDraw_DMP_v1(FH,0, Gr,Gz, x,dx)
    Vt = get_TranslateVector_v1(x,dx);
    text(x(1)+Vt(1)*Gr,x(2)+Vt(2)*Gr,Vt(3)*Gz,strcat('$$J_{Z}\dot{z}$$ with $$\dot{z}$$=',num2str(dx(3))) ,'interpreter','latex', 'fontsize',20)
    
end

fRoboAnimation(FH, q,x(1:N-1),0)


set(gca, 'FontSize',20)
xlim([-0.0, 0.25])
% xlim([-1, 1]*0.06+ones(1,2)*0.2)
ylim(xlim)
zlim([-1, 1]*0.05)
view([-25,20])
% view([0,90])

xlabel('$$\ddot{r}_1[m/s^2]$$', 'FontSize',20, 'interpreter','latex')
ylabel('$$\ddot{r}_2[m/s^2]$$', 'FontSize',20, 'interpreter','latex')
zlabel('$$\ddot{z}[s^{-2}]$$', 'FontSize',20, 'interpreter','latex')

% plot3(xlim,[x(2),x(2)],[0,0],'k')
% plot3([x(1),x(1)],ylim,[0,0],'k')
% plot3([x(1),x(1)],[x(2),x(2)],zlim,'k')

box off
grid off

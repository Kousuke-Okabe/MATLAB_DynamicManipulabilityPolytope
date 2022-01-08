%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ç’·‘¬“xdz‚ð•Ï‰»‚³‚¹‚½‚Æ‚«‚É•ÀiƒxƒNƒgƒ‹‚ª•`‚­‹OÕ
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

FH = 1;
Gr = 1/150;
Gz = Gr/20;

r = [0.2; 0.2];
z = pi/6;
x = [r;z];

q = zeros(3,1);
dq = zeros(3,1);
q(:,1) = fIKinematics(r(:,1),z);

drawDMP3d = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param = get_parameter();
    Tlim = param.Tlim;
    N = param.Link;

Gr = 1/150;
Gz = Gr/20;

nState = 101;
nState_r = 3;
nState_z = 5;

dxref = [ linspace(0,2,nState);
          linspace(0,4,nState);
          linspace(-25,25,nState) ];


Vt = nan(3,nState,nState,nState);

for alpha = 1:nState
for beta = 1;%linspace(1,nState,nState_z)
for gamma = round(linspace(1,nState,nState_z))
    
    dx = [ dxref(1,alpha);
           dxref(2,beta);
           dxref(3,gamma) ];
    Vt(:,alpha,beta,gamma) = get_TranslateVector(x,dx);
    
end
end
end

% for alpha = linspace(1,nState,nState_z)
% for beta = 1:nState
% for gamma = linspace(1,nState,nState_z)
%     
%     dx = [ dxref(1,alpha);
%            dxref(2,beta);
%            dxref(3,gamma) ];
%     Vt(:,alpha,beta,gamma) = get_TranslateVector(x,dx);
%     
% end
% end
% end

for alpha = round(linspace(1,nState,nState_r))
for beta = 1;%linspace(1,nState,nState_z)
for gamma = 1:nState
    
    dx = [ dxref(1,alpha);
           dxref(2,beta);
           dxref(3,gamma) ];
    Vt(:,alpha,beta,gamma) = get_TranslateVector(x,dx);
    
end
end
end



%% •`‰æ
clf(FH)
fRoboAnimation(FH,q, r,0)
hold on

for beta = 1;%linspace(1,nState,nPState)
for gamma = round(linspace(1,nState,nState_z))
    plot3(r(1)+Vt(1,:,beta,gamma)*Gr,r(2)+Vt(2,:,beta,gamma)*Gr,Vt(3,:,beta,gamma)*Gz,'r', 'LineWidth',2.0)
    
    if sum(gamma == round(linspace(1,nState,nState_z))) > 0
        text(r(1)+Vt(1,1,beta,gamma)*Gr,r(2)+Vt(2,1,beta,gamma)*Gr,Vt(3,1,beta,gamma)*Gz, ...
            strcat('$$\dot{x}$$=(', num2str(dxref(1,1)),',',num2str(0),',',num2str(dxref(3,gamma)) ,')'), 'interpreter','latex', 'fontsize',20)
        text(r(1)+Vt(1,nState,beta,gamma)*Gr,r(2)+Vt(2,nState,beta,gamma)*Gr,Vt(3,nState,beta,gamma)*Gz, ...
            strcat('$$\dot{x}$$=(', num2str(dxref(1,nState)),',',num2str(0),',',num2str(dxref(3,gamma)) ,')'), 'interpreter','latex', 'fontsize',20)
    end
end
end

% for alpha = 1;%linspace(1,nState,nPState)
% for gamma = linspace(1,nState,nPState)
%     plot3(r(1)+permute(Vt(1,alpha,:,gamma),[1,3,2,4])*Gr,r(2)+permute(Vt(2,alpha,:,gamma),[1,3,2,4])*Gr,permute(Vt(3,alpha,:,gamma),[1,3,2,4])*Gz,'g')
% %     text(r(1)+Vt(1,alpha,beta,1)*Gr,r(2)+Vt(2,alpha,beta,1)*Gr,Vt(3,alpha,beta,1)*Gz, strcat('dr1=',num2str(dxref(1,alpha)),',dz=',num2str(dxref(3,1))))
% %     text(r(1)+Vt(1,alpha,beta,nState)*Gr,r(2)+Vt(2,alpha,beta,nState)*Gr,Vt(3,alpha,beta,nState)*Gz, strcat('dr1=',num2str(dxref(1,alpha)),',dz=',num2str(dxref(3,nState))))
% 
% end
% end

for alpha = round(linspace(1,nState,nState_r))
for beta = 1;%linspace(1,nState,nPState)
    plot3(r(1)+permute(Vt(1,alpha,beta,:),[1,4,2,3])*Gr,r(2)+permute(Vt(2,alpha,beta,:),[1,4,2,3])*Gr,permute(Vt(3,alpha,beta,:),[1,4,2,3])*Gz,'b', 'LineWidth',2.0)
    
	if sum(alpha == round(linspace(1,nState,nState_r))) > 0
        text(r(1)+Vt(1,alpha,beta,1)*Gr,r(2)+Vt(2,alpha,beta,1)*Gr,Vt(3,alpha,beta,1)*Gz, ...
            strcat('$$\dot{x}$$=(', num2str(dxref(1,alpha)),',',num2str(0),',',num2str(dxref(3,1)) ,')'), 'interpreter','latex', 'fontsize',20)
        text(r(1)+Vt(1,alpha,beta,nState)*Gr,r(2)+Vt(2,alpha,beta,nState)*Gr,Vt(3,alpha,beta,nState)*Gz, ...
            strcat('$$\dot{x}$$=(', num2str(dxref(1,alpha)),',',num2str(0),',',num2str(dxref(3,nState)) ,')'), 'interpreter','latex', 'fontsize',20)
    end
end
end

set(gca, 'FontSize',20)
xlabel('$$\ddot{r}_1[m/s^2]$$', 'FontSize',20)
ylabel('$$\ddot{r}_2[m/s^2]$$', 'FontSize',20)
zlabel('$$\ddot{z}[s^{-2}]$$', 'FontSize',20)

view([-20,20])
xlim([-0.05, 0.2])
ylim([-0.1, 0.25])
zlim([-0.1, 0.05])

% plot3(xlim,[0,0],[0,0],'k')
% plot3([0,0],ylim,[0,0],'k')
% plot3(xlim,[1,1]*max(ylim),[0,0],'k')
% plot3(xlim,[1,1]*min(ylim),[0,0],'k')
% plot3([1,1]*max(xlim),ylim,[0,0],'k')
% plot3([1,1]*min(xlim),ylim,[0,0],'k')
% plot3([0,0],[0,0],zlim,'k')

box off


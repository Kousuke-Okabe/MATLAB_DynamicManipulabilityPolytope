%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   冗長速度dzを変化させたときに並進ベクトルが描く軌跡
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

FH = 1;
Gr = 1/150;
Gz = Gr/20;

r = [0.2; 0.2];
z = pi/6;


q = zeros(3,1);
dq = zeros(3,1);
q(:,1) = fIKinematics(r(:,1),z);
fRoboAnimation(FH,q, r,0)

drawDMP3d = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(FH)
hold on

param = get_parameter();
    Tlim = param.Tlim;
    N = param.Link;

Gr = 1/150;
Gz = Gr/20;

nState = 100;
nVector = 5;
dr = [2;0];
dz = linspace(-25,25,nState);

V = nan(3,nState);
for alpha = 1:nState
    x = [r;z];
    dx = [dr; dz(alpha)];

    [J,Jp,Jd,U,Ud,err] = fJacobi_tensor(x(1:2),dx(1:2), x(3),dx(3));
    Je = [J;U'];
    Jei = [Jp,U];

    q = fIKinematics(x(1:2),x(3));
    % dq = Jp*dr(:,1)+U*dz;
    dq = Jei*dx;

    [M,h,g,hd] = get_matrix(q,dq);

    % Transpose Vector of DMP
    Jed = Jd;
    Jed(3,:,:) = Ud;

    Jck = zeros(N,N,N);
    Jch = zeros(N,N,N);
    Mi =inv(M);
    for i = 1:N
        for j = 1:N
            for k = 1:N
                for n = 1:N
                    for o = 1:N
                        Jck(i,n,o) = Jck(i,n,o) + Jed(i,j,k)*Jei(j,n)*Jei(k,o);
                    end
                end
            end
        end
    end
    for i = 1:N
        for j = 1:N
            for k = 1:N
                for l = 1:N
                    for m = 1:N
                        for n = 1:N
                            for o = 1:N
                                Jch(i,n,o) = Jch(i,n,o) - Je(i,j)*Mi(j,k)*hd(k,l,m)*Jei(l,n)*Jei(m,o);
                            end
                        end
                    end
                end
            end
        end
    end

    DMP_t = zeros(N,1);
    for i = 1:N
        for j = 1:N
            for k = 1:N
                DMP_t(i) = DMP_t(i) + Jck(i,j,k)*dx(j)*dx(k) + Jch(i,j,k)*dx(j)*dx(k);
            end
        end
    end

    V(:,alpha) = DMP_t;
    
%     % Drawing transpose vector
%     if alpha == 1 || alpha == nState || alpha == round(nState/2)
%         quiver3(r(1),r(2),0, DMP_t(1)*Gr,DMP_t(2)*Gr,DMP_t(3)*Gz,'g', 'AutoScaleFactor',1)
%     end
end


plot3(r(1)+V(1,:)*Gr,r(2)+V(2,:)*Gr,V(3,:)*Gz,'r')
t = strcat('$$\dot{z}$$=',num2str(dz(1)));
% text(r(1)+V(1,1)*Gr,r(2)+V(2,1)*Gr,V(3,1)*Gz, strcat('$$\leftarrow \dot{z}$$=',num2str(dz(1))), 'Interpreter','latex', 'FontSize',20)
% text(r(1)+V(1,nState)*Gr,r(2)+V(2,nState)*Gr,V(3,nState)*Gz, strcat('$$\leftarrow \dot{z}$$=',num2str(dz(nState))), 'Interpreter','latex', 'FontSize',20)

% title(strcat('$$\dot{\bf r}$$=[',num2str(dx(1)),',',num2str(dx(2)),']$$^T$$ [m/s]'), 'interpreter','latex')
view([-20,20])
xlim([-0.1, 0.3])
ylim(xlim)
zlim([-0.1, 0.1])

plot3(xlim,[0,0],[0,0],'k')
plot3(xlim,[1,1]*max(ylim),[0,0],'k')
plot3(xlim,[1,1]*min(ylim),[0,0],'k')
plot3([0,0],ylim,[0,0],'k')
plot3([1,1]*max(xlim),ylim,[0,0],'k')
plot3([1,1]*min(xlim),ylim,[0,0],'k')
% plot3([0,0],[0,0],zlim,'k')
 
box off

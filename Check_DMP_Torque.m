%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   動作加速度ベクトルとDMPの関係を駆動トルクにより確認
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

param = get_parameter;
    Tlim = param.Tlim;
    N = param.Link;

% DMP表示倍率
Gr = 1/250;
Gz = Gr/20;

% Motion
r = [0.2; 0.2];
dr = [0;0];
ddr = [-20;-10];

z = 0;
% dz =0;
dz = 25;
ddz =0/Gz;

x = [r;z];
dx = [dr;dz];
ddx = [ddr;ddz];


% Calculating Driven Torque
q = fIKinematics(r,z);
[J,Jp,Jd,U,Ud,err] = fJacobi_tensor(r,dr, z,dz);
dq = Jp*dr + U*dz;

Jei = [Jp,U];
Jed = Jd;
Jed(3,:,:) = Ud;


ddq = Jei*ddx;
for i = 1:N
    for j = 1:N
        for k = 1:N
            for l = 1:N
                for m = 1:N
                    for n = 1:N
                        ddq(i) = ddq(i) - Jei(i,j)*Jed(j,k,l)*Jei(k,m)*Jei(l,n)*dx(m)*dx(n);
                    end
                end
            end
        end
    end
end

[M,H,G,hd] = get_matrix_plus(q,dq);

T = M*ddq + H + G
T_n= sqrt(T'*T)
Tlim

FH = 1;
figure(FH)
clf(FH)
hold on
fRoboAnimation(FH, q,r,0)
% fDraw_DMP(FH,1, Gr,Gz, x,dx)
% fDraw_DMP(FH,1, Gr,Gz, x,[dr;0])
fDraw_DMP(FH,1, Gr,Gz, x,dx)
fDraw_DMP(FH,1, Gr,Gz, x,[dr;0])

quiver3(r(1),r(2),0, ddx(1)*Gr,ddx(2)*Gr,ddx(3)*Gz,'b', 'AutoScaleFactor',1, 'LineWidth',1.5)


xlim([-0, 0.3])
ylim(xlim)
zlim([-1, 1]*0.1)
view([0,90])

set(gca, 'FontSize',20)
xlabel('$$\ddot{r}_1[m/s^2]$$', 'FontSize',20)
ylabel('$$\ddot{r}_2[m/s^2]$$', 'FontSize',20)
zlabel('$$\ddot{z}[s^{-2}]$$', 'FontSize',20)

legend('Acceleration Vector','Translating Vector')

% plot3(xlim,[0,0],[0,0],'k')
% plot3([0,0],ylim,[0,0],'k')
% plot3([0,0],[0,0],zlim,'k')

grid off

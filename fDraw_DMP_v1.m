function[] = fDraw_DMP_v1(FH,drawDMP3d, Gr,Gz, x,dx)

%**************************************************************************
%
%   [] = fDraw_DMP(FH,drawDMP3d, Gr,Gz, x,dx)
%
%       FH  : フィギュアーハンドル
%       drawDMP3d : Flag of drawingDMP on 3Dspace. yes:1 no:0
%       Gr  : r空間描画サイズ
%       Gz  : z空間描画サイズ
%       x   : Position on ExOS
%       dx  : Velocity on ExOS
%
%       v2 : デティールアップバージョン
%       v3 : 直動系対応
%       v4 : 動的可操作性多面体描画
%
%**************************************************************************

%%
clear

FH = 1;
Gr = 1/150;
Gz = Gr/20;

x = [0.2; 0.1; 0];
dx = [0; 0; 0];

q = zeros(3,1);
dq = zeros(3,1);
q(:,1) = fIKinematics(x(1:2),x(3));

figure(FH)
clf(FH)
fRoboAnimation(FH,q, x(1:2),0)

drawDMP3d = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(FH)
hold on

param = get_parameter();
    Tlim = param.Tlim;
    N = param.Link;

% [J,Jp,Jd,U,Ud,err] = fJacobi_tensor(r,dr, z,dz);
[J,Jp,dJ,U,error] = fJacobi_plus(x(1:N-1),dx(1:N-1), x(N),dx(N));
Je = [J;U'];
Jei = [Jp,U];

q = fIKinematics(x(1:N-1),x(N));
% dq = Jp*dr(:,1)+U*dz;
dq = Jei*dx;

[M,h,g] = get_matrix(q,dq);


% Vertex of DMP
Jt = Je/M;
DMP_v = nan(size(Jt,1),factorial(size(Jt,1)));

a = 0;
for i = -1:2:1
    for j = -1:2:1
        for k = -1:2:1
            a = a+1;
            % Previous : h & g vector in the T'.
            DMP_v(:,a) = Jt(:,1)*(i*Tlim(1)-g(1)-h(1)) + Jt(:,2)*(j*Tlim(2)-g(2)-h(2)) + Jt(:,3)*(k*Tlim(3)-g(3)-h(3));
        end
    end
end

% Calclating Basis Vector of DMP
% DMP_t = [ dJ*Jp, dJ*U; -U'*dJ'/(J*J'), 0 ]*dx;
DMP_t = get_TranslateVector_v1(x,dx);

%         % Drawing Basis vector of Transpose Vector
%         figure(2)
%         hold on
%         i = 3;
%         DMP_te = [ dJ*Jp, dJ*U; -U'*dJ'/(J*J'), 0 ];
%         DMP_te = DMP_te(:,3);
%         temp = DMP_te/norm(DMP_te)
%         quiver3(0,0,0, DMP_te(1)*Gr,DMP_te(2)*Gr,DMP_te(3)*Gz,'b', 'LineWidth',2,'AutoScaleFactor',1)

% Adding transpose vector to DMP
for i = 1:size(DMP_v,2)
    DMP_v(:,i) = DMP_v(:,i)+DMP_t;
end

% Drawing transpose vector
figure(FH)
quiver3(x(1),x(2),0, DMP_t(1)*Gr,DMP_t(2)*Gr,DMP_t(3)*Gz,'b', 'LineWidth',2,'AutoScaleFactor',1)

% DMP drawing 3dim space   
if drawDMP3d
    V = convhull(DMP_v(1,:),DMP_v(2,:),DMP_v(3,:));
    for i = 1:size(V,1)
        patch(x(1)+DMP_v(1,V(i,:))*Gr,x(2)+DMP_v(2,V(i,:))*Gr,DMP_v(3,V(i,:))*Gz,'k', 'FaceAlpha',0 ,'EdgeColor','b','EdgeAlpha',1)
    end
    
    figure(FH)
    zlabel('$$\dot{z}$$', 'FontSize',21, 'interpreter','latex')
    view([-20,40])
    rotate3d on
end




% DMP drawing 2dim space
if drawDMP3d == 0
% for i = 1:size(DMP_v,2)
%     text(r(1)+DMP_v(1,i)*Gr,r(2)+DMP_v(2,i)*Gr,DMP_v(3,i)*Gz,num2str(i), 'FontSize',12)
% end

DMP0_v = nan(3,12);
EdgeNumber = [1,2; 1,3; 1,5; 2,4; 2,6; 3,4; 3,7; 4,8; 5,6; 5,7; 6,8; 7,8];

l = 0;
for k = 1:size(EdgeNumber,1)
    i = EdgeNumber(k,1);
    j = EdgeNumber(k,2);
    
    if DMP_v(3,i)*DMP_v(3,j) < 0
        l = l+1;
        DMP0_v(:,l) = ( abs(DMP_v(3,j))*DMP_v(:,i)+abs(DMP_v(3,i))*DMP_v(:,j) )/( abs(DMP_v(3,i))+abs(DMP_v(3,j)) );
    end
end

DMP0 = convhull(DMP0_v(1,1:l),DMP0_v(2,1:l));
plot(x(1)+DMP0_v(1,DMP0)*Gr,x(2)+DMP0_v(2,DMP0)*Gr,'g', 'LineWidth',2)

end
end

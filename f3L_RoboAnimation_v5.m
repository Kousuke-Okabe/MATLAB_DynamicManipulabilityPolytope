function[] = f3L_RoboAnimation_v5(FH,dim, Gr,Gz, r,dr,iz,z)

%**************************************************************************
%
%   [] = f3L_RoboAnimation_v5(FH,dim, Gr,Gz, r,dr,iz,z)
%
%       FH  : フィギュアーハンドル
%       dim : DMP描画次元
%       Gr  : r空間描画サイズ
%       Gz  : z空間描画サイズ
%       r   : 手先位置ベクトル
%       dr  : 手先速度ベクトル
%       iz  : Integral(z)
%       z   : 冗長速度
%
%       v2 : デティールアップバージョン
%       v3 : 直動系対応
%       v4 : 動的可操作性多面体描画
%
%                                                       18.03.01 by OKB
%**************************************************************************

%%
clear

FH = 1;
Gr = 1/150;
Gz = Gr/20;

r = [ 0.1654; 0.1362];
dr = [ 0; 0 ];

iz = 3.8445;
z = 18.0000;


q = zeros(3,1);
dq = zeros(3,1);
q(:,1) = f3L_iKinematics(r(:,1),iz);
f3L_RoboAnimation_v3(FH,q, r,0)

dim = 3; % 3:Drawing 3dim  2:Drawing 2dim

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(FH)
hold on

ddr_d = zeros(3,6);
d = zeros(12,3);
d2 = zeros(6,3);
d2_t = zeros(3,6);
t = zeros(1,3);

param = get_parameter();
    Tlim = param.Tlim;
     
[J,Jp,dJ,u] = f3L_Jacobi_plus(r,dr,iz,z);
Je = [J;u'];

q = f3L_iKinematics(r,iz);
dq = Jp*dr(:,1)+u*z;

[M,h,g] = get_matrix(q,dq);

Ut = Je/M;
Ur = [ dJ*Jp; u'*dJ'/(J*J') ];
Uz = [ dJ*u; zeros(1,size(u,2)) ];

switch dim
    % DMP drawing 2dim space
    case 2
        for p = 1:3
            v = diag([1,1,1]);
            a(1,p) = u.'/M*(v(p,:)).';
        end
        for n = 1:3
            i = 0;
            for k = -1:2:1
                for l = -1:2:1
                    i = i+1;
                    if n == 1
                        Tc = diag([0,k,l]);
                        t(1,n) = (u.'/M*(h+g-Tc*Tlim))/a(1,n);
                        d(i,:) = (Tc*Tlim-h-g).' + v(1,:)*t(1,1);
                    elseif  n == 2
                        Tc = diag([k,0,l]);
                        t(1,n) = (u.'/M*(h+g-Tc*Tlim))/a(1,n);
                        d(i+4,:) = (Tc*Tlim-h-g).' + v(2,:)*t(1,2);
                    else 
                        Tc = diag([k,l,0]);
                        t(1,n) = (u.'/M*(h+g-Tc*Tlim))/a(1,n);
                        d(i+8,:) = (Tc*Tlim-h-g).' + v(3,:)*t(1,3);
                    end
                end
            end
        end
        
        b = 1;
        for k = 1:4
            if d(k,1) <= Tlim(1,1) && d(k,1) >= -Tlim(1,1)
                d2(b,:) = d(k,:);
                b = b+1;
            end
        end
        for k = 5:8
            if d(k,2) <= Tlim(2,1) && d(k,2) >= -Tlim(2,1)
                d2(b,:) = d(k,:);
                b = b+1;
            end
        end
        for k = 9:12
            if d(k,3) <= Tlim(3,1) && d(k,3) >= -Tlim(3,1)
                d2(b,:) = d(k,:);
                b = b+1;
            end
        end

        d2_t = d2'; 

        for k = 1:6
            ddr_d(:,k) = Ut*d2_t(:,k) + Ur*dr + Uz*z;
        end

        Ur2 = dJ*u*z;
        % for k = 1:6
        %     quiver(r(1,1)+Ur2(1,1)*G, r(2,1)+Ur2(2,1)*G, ddr_d(1,k)*G-Ur2(1,1)*G ,ddr_d(2,k)*G-Ur2(2,1)*G,'g','AutoScaleFactor',1.0, 'LineWidth',1);
        % end

        N = convhull(ddr_d(1,:)+Ur2(1,1)*Gr ,ddr_d(2,:)+Ur2(2,1)*Gr);
        plot(r(1,1)+ddr_d(1,N)*Gr, r(2,1)+ddr_d(2,N)*Gr,'g');

        % quiver(r(1,1),r(2,1),ddr(1,1)*G,ddr(2,1)*G,'r','AutoScaleFactor',1.0,'LineWidth',1.5);
        quiver(r(1,1),r(2,1),Ur2(1,1)*Gr,Ur2(2,1)*Gr,'b','AutoScaleFactor',1.0,'LineWidth',1.5);
     
        
        
    % DMP drawing 3dim space   
    case 3
        Urdr = Ur*dr;
        
        % Drawing Jef Mi Td
        Td = Tlim-h-g;
        for i = 1:size(Tlim,1)
            quiver3(r(1),r(2),0, Ut(1,i)*Td(i)*Gr,Ut(2,i)*Td(i)*Gr,Ut(3,i)*Td(i)*Gz,'g', 'AutoScaleFactor',1)
        end
        
        % Drawing transpose vector
        quiver3(r(1),r(2),0, (Urdr(1)+Uz(1)*z)*Gr,(Urdr(2)+Uz(2)*z)*Gr,(Urdr(3)+Uz(3)*z)*Gz,'b', 'AutoScaleFactor',1)
        
        % Drawing DMP
        dmp_v = nan(size(Ut,1),factorial(size(Ut,1)));
        a = 0;
        for i = -1:2:1
            for j = -1:2:1
                for k = -1:2:1
                    a = a+1;
                    dmp_v(:,a) = Ut(:,1)*(i*Tlim(1)-h(1)-g(1)) + Ut(:,2)*(j*Tlim(2)-h(2)-g(2)) + Ut(:,3)*(k*Tlim(3)-h(3)-g(3));
                end
            end
        end
        
        N = convhull(dmp_v(1,:),dmp_v(2,:),dmp_v(3,:));
        for i = 1:size(N,1)
            patch(r(1)+(Urdr(1)+Uz(1)*z+dmp_v(1,N(i,:)))*Gr,r(2)+(Urdr(2)+Uz(2)*z+dmp_v(2,N(i,:)))*Gr,(Urdr(3)+Uz(3)*z+dmp_v(3,N(i,:)))*Gz, 'r', 'FaceAlpha',0 ,'EdgeColor','r','EdgeAlpha',1)
        end
        
        zlabel('$$\dot{z}$$', 'FontSize',21, 'interpreter','latex')
end

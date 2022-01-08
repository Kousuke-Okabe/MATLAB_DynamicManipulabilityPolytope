function[Vt] = get_TranslateVector(x,dx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Vt = get_TranslateVector(x,dx)
%
%       Vt : Transpose Vector
%       x  : Position on ExOS
%       dx : Velocity on ExOS
%                                                         19.06.11 by OKB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%
% clear
% 
% x = [ 0.2; 0.2; pi/6 ];
% dx = [ -2; 0; -25 ]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param = get_parameter();
    N = param.Link;

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

Vt = DMP_t;

%         % Drawing Basis vector of Transpose Vector
%         Gr = 1/500;
%         Gz = Gr/20;
%         figure(2)
%         hold on
%         i = 3;
%         DMP_te = (Jck(:,i,i)+Jch(:,i,i))
% %         DMP_te = DMP_t/dx(3)^2;
%         quiver3(0,0,0, DMP_te(1)*Gr,DMP_te(2)*Gr,DMP_te(3)*Gz,'m', 'LineWidth',2,'AutoScaleFactor',1)

end
    

function[Vt] = get_TranslateVector_v1(x,dx)

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
% dx = [ 2; 0; -25 ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param = get_parameter();
    N = param.Link;
    
[J,Jp,dJ,U,error] = fJacobi_plus(x(1:2),dx(1:2), x(3),dx(3));

DMP_t = [ dJ*Jp, dJ*U; -U'*dJ'/(J*J'), 0 ]*dx;

Vt = DMP_t;

%         % Drawing Basis vector of Transpose Vector
%         Gr = 1/500;
%         Gz = Gr/20;
%         figure(2)
%         hold on
%         i = 3;
%         DMP_te = (Jck(:,i,i)+Jch(:,i,i));
% %         DMP_te = DMP_t/dx(3)^2;
%         quiver3(0,0,0, DMP_te(1)*Gr,DMP_te(2)*Gr,DMP_te(3)*Gz,'m', 'LineWidth',2,'AutoScaleFactor',1)

end
    

function[Vt] = DMPTranslate(type1,type2,g,dr,f,q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   MFPの並進ベクトル
%
%   [Vt] = DMPTranslate(type1,type2,g,dr,f,q)
%       x   : 拡張作業空間上の位置
%       dx  : 拡張作業空間上の速度
%       g   : 重力ベクトル
%       f   : 拡張作業空間上の力
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%
%  clear all
% % % 
%  type1 = '2d_RRRR';
%  type2 = '2d_RR';
%  r = [0.2; 0.2 ;0.1 ;0.15];  %%[r1(メインタスクの手先位置) ; r2(サブタスクの手先位置)]
%  dr = [5; 5; 5; 5];  %%[dri(メインタスクの手先速度) ; dr2(サブタスクの手先速度)]
% %  g = [0;0;0;0];
% % % f = [0; 0; 0; 0];
% % draw = 0;
%  q = fIKinematics_msDMP(type1, r);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param = get_parameter(type1);
    Link = param.Link;    
    Vt = zeros(size(dr));
    dr = [-1; 1; 0; 0];
    [J1,Jp1] = fJacobi_minus_q(type1,q);
    [J2,Jp2] = fJacobi_minus_q(type2,q);
    J2 = [J2,zeros(2,2)];
    Jp2 = J2'/(J2*J2');
    I=eye(4);
    JI2 = (I-Jp1*J1)*Jp2;
    JI = [Jp1,JI2];
    JpI = inv(JI);
    [JIid] = fJacobi_tensor_4(type1,type2,q);
    
    
    dq = JI*dr;
    [M,H,G,Hd] = get_matrix_plus(type1,q,dq);
    Mi = inv(M);

for i = 1:Link
for j = 1:Link
for l = 1:Link
for n = 1:Link
for o = 1:Link
    Vt(i) = Vt(i) - JpI(i,j)*JIid(j,n,l)*dr(n)*JI(l,o)*dr(o);
end
end
end
end
end

for i = 1:Link
for j = 1:Link
for k = 1:Link
for l = 1:Link
for m = 1:Link
for n = 1:Link
for o = 1:Link
    Vt(i) = Vt(i) - JpI(i,j)*Mi(j,k)*Hd(k,l,m)*JI(l,n)*dr(n)*JI(m,o)*dr(o);
end
end
end
end
end
end
end






    

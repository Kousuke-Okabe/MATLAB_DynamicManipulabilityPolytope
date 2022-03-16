function[Vt] = DMPTranslate(type1,type2,g,dr,f,q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   MFP�̕��i�x�N�g��
%
%   [Vt] = DMPTranslate(type1,type2,g,dr,f,q)
%       x   : �g����Ƌ�ԏ�̈ʒu
%       dx  : �g����Ƌ�ԏ�̑��x
%       g   : �d�̓x�N�g��
%       f   : �g����Ƌ�ԏ�̗�
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%
%  clear all
% % % 
%  type1 = '2d_RRRR';
%  type2 = '2d_RR';
%  r = [0.2; 0.2 ;0.1 ;0.15];  %%[r1(���C���^�X�N�̎��ʒu) ; r2(�T�u�^�X�N�̎��ʒu)]
%  dr = [5; 5; 5; 5];  %%[dri(���C���^�X�N�̎�摬�x) ; dr2(�T�u�^�X�N�̎�摬�x)]
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






    

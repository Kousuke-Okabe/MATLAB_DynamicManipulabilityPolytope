clear all

type1 = '2d_RRRR'
type2 = '2d_RR'

% set parameter
param = get_parameter(type1);
    Link = param.Link; 
r = [ 0.2;
     0.2;
     0.1;
     0.15 ];
dr =  [-1; 1; 0; 0];
g = [0; 0; 0; 0];
f = [0; 0; 0; 0];
q = fIKinematics_msDMP(type1, r);
draw =1;


% DMPåvéZ
[Graph,V] = DMPGraph(type1,type2,q); %DMPÇÃí∏ì_ÇãÅÇﬂÇÈ
V;
[Vt] = DMPTranslate(type1,type2,dr,g,f,q); %DMPÇÃï¿êi
Vt
for i = 1:size(V,2)
    V(:,i) = V(:,i) + Vt; % + [x(1:2);0];
end
V


[J1,Jp1] = fJacobi_minus_q(type1,q);
[J2,Jp2] = fJacobi_minus_q(type2,q);
    J2 = [J2,zeros(2,2)];
    Jp2 = J2'/(J2*J2');
    I=eye(4);
JI2 = (I-Jp1*J1)*Jp2 ;
JI = [Jp1,JI2];
JpI = inv(JI);
dq = JI*dr;
[JIid] = fJacobi_tensor_4(type1,type2,q);
    [M,H,G,Hd] = get_matrix_plus(type1,q,dq);
    
    Out =[  -23.2282;     4.3032 ;  208.8060 ;   14.6839  ]
q = zeros(size(dr));
for i = 1:Link
for j = 1:Link
  q(i) =q(i)+JI(i,j)*Out(j);
end
end

for i = 1:Link
for j = 1:Link
for l = 1:Link
for k = 1:Link
    q(i) =q(i)+JIid(i,j,k)*dr(j)*JI(k,l)*dr(l);
end
end
end
end

tau =M*q +H
     
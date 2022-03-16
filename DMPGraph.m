 function[G,V] = DMPGraph(type1,type2,q)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %   DMPの頂点グラフをReturnする関数
% %
% %   [G,V] = get_DMPgraph(type1,type2,q)
% %       G : グラフ
% %       V : 頂点データ
% %       x : 作業空間上の位置
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%
% clear all

% type1 = '2d_RRRR'
% type2 = '2d_RR'
% r = [ 0.2;
%      0.2;
%      0.1;
%      0.15 ];
% dr = [0; 0; 0; 0];
% q = fIKinematics_msDMP(type1, r);
% draw =0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 変数宣言
param = get_parameter(type1);
    Link = param.Link;
    Tlim = param.Tlim;

% 頂点データ作成
V = nan(Link,2^Link);
for i = 1:Link
    for j = 1:size(V,2)
        if j-floor((j-1)/(2^(Link+1)/2^i))*(2^(Link+1)/2^i) <= 2^Link/2^i
            V(i,j) = Tlim(i);
        else
            V(i,j) = -Tlim(i);
        end
    end
end

% グラフ作成
G = zeros(2^Link);    % グラフ

for i = 1:size(G,2)
    for j = i+1:size(G,2)
        temp = 0;
        for k = 1:Link
            if V(k,i)*V(k,j) < 0
                temp = temp+1;
            end
        end
        
        if temp == 2
            G(i,j) = 1;
            G(j,i) = 1;
        end
    end
end

% MFP頂点作成
 [M,g] = get_matrix_minus(type1, q);

for i = 1:size(V,2)
    V(:,i) = V(:,i)-g;
end

[J1,Jp1] = fJacobi_minus_q(type1,q);
[J2,Jp2] = fJacobi_minus_q(type2,q);
    J2 = [J2,zeros(2,2)];
    Jp2 = J2'/(J2*J2');
    I=eye(4);
JI2 = (I-Jp1*J1)*Jp2 ;
JI = [Jp1,JI2];
JpI = inv(JI);

V = JpI/M*V;


% ah =size(4,16);
% for k=1:16
%     ah(1,k) =[V(1,k)];
%     ah(2,k) =[V(2,k)];
%     ah(3,k) =[0];
%     ah(4,k) =[0];
% end
% 
% ad =size(4,16);
% ad =V-ah;
% 
% Ans = size(4,16);
% counter = 1;
% for i = 1:16 %メインループ
%     for j = 1:16 %サブループ
%         if i~=j
%             a = norm(ad(:,i));
%             b = norm(ad(:,j));
%             c = b/(a+b);
%             d = V(:,i)-ah(:,i); 
%             e = V(:,j)-ah(:,i); 
%             if dot(d,e)<0
%                 Z =V(:,j)+((V(:,i)-V(:,j))*c);
%                 Ans(1,counter) = Z(1,1);
%                 Ans(2,counter) = Z(2,1);
%                 Ans(3,counter) = Z(3,1);
%                 Ans(4,counter) = Z(4,1);
%                 
%                 counter = counter+1;
%             end
%         end
%     end
% end 
%  Ans
%  Ans2=(Ans(1:2,:));
%  Vd=(V(1:2,:));
%  Ans3=Ans2.';
%  
%  
%  
%  p1=[0 0];
%  p2=[5 ;10.1];
%  p3=[9 8];
%  p4=[7 -2];
%  p5=[8 -12];
%  
%  
%  
%  plot(Vd(1,:),Vd(2,:),'mo')
%  quiver(p1(1), p1(2), p2(1)-p1(1), p2(2)-p1(2), 0)
%  [Ans4] = convhull(Ans3);
%  plot(Ans3(:,1),Ans3(:,2),'*')
%  hold on
%  plot(Ans3(Ans4,1),Ans3(Ans4,2))
%  quiver(p1(1), p1(2), p2(1)-p1(1), p2(2)-p1(2), 0, 'autoscalefactor', 1)
%  quiver(p1(1), p1(2), p3(1)-p1(1), p3(2)-p1(2), 0, 'autoscalefactor', 1)
%  quiver(p1(1), p1(2), p4(1)-p1(1), p4(2)-p1(2), 0, 'autoscalefactor', 1)
%  quiver(p1(1), p1(2), p5(1)-p1(1), p5(2)-p1(2), 0, 'autoscalefactor', 1)
%  xlim([-30 30]);
%  ylim([-30 30]);
%  
%             a = norm(ad(:,1));
%             b = norm(ad(:,5));
%             c = b/(a+b);
%             d = V(:,1)-ah(:,1); 
%             e = V(:,5)-ah(:,1);
%             f = size(5,1);
%             if dot(d,e)<0
%                 Z =V(:,5)+((V(:,1)-V(:,5))*c);
%                 f(1,1) = Z(1,1)
%                 f(2,1) = Z(2,1)
%                 f(3,1) = Z(3,1)
%                 f(4,1) = Z(4,1)
%             end
%             f
     
 ddq = Jp1*[ -20.52072 ; 5.41368]+JI2*[207.35217; 22.20507];
%  ddq1 = Jp1*[4.9 ;10.1]+JI2*[0; 0];
%  ddq2 = Jp1*[9 ;8]+JI2*[0;0];
%  ddq3 = Jp1*[7 ;-2]+JI2*[0;0];
%  ddq4 = Jp1*[8 ;-12]+JI2*[0;0];
  t =M*ddq;
%  t1 =M*ddq1
%  t2 =M*ddq2
%  t3 =M*ddq3
%  t4 =M*ddq4
 
% if draw == 1
%     FH = 1;
%     clf(FH)
%     fDraw_Graph3(FH,G,V)
% end

 end
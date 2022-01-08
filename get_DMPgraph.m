function[G,V] = get_DMPgraph(type,x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   DMPの頂点グラフをReturnする関数
%
%   [G,V] = get_DMPgraph(type, x)
%       G : グラフ
%       V : 頂点データ
%       x : 作業空間上の位置
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

draw = 0;

% %%
% clear all
% 
% x = [0.1;0.15;0];
% draw = 1
% type = '2d_RRR'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 引数の再配置
r = x(1:2);
Rx = x(3);

% 変数宣言
param = get_parameter(type);
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
        
        if temp == 1
            G(i,j) = 1;
            G(j,i) = 1;
        end
    end
end

% MFP頂点作成
q = fIKinematics(type, r,Rx);
[M,g] = get_matrix_minus(type, q);

for i = 1:size(V,2)
    V(:,i) = V(:,i)-g;
end

[J,Jp,U] = fJacobi(type, r,Rx);
Je = [J;U'];
V = Je/M*V;

% if draw == 1
%     FH = 1;
%     clf(FH)
%     fDraw_Graph3(FH,G,V)
% end

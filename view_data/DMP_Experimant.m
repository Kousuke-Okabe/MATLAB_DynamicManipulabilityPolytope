%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   DMP計測実験結果（まとめ）表示
%
%                                               '20.03.02 by.OKB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% データ読み込み
data = importdata('data.csv');

dz00 = [data.data(:,1),data.data(:,2)]';
dz10 = [data.data(:,3),data.data(:,4)]';

% パラメータ
FH = 1;
x00 = [0.2;0.1;0];
dx00 = [0;0;0];
x10 = [0.2;0.1;0];
dx10 = [0;0;10];

% 並進ベクトル
Vt00 = zeros(2,1);
Vt10 = zeros(2,1);
for i = 1:size(dz00,2)
    Vt00 = Vt00 + dz00(:,i);
    Vt10 = Vt10 + dz10(:,i);
end
Vt00 = Vt00/i;
Vt10 = Vt10/i-Vt00;

%% 描画
figure(FH)
clf(FH)
hold on
range = 4;
    % 実験結果描画
    plot(dz00(1,:),dz00(2,:),'b.')
    plot(Vt00(1),Vt00(2),'bo')
    
    plot(dz10(1,:),dz10(2,:),'r.')
    plot(Vt10(1)+Vt00(1),Vt10(2)+Vt00(2),'ro')
    quiver(0,0, Vt10(1),Vt10(2),'r', 'AutoScaleFactor',1.0)
    quiver(Vt00(1),Vt00(2), Vt10(1),Vt10(2),'r', 'AutoScaleFactor',1.0)
    
    % シミュレーション結果描画
    AH00 = fDraw_DMP(FH,0, 1,1, x00,dx00);
    set(AH00, 'Color','b', 'LineStyle','-')
    
    AH10 = fDraw_DMP(FH,0, 1,1, x10,dx10);
    set(AH10, 'Color','r', 'LineStyle','-')
    
    % 描画設定
    axis equal
    xlim([-1,1]*range)
    ylim([-1,1]*range)
    area = axis;
    plot(area(1:2),zeros(1,2),'k', zeros(1,2),area(3:4),'k')
    
    box on
    xlabel('$$\ddot{r}_{1}$$ [m/s$$^{2}$$]', 'Interpreter','latex')
    ylabel('$$\ddot{r}_{2}$$ [m/s$$^{2}$$]', 'Interpreter','latex')
%     legend({'Exp:$$\dot{z}=0$$[m/s]','Exp:$$\dot{z}=10$$[m/s]','Sim:$$\dot{z}=0$$[m/s]','Sim:$$\dot{z}=10$$[m/s]'}, 'Interpreter','latex')

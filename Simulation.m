clear all

param = get_parameter();
    Ts = param.Ts;
    Fs = param.Fs;
    Link = param.Link;

time = 0.1;
Ns = time*Fs;
t = (1:Ns)*Ts;

% Reference
r_ref = [ 0.2*ones(1,Ns);
          0.1*ones(1,Ns) ];
dz_ref = 50*ones(1,Ns);

% Valiable
q = nan(Link,Ns);
dq = nan(Link,Ns);

r = nan(Link-1,Ns);
dr = nan(Link-1,Ns);
dz = nan(1,Ns);

% Initialize
r(:,1) = [0.15;0.15];
dr(:,1) = zeros(Link-1,1);
q(:,1) = fIKinematics(r(:,1),0);
dq(:,1) = zeros(Link,1);
dz(:,1) = 0;

% Control Parameter
Kpr = diag([1,1])*50;
% Kpz = 0.8;
Gr = 1/500;
Gz = 1/10;

% Simulation
FH = 1;
figure(FH)
    clf(FH)
    fRoboAmination_Figure(FH)
for k = 1:Ns
    % Model
    q(:,k+1) = q(:,k) + dq(:,k)*Ts;
    
    % Controler
    [J,Jp,U,error] = fJacobi_q(q(:,k));
    r(:,k) = fKinematics(q(:,k));
    
    dr(:,k) = J*dq(:,k);
    dz(:,k) = U'*dq(:,k);
    
    dq(:,k+1) = Jp*Kpr*(r_ref(:,k)-r(:,k)) + U*dz_ref(:,k);
    
    % Translation Vector
    Vt = get_TranslateVector(q(:,k),[dr(:,k);dz(:,k)]);
    
    % Animation
    cla(FH)
    title(strcat('t=',num2str(k*Ts,3),', dz=',num2str(dz(:,k),3)));
    plot3(r_ref(1,k),r_ref(2,k),0,'ro')
    quiver3(r(1,k),r(2,k),0, Vt(1)*Gr,Vt(2)*Gr,Vt(3)*Gz,1,'b', 'ShowArrowHead','off')
    fRoboAnimation_Robo(FH,q(:,k), r(:,1:k));
end


figure(2)
clf(2)
subplot(2,1,1)
    hold on
    plot(t,r_ref(1,:),'r', t,r_ref(2,:),'r')
    plot(t,r(1,:),'b', t,r(2,:),'b')
    ylim([0,max(max(r))*1.2])
    xlabel('time [sec]')
    ylabel('q [rad]')

subplot(2,1,2)
    hold on
    plot(t,dz_ref(1,:),'r', t,dz(1,:),'b')
    ylim([0,max(dz)*1.2])
    xlabel('time [sec]')
    ylabel('dz [1/s]')







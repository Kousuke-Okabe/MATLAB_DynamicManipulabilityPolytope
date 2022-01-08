clear all

time = 1;
Fs = 1000;
Ts = 1/Fs;
N = time*Fs;

t = linspace(0,time-Ts,N);

param = get_parameter();
    Link = param.Link;

q = zeros(Link,N);
dq = zeros(Link,N);

x = nan(Link,N);
dx = nan(Link,N);

tau = nan(Link,N);
dq_ref = nan(Link,N);

% dq_ref(1,:) = pi/2*sin(2*pi*2*t);
% dq_ref(2,:) = pi/2*sin(2*pi*2*t);
% dq_ref(3,:) = pi/2*sin(2*pi*2*t);

% Reference
r_ref = nan(Link-1,N);
dz_ref = zeros(1,N);

% r_ref(1,1:N) = linspace(0.1,0.2,N);
r_ref(1,:) = 0.2*ones(1,N);
r_ref(2,:) = 0.1*ones(1,N);

dz_ref = 10*ones(1,N);

% Initial value
x(:,1) = [r_ref(:,1); 0];
% q(:,1) = fIKinematics(r_ref(:,1),0);
q(:,1) = [0;pi/2;pi/2];
dq(:,1) = zeros(Link,1);

% Gain
Kr = 1000*eye(2);

% Simulation
for k = 1:N
    [J,Jp,U] = fJacobi_q(q(:,k));
    [M,H,G] = get_matrix(q(:,k),dq(:,k));
    x(1:2,k) = fKinematics(q(:,k));
    
    Je = [J;U'];
    dx(:,k) = Je*dq(:,k);
    
    % Controler
    dq_ref(:,k) = Jp*Kr*(r_ref(:,k)-x(1:2,k))+ U*dz_ref(:,k);
    
    % Velocity Simulation
    if k < N
        q(:,k+1) = q(:,k) + dq_ref(:,k)*Ts;
        dq(:,k+1) = dq_ref(:,k);
    end

%     % Controller
%     tau(:,k) = M*Kq*(dq_ref(:,k)-dq(:,k));
% %     tau(:,k) = Je'*Kp*[r_ref(:,k)-x(1:2,k);dz_ref(:,k)];
% 
%     % Dynamics Simulation
%     if k < N
%     A = [eye(Link),Ts*eye(Link); zeros(Link),eye(Link)];
%     B = [Ts^2/2*inv(M); inv(M)*Ts];
%     X = A*[q(:,k);dq(:,k+1)] + B*(tau(:,k)-H);
%     
%     q(:,k+1) = X(1:Link);
%     dq(:,k+1) = X(Link+1:2*Link);
%     end
end

% View
FH = 2;
figure(FH)
clf(FH)
subplot(2,1,1)
hold on
    plot(t,r_ref(1,:),'r', t,x(1,:),'b')
    plot(t,r_ref(2,:),'r', t,x(2,:),'g')
    ylabel('Position[m]')
subplot(2,1,2)
    plot(t,dz_ref(1,:),'r', t,dx(3,:),'b')
    xlabel('time[s]')
    ylabel('dz[1/s]')

% FH = 3;
% figure(3)
% clf(FH)
% hold on
%     plot(t,dq_ref(1,:),'k', t,dq_ref(2,:),'k', t,dq_ref(3,:),'k')
%     plot(t,dq(1,:),'r', t,dq(2,:),'g', t,dq(3,:),'b')
%     xlabel('time[s]')
%     ylabel('Joint angle[rad]')

%% Animation
FH = 1;
figure(FH)
clf(FH)
    fRoboAnimation(FH,q,r_ref,Ts);
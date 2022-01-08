%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ó—é¿å±ÉfÅ[É^ï`âÊ
%
%                                                      '20.03.02 by OKB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

data = importdata('tessxx0.csv');

q = data(:,15:17)';
dq = data(:,18:20)';
ddr = data(:,30:32)';
tau = data(:,21:23)';

N = size(q,2);
t = (1:N)/1000;

r = nan(2,N);
dr = nan(2,N);
dr_norm = nan(1,N);
dz = nan(1,N);
ddr_norm = nan(1,N);
tau_norm = nan(1,N);

for k = 2:N
    r(:,k) = fKinematics(q(:,k));
    dq(:,k) = (q(:,k)-q(:,k-1))*1000;
    
    [J,Jp,U] = fJacobi(r(:,k),sum(q(:,k)));
    dr(:,k) = J*dq(:,k);
    dr_norm(:,k) = norm(dr(:,k));
    dz(:,k) = U'*dq(:,k);
    
    ddr_norm(1,k) = norm(ddr(:,k));
    tau_norm(1,k) = norm(tau(:,k));
end
%%
area = 600:1000;
figure(1)
clf(1)
subplot(4,1,1)
    plot(t(1,area),dr_norm(1,area))
    grid on
subplot(4,1,2)
    plot(t(1,area),dz(1,area))
    grid on
subplot(4,1,3)
    plot(t(1,area),ddr_norm(1,area))
    grid on
subplot(4,1,4)
    plot(t(1,area),tau_norm(1,area))
    grid on

clf(2)
fRoboAnimation(2,q(:,area), r(:,area),0)

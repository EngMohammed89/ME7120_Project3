% Newmark method
load K_M.mat % Load M and K matrieces
[Kr, Mr, C, wmax] = find_C(K,M); 
C=0
Inv=eye(150,150); % Identity matrix will be used to take inverse of [150x150] matrix

% % Average acceleration
% beta = 1/4; % beta
% gamma = 1/2; % gamma

% % Linear acceleration
% beta = 1/6; % beta
% gamma = 1/2; % gamma

% % Fox-Goodwin
% beta = 1/12; % beta
% gamma = 1/2; % gamma

% Algorithmically damped
gamma = 0.6; % gamma
beta = (1/4)*(gamma+0.5)^2; % beta

% % Hilber-hughes-Taylor (alpha-method)
% alpha =-0.2; % -1/3 <= alpha <=0
% beta = (1/4)*(1-alpha)^2; % beta
% gamma = 0.7; % gamma

% % Finding dt using omega critical. use this for Linear acceleration and
% % Fox-Goodwin methods.
% Z=0.02; % Damping ratio
% Ocrit=(Z*(gamma-0.5)+sqrt((gamma/2)-beta+(Z^2)*(gamma-0.5)^2))/(gamma/2-beta);
% dt=Ocrit/wmax;

dt = 0.0001; % delta t. use this for above methods (1st, 4th and 5th).
tf = 0.15; % Final t
n=floor(tf/dt); % Steps
t=zeros(n,1); % Time

D=zeros(150,n); % Displacement
DD=zeros(150,n); % Velocity
DDD=zeros(150,n); % Acceleration

% Initial conditions for Disp., vel., and acc.
R0=zeros(150,1); % Force, Rt
R0(149,1) = 100000;
D(:,1) = zeros(150,1);
DD(:,1) = zeros(150,1);
DDD(:,1) = Mr\R0;

Rt=zeros(150,1); % Force, Rt

for i=1:n

  % Impulse loading applied at node 51 starts form t=0 to t=0.01s.       
      if t(i)<=0.01
            Rt(149,1) = 100000;
      else
            Rt(149,1) = 0;
      end
      
   D(:,i+1) = (((1/(beta*dt^2))*Mr +(gamma/(beta*dt))*C+ Kr )\Inv )*(Rt +...
            Mr*( (1/(beta*dt^2))*D(:,i) + (1/(beta*dt))*DD(:,i) + (1/(2*beta)-1)*DDD(:,i) )...
            + C*((gamma/(beta*dt))*D(:,i)+(gamma/beta-1)*DD(:,i)+(gamma/beta-2)*(dt/2)*DDD(:,i)));
        
   DDD(:,i+1) = (1/(beta*dt^2))*( D(:,i+1) - D(:,i) - dt*DD(:,i) ) - (1/(2*beta) - 1)*DDD(:,i);
        
   DD(:,i+1) = (gamma/(beta*dt))*(D(:,i+1) - D(:,i)) - (gamma/beta - 1)*DD(:,i) -...
            dt*(gamma/(2*beta) - 1)*DDD(:,i);     
           
  t(i+1)=t(i)+dt;      
end

figure(1)
grid on; hold on
plot(t,D(121,:)) % Displacement of the rotational DOF at node 41
ylabel('$\theta z41$ (rad/s)','interpreter','latex')
xlabel('t (s)','interpreter','latex')
title('Hilber-hughes-Taylor (\alpha-method)')
figure(2)
grid on; hold on
plot(t,DD(121,:)) % Velocity of the rotational DOF at node 41
ylabel('$\dot{\theta}z41$ (rad/s)','interpreter','latex')
xlabel('t (s)','interpreter','latex')
title('Hilber-hughes-Taylor (\alpha-method)')
figure(3)
grid on; hold on
plot(t,DDD(121,:)) % Acceleration of the rotational DOF at node 41
ylabel('$\ddot{\theta}z41$ (rad/s)','interpreter','latex')
xlabel('t (s)','interpreter','latex')
title('Hilber-hughes-Taylor (\alpha-method)')

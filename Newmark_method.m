% Newmark method
load K_M_Newmark.mat
load C.mat

Inv=eye(150,150); % Identity matrix will be used to take inverse of [150x150] matrix
beta = 1/4; % beta
gamma = 1/2; % gamma
dt = 0.0001; % delta t
tf = 0.13; % Final t
t=0:dt:tf; % Time
n=tf/dt;

D=zeros(150,n);
DD=zeros(150,n);
DDD=zeros(150,n);

% Initial conditions
D(:,1) = zeros(150,1);
DD(:,1) = zeros(150,1);
DDD(:,1) = zeros(150,1);

Rt=zeros(150,1); % Force, R(t)



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
           
        
end

figure(1)
plot(t,D(121,:)) % Displacement of the rotational DOF at node 41
ylabel('$\theta z41$ (rad/s)','interpreter','latex')
xlabel('t (s)','interpreter','latex')
figure(2)
plot(t,DD(121,:)) % Velocity of the rotational DOF at node 41
ylabel('$\dot{\theta}z41$ (rad/s)','interpreter','latex')
xlabel('t (s)','interpreter','latex')
figure(3)
plot(t,DDD(121,:)) % acceleration of the rotational DOF at node 41
ylabel('$\ddot{\theta}z41$ (rad/s)','interpreter','latex')
xlabel('t (s)','interpreter','latex')

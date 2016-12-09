clear all
close all
load Project3_K_M.mat

NumNode = 51;       % Number of Nodes

%% K and M in 3D (306*306)
K_3D = full(K);      % Global Stiffness Matrix
M_3D = full(M);      % Global Mass Matrix

%% Reduce K and M from 3D to 2D
dof_2D = [];    % Dof Number in 2D case
for i = 1: NumNode
    dof_2D = [dof_2D, [1,2,6]+6*(i-1)];    
end
dof_2D;
% K and M in 2D (153*153)
K_2D = K_3D(dof_2D,dof_2D);
M_2D = M_3D(dof_2D,dof_2D);

%% Apply BCs (DOF1 = DOF2 = DOF152 = 0)
dof_BC = [1,2,152];
dof_all = 1:length(K_2D);
for i = 1:length(dof_BC)
    index = find(dof_all == dof_BC(i));
    dof_all(index) = [];
end
% Reduced K and M
Kr = K_2D(dof_all,dof_all);
Mr = M_2D(dof_all,dof_all);

%% Calculate Phi
L = chol(Mr)';
K_h = inv(L)*Kr*inv(L');
[vectors_V, values]=eig(K_h);

% Sort the eigenvalues and eigenvecgtors in ascending order
[values, index] = sort(diag(values));
vectors_V = vectors_V(:,index);

w = sqrt(values); % natural frequency

Phi = (L')\vectors_V;

%% Calculate C
zeta = 0.02;    % damping ratio
C = Phi'\diag(2*zeta*w)*inv(Phi);

% Testing Phi
% Phi'*Mr*Phi
% Phi'*Kr*Phi
% Phi'*C*Phi

%% Newmark beta method

for caseID = 1:3
    if caseID==1
        % Average acceleration
        beta = 1/4;
        gama = 1/2;
    elseif caseID==2
        % Algorithmically damped
        gama = 1/2+0.1;
        beta = 1/4*(gama+1/2)^2+0.1;
    elseif caseID==3
        % Hilber-Hughes-Taylor
        alpha = -0.1;
        beta = 1/4*(1-alpha)^2;
        gama = 1/2*(1-2*alpha);
    elseif caseID==4
        % Linear acceleration
        beta = 1/6;
        gama = 1/2;
    elseif caseID==5
        % Fox-Goodwin
        beta = 1/12;
        gama = 1/2;
    end

%     Omega_crit = (zeta*(gama-1/2)+sqrt(gama/2-beta+zeta^2*(gama-1/2)^2))/(gama/2-beta);
%     dt = Omega_crit/max(w);
%     
%     if isnan(dt)
%         dt = 0.0001;
%     end
    dt = 0.0001;
    t = 0.15;
    step = floor(t/dt);

    D = zeros(length(Kr),1);    % Displacement
    Dd = D;                     % Velocity
    Ddd = D;                    % Acceleration

    R = zeros(length(Kr),1);
    time = [0];
    for i = 1:step
        ti = dt*i;
        if ti<=0.01
            R(end-1) = 100000;  % Load vector
        else
            R = zeros(length(Kr),1);
        end

        Di = (1/(beta*dt^2)*Mr + gama/(beta*dt)*C +Kr)\...
            (R+Mr*(1/(beta*dt^2)*D(:,i) + 1/(beta*dt)*Dd(:,i) + (1/(2*beta)-1)*Ddd(:,i)) ...
            + C*( gama/(beta*dt)*D(:,i) + (gama/beta-1)*Dd(:,i) + (gama/beta-2)*dt/2*Ddd(:,i)));
        Ddi = gama/(beta*dt)*(Di-D(:,i))-(gama/beta-1)*Dd(:,i)-dt*(gama/(2*beta)-1)*Ddd(:,i);
        Dddi = 1/(beta*dt^2)*(Di-D(:,i)-dt*Dd(:,i))-(1/(2*beta)-1)*Ddd(:,i);
        D = [D,Di];
        Dd = [Dd,Ddi];
        Ddd = [Ddd,Dddi];
        time = [time,ti];

    end
    D;
    Dd;
    Ddd;

    % plotting
    plot(time,D(121,:),'LineWidth',3);
    xlabel('t(s)');
    ylabel('\theta_{z41}(rad)');
    title('Modal analysis  \zeta=0.02')
    hold on
end









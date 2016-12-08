function [Kr Mr C] = find_C(K,M)

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
zeta = 0.02;
C = Phi'\diag(2*zeta*w)*inv(Phi);
end
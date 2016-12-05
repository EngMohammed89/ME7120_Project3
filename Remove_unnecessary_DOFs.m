% This code removes the unnecessary DOFs.
% For node 1, there is only 1 DOF which is rotation about z-axis.
% For node 51, there are 2 DOFs which are displacement in x-axis 
% and rotation about z-axis.
% For other nodes, there are 3 DOFs which are displacement in x-axis, 
% displacement in y-axis and rotation about z-axis. 
% The size of M and K will be reduced from [306x306] to [150x150].

load frame_structure.mat
clearvars -except K M % Clear unnecessary variables
K=K(1:306,1:306); % [306x306] Stiffness matrix
M=M(1:306,1:306); % [306x306] mass matrix

c1=7:6:301; % DOFs: displacements in x-axis
c2=8:6:296; % DOFs: displacements in y-axis
c3=6:6:306; % DOFs: rotaions about z-axis
columns=sort([c1 c2 c3]); % 150 DOFs
rows=columns;
% Deleting the unnecessary DOFs
Kr=full(K(rows,columns)); % [150x150] Stiffness matrix
Mr=full(M(rows,columns)); % [150x150] Mass matrix

% Saving M and K to *.mat file. 
save('K_M_Newmark.mat','Kr','Mr')

%%%%%%%%%%%%%%%%%%%%
% PRADMM for cube data 
%%%%%%%%%%%%%%%%%%%%
clc;
clear;
run('addpath_myslam.m');
%% Settings
Setting.g2omat_information.kappa = 0.03;
Setting.g2omat_information.tau = 0.12;

%% read .mat file
solution_index = "09_6";
cube_filename = "cube_"+ solution_index + ".mat";
load(cube_filename)
%% Deal with .g2o file
vertex_true=[vertex_true(:,1:4),vertex_true(:,8),vertex_true(:,5:7)];
edge_true=[edge_true(:,1:5),edge_true(:,9),edge_true(:,6:8),edge_true(:,10:30)];
vertex_true(:,1)=vertex_true(:,1)+1;
edge_true(:,1:2)=edge_true(:,1:2)+1;

vertex=[vertex_raw(:,1:4),vertex_raw(:,8),vertex_raw(:,5:7)];
edge=[edge_raw(:,1:5),edge_raw(:,9),edge_raw(:,6:8),edge_raw(:,10:30)];
vertex(:,1)=vertex(:,1)+1;
edge(:,1:2)=edge(:,1:2)+1;

num_edge = size(edge,1);
edge(:,10:30) = ones(num_edge,1)*[Setting.g2omat_information.tau,0,0,0,0,0,...
                                       Setting.g2omat_information.tau,0,0,0,0,...
                                       Setting.g2omat_information.tau,0,0,0,...
                                       Setting.g2omat_information.kappa,0,0,...
                                       Setting.g2omat_information.kappa,0,...
                                       Setting.g2omat_information.kappa];

%% data storage 
fprintf('Data Storage Type:  Adjacency matrix ...');
Ad_Mat = Adjacency_Mat_sp(vertex,edge); 
fprintf('done.\n\n');

%% Initial
disp('==================');
disp('Initialization starting...');
fprintf('\t\t Chordal + LeastSquare\n');
num_v = size(vertex,1); 
num_edge = size(edge,1);      
[opts.p0,tchordal] = initialization_chordal_ve(edge);
opts.q0=opts.p0;
if exist('tchordal', 'var')
    opts.t0 = tchordal;
else
    opts.t0 = vertex(:,2:4)';
end
opts.lambda0 = zeros(4,num_v);
opts.u0 = opts.t0; 
opts.z0 = zeros(3,num_v);

%% parameter
opts.beta1 = 10;  
opts.beta2 = 0.1;
opts.H1=100;   
opts.H2=0.001;
opts.H3=0.001;
opts.H4=0.001;
opts.MaxIter = 210;
opts.tol = 1e-3;
opts.max_beta1 = 100;
opts.max_beta2 = 100;
opts.rho = 1.01;
[opts.kappa,opts.tau] = extract_information(edge);
[opts.Sigma1,opts.Sigma2]=generate_sigma(opts.tau,opts.kappa);
opts.vertex_true=vertex_true;
opts.edge_true=edge_true;
fprintf('Finished.\n\n');

%% ground true and noised trajectory
figure(1)
subplot(1,2,1)
plot_trajectory_line(vertex_true(:,2:4))
title('Ground True')
subplot(1,2,2)
plot_trajectory_line(vertex(:,2:4))
title('Noised Trajectory')

%% PRADMM
disp('==================');
disp('PRADMM starting...');
result = PGO_PRADMM_cube(vertex,edge,Ad_Mat,opts); 
result.NRMSE = result.NRMSE/sqrt(num_v);
fprintf('Finished.\n\n');

%% result
figure(3)     
plot_trajectory_line(result.pose7n_new(:,1:3))
title('Recovered Trajectory')
axis equal
fprintf('RE: %12.3e,  NRMSE:%12.3e \n', result.RE(end),result.NRMSE(end));

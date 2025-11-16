%%%%%%%%%%%%%%%%%%%%
% PRADMM for benchmark data in SLAM
%%%%%%%%%%%%%%%%%%%%

clc;
clear;
run('addpath_myslam.m');


%% Read in .g2o file
g2o_file='parking-garage.g2o';
[vertex_raw,edge_raw] = parse_g2o(g2o_file);

vertex=[vertex_raw(:,1:4),vertex_raw(:,8),vertex_raw(:,5:7)];
edge=[edge_raw(:,1:5),edge_raw(:,9),edge_raw(:,6:8),edge_raw(:,10:30)];
vertex(:,1)=vertex(:,1)+1;
edge(:,1:2)=edge(:,1:2)+1;

Ad_Mat = Adjacency_Mat_sp(vertex,edge);
fprintf('done.\n\n');
[opts.kappa,opts.tau] = extract_information(edge);

%% initial
disp('==================');
disp('Initialization starting...');
fprintf('\t\t Chordal + LeastSquare\n');
num_v = size(vertex,1); 
num_edge = size(edge,1); 
[opts.p0,tchordal] = initialization_chordal(load_g2o_data(g2o_file));
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
opts.beta1 = 2;  
opts.beta2 = 2;
opts.MaxIter = 40;
opts.tol = 1e-4;
opts.max_beta1 = 100;
opts.max_beta2 = 100;
opts.rho = 1;
opts.H1=5;   
opts.H2=5;
opts.H3=0.1;
opts.H4=0.1;
[opts.Sigma1,opts.Sigma2]=generate_sigma(opts.tau/2000,opts.kappa/5000); 
fprintf('Finished.\n\n');

%% noised trajectory
figure(1)
plot_trajectory(vertex(:,2:4))
title('Noised Trajectory')

disp('==================');
disp('PRADMM starting...');
result = PGO_PRADMM(vertex,edge,Ad_Mat,opts);
fprintf('Finished.\n\n');

%% recovered trajectory
figure(2)  
plot_trajectory(result.pose7n_new(:,1:3))
title('Recovered Trajectory')
% axis equal

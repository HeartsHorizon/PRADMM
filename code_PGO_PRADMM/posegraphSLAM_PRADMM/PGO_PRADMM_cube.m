%% PRADMM for pose graph optimization

function result = PGO_PRADMM_cube(vertex,edge,Ad_Mat,opts)

%% pose graph optimization
                
check_tau = 1.4;

% check_p_subproblem
       %  1  -  Exact Solution        
       %  2  -  Generalized Eigenvalue
check_p_subproblem = 1;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
%       motion         -    the known motion
%       opts    -    Structure value in Matlab. The fields are
%           opts.p0         -   starting ponit
%           opts.q0         -   starting ponit
%           opts.t0         -   starting ponit
%           opts.u0         -   starting ponit
%           opts.lambda0    -   starting ponit
%           opts.z0         -   starting ponit
%           opts.tol        -   termination tolerance
%           opts.MaxIter    -   maximum number of iterations
%           opts.beta1       -   stepsize for dual variable updating in ADMM
%           opts.max_beta1   -   maximum stepsize
%           opts.beta2       -   stepsize for dual variable updating in ADMM
%           opts.max_beta2   -   maximum stepsize
%           opts.rho        -   rho>=1, ratio used to increase beta1
%           opts.Sigma1     -   the weighting coefficient in the objective function
%           opts.Sigma2     -   the weighting coefficient in the objective function
%           opts.H1         -   the prox coefficient in the p-subproblem
%           opts.H2         -   the prox coefficient in the q-subproblem
%           opts.H3         -   the prox coefficient in the t-subproblem
%           opts.H4         -   the prox coefficient in the t-subproblem
% 
% Output:
%       pose7n_new     -    the pose from origin to here
%       k              -    the number of iterations
%       epsm           -    MaxIter*1 vector of iterations residual
%       t_solve        -    running time 
%       obj_admm       -    value of ADMM
%       obj_Lag        -    value of Lag
%       RE             -    relative error
%       NRMSE          -    Normalized Root Mean Square Error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% parameter
num_v = size(vertex,1); % 顶点数

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');         tol = opts.tol;                end
if isfield(opts, 'MaxIter');     MaxIter = opts.MaxIter;        end
if isfield(opts, 'beta1');        beta1 = opts.beta1;              end
if isfield(opts, 'max_beta1');    max_beta1 = opts.max_beta1;      end
if isfield(opts, 'beta2');        beta2 = opts.beta2;              end
if isfield(opts, 'max_beta2');    max_beta2 = opts.max_beta2;      end
if isfield(opts, 'rho');         rho = opts.rho;                end
if isfield(opts, 'p0');          p0 = opts.p0;                  end
if isfield(opts, 'q0');          q0 = opts.q0;                  end
if isfield(opts, 't0');          t0 = opts.t0;                  end
if isfield(opts, 'u0');          u0 = opts.u0;                  end
if isfield(opts, 'lambda0');     lambda0 = opts.lambda0;        end
if isfield(opts, 'z0');          z0 = opts.z0;                  end
if isfield(opts, 'Sigma1');      Sigma1 = opts.Sigma1;          end
if isfield(opts, 'Sigma2');      Sigma2 = opts.Sigma2;          end
if isfield(opts, 'H1');          H1 = opts.H1;                  end
if isfield(opts, 'H2');          H2 = opts.H2;                  end
if isfield(opts, 'H3');          H3 = opts.H3;                  end
if isfield(opts, 'H4');          H4 = opts.H4;                  end
vertex_true=opts.vertex_true;

%% Initialization
% output initialization
RE = zeros(MaxIter,1);
NRMSE = zeros(MaxIter,1);
t_solve = zeros(MaxIter,1);
stopc = 1;  
k=0;


% point initialization
p = p0;
q = q0;
t = t0;
u = u0;
lambda = lambda0;
z = z0;

% time initialization
times.ppp=zeros(num_v-1,1);
times.qqq=zeros(num_v-1,1);
times.ttt=zeros(num_v-1,1);
times.uuu=zeros(num_v-1,1);

t_start = tic;
%% Update
while(stopc > tol && k <= MaxIter)
    k = k+1;
    p0 = p;     
    q0 = q;     
    t0 = t;
    u0 = u;
    lambda0 = lambda;
    z0 = z;  
    
    %% update p^k+1
    if check_p_subproblem == 1  
        for i = 2 : num_v
            times.tttp=tic;
            [row1,row2]=Adjacency_Mat_row(Ad_Mat,i);   
            [col1,col2]=Adjacency_Mat_col(Ad_Mat,i);   

            p_s=zeros(4,1);
            num_col=length(col1);
            num_row=length(row1);
            for j=1:num_row
                M1=matrixM(q0(:,i))*matrixM([0,edge(row2(j),3:5)]).*[1,-1,-1,-1];
                s1=[0;t0(:,row1(j))-u0(:,i)];                
                p_s = p_s + M1'*Sigma1{row2(j)}*s1;
            end
            for l=1:num_col
                M2=matrixW(edge(col2(l),6:9))*matrixW(q0(:,col1(l))).*[1,-1,-1,-1];
                s2=[1,0,0,0]';
                p_s = p_s + M2'*Sigma2{col2(l)}*s2;
            end
            s3 = q0(:,i)+1/beta1*lambda0(:,i);
            p_s = -2*(p_s + beta1/2*s3+1/2*H1*p0(:,i));
            p(:,i) = -p_s/norm(p_s);
            times.ppp(i-1)=times.ppp(i-1)+toc(times.tttp);
        end
    else
        for i = 2 : num_v
            times.tttp=tic;
            [row1,row2]=Adjacency_Mat_row(Ad_Mat,i);   
            [col1,col2]=Adjacency_Mat_col(Ad_Mat,i);   

            p_M=zeros(4,4);
            p_s=zeros(4,1);
            num_col=length(col1);
            num_row=length(row1);
            for j=1:num_row
                M1=matrixM(q0(:,i))*matrixM([0,edge(row2(j),3:5)]).*[1,-1,-1,-1];
                s1=[0;t0(:,row1(j))-u0(:,i)]; 
                p_M = p_M + M1'*Sigma1{row2(j)}*M1;
                p_s = p_s + M1'*Sigma1{row2(j)}*s1;
            end
            for l=1:num_col
                M2=matrixW(edge(col2(l),6:9))*matrixW(q0(:,col1(l))).*[1,-1,-1,-1];
                s2=[1,0,0,0]';
                p_M = p_M + M2'*Sigma2{col2(l)}*M2;
                p_s = p_s + M2'*Sigma2{col2(l)}*s2;
            end
            s3 = q0(:,i)+1/beta1*lambda0(:,i);
            p_M = 2*(p_M + beta1/2*eye(4)+1/2*H1*eye(4));
            p_s = -2*(p_s + beta1/2*s3+1/2*H1*p0(:,i));
            
            p(:,i) = TRS_gen_eig(p_M,p_s);   
            times.ppp(i-1)=times.ppp(i-1)+toc(times.tttp);
        end
    end

    % update q^k+1
    for i = 2 : num_v
        times.tttq = tic;
        [row1,row2]=Adjacency_Mat_row(Ad_Mat,i);   

        q_W=zeros(4,4);
        q_u=zeros(4,1);
        num_row=length(row1);
        for j=1:num_row
            W1 = matrixW(p(:,i))'*matrixW([0,edge(row2(j),3:5)]);
            u1 = [0;t0(:,row1(j))-u0(:,i)];  
            W2 = matrixW(edge(row2(j),6:9))*matrixM(p(:,row1(j)))';
            u2 = [1,0,0,0]';

            q_W = q_W + W1'*Sigma1{row2(j)}*W1 + W2'*Sigma2{row2(j)}*W2;
            q_u = q_u + W1'*Sigma1{row2(j)}*u1 + W2'*Sigma2{row2(j)}*u2;
        end
        u3 = p(:,i)-1/beta1*lambda0(:,i);
        q_W = q_W + beta1/2*eye(4) + 1/2*H2*eye(4);
        q_u = q_u + beta1/2*u3 + 1/2*H2*q0(:,i);

        q(:,i) = eye(4)/q_W*q_u;
        times.qqq(i-1)=times.qqq(i-1)+toc(times.tttq);
    end


    % update t^k+1
    for i = 2 : num_v
        times.tttt = tic;
        [col1,col2]=Adjacency_Mat_col(Ad_Mat,i);   
        p_M=zeros(3,3);
        p_s=zeros(3,1);
        num_col=length(col1);
        for l=1:num_col
            s1 = u0(:,col1(l))+Q_to_Vector(matrixW(p(:,col1(l)))'*matrixM(q(:,col1(l)))*...
                [0,edge(col2(l),3:5)]' )';
            Sigma1_hat = Sigma1{col2(l)}(2:4,2:4);     
            p_M = p_M + Sigma1_hat;
            p_s = p_s + Sigma1_hat*s1;
        end
        s2 = u0(:,i) + 1/beta2*z0(:,i);
        p_M = p_M + beta2/2*eye(3)+1/2*H3*eye(3);
        p_s = p_s + beta2/2*s2 + 1/2*H3*t0(:,i);
        t(:,i) = p_M\p_s;
        times.ttt(i-1)=times.ttt(i-1)+toc(times.tttt);
    end
    
    % update u^k+1
    for i = 2 : num_v
        times.tttu = tic;
        [row1,row2]=Adjacency_Mat_row(Ad_Mat,i);  
        p_M=zeros(3,3);
        p_s=zeros(3,1);
        num_row=length(row1);
        for j=1:num_row
            s1 = t(:,row1(j))-Q_to_Vector(matrixW(p(:,i))'*matrixM(q(:,i))*...
                [0,edge(row2(j),3:5)]')';
            Sigma1_hat = Sigma1{row2(j)}(2:4,2:4);     
            p_M = p_M + Sigma1_hat;
            p_s = p_s + Sigma1_hat*s1;
        end
        s2 = t(:,i) - 1/beta2*z0(:,i);
        p_M = p_M + beta2/2*eye(3)+1/2*H4*eye(3);
        p_s = p_s + beta2/2*s2 + 1/2*H4*u0(:,i);
        u(:,i) = p_M\p_s;
        times.uuu(i-1)=times.uuu(i-1)+toc(times.tttu);
    end       

    % update lambda^k+1
    lambda = lambda0 - check_tau*beta1*(p-q);
    z = z0 - check_tau*beta2*(t-u);

    % update stepsize
    beta1 = min(rho*beta1,max_beta1);
    beta2 = min(rho*beta2,max_beta2);

    % update residuals
    ex.q=q-q0;
    ex.t=t-t0;
    ex.u=u-u0;
    ex.lambda=lambda-lambda0;
    ex.z=z-z0;
    eps.pri = norm(ex.lambda,'fro')/beta1 + norm(ex.z,'fro')/beta2;
    eps.dual = norm(ex.q,'fro')*beta1 + norm(ex.t,'fro')*beta2 + norm(ex.u,'fro')*beta2;
    stopc = eps.pri+eps.dual;
    
    % update CPU time
    t_solve(k) = toc(t_start)-sum(times.ppp)-sum(times.qqq)-sum(times.ttt)-sum(times.uuu)...
    +max(times.ttt)+max(times.uuu)+max(times.ppp)+max(times.qqq);

    % update RE, NRMSE
    RMSE = norm_quater(p',vertex_true(:,5:8))+norm(t'-vertex_true(:,2:4),'fro');
    RE(k) = RMSE/norm(vertex_true(:,2:8),'fro');
    NRMSE(k) = RMSE /(max(max(vertex_true(:,2:4)))-min(min(vertex_true(:,2:4))));


end
t_solve = t_solve(1:k);
RE = RE(1:k);
NRMSE = NRMSE(1:k);
pose7n_new = [t',p'];

%%
result.pose7n_new=pose7n_new;
result.k=k;
result.t_solve=t_solve;
result.RE=RE;
result.NRMSE=NRMSE;
end
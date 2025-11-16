function [Sigma1,Sigma2]=generate_sigma_nomatrix(tau,kappa)
m = length(tau);
Sigma1 = cell(m,1);
Sigma2 = cell(m,1);
for edge_id = 1:m
    Sigma1{edge_id} = tau(edge_id);
    Sigma2{edge_id} = kappa(edge_id);
end
end
function [Sigma1,Sigma2]=generate_sigma(tau,kappa)
m = length(tau);
Sigma1 = cell(m,1);
Sigma2 = cell(m,1);
for edge_id = 1:m
    tau_i = tau(edge_id);
    kappa_i = kappa(edge_id);
    Sigma1{edge_id} = diag([tau_i,tau_i,tau_i,tau_i]);
    Sigma2{edge_id} = diag([kappa_i,kappa_i,kappa_i,kappa_i]);
end
end
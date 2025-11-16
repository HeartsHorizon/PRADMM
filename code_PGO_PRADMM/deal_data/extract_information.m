function [kappa,tau] = extract_information(edge)
m=size(edge,1);
kappa = zeros(m,1);
tau = zeros(m,1);
const = 1;
for edge_id=1:m
    v = edge(edge_id,10:30);
    I=[v(1:6);
        0,v(7:11);
        0,0,v(12:15);
        0,0,0,v(16:18);
        0,0,0,0,v(19:20);
        0,0,0,0,0,v(21)];
    I=triu(I,1)'+I;
    tau(edge_id) = 3 / trace(inv(I(1:3, 1:3)) );
    kappa(edge_id) = 3 / ( const*trace(inv(I(4:6, 4:6))));  
end

end
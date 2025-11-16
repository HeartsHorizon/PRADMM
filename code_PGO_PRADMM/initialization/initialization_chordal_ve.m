function [p0, tchordal] = initialization_chordal_ve(e)
edges = e(:,1:2);
m = size(edges,1);
R = cell(m,1);
t = cell(m,1);
kappa = cell(m,1);
tau = cell(m,1);
for edge_id = 1:m
    R{edge_id} = Q_to_SO3(e(edge_id,6:9));
    t{edge_id} = e(edge_id,3:5)';
    v = e(edge_id,10:30);
    I=[v(1:6);
        0,v(7:11);
        0,0,v(12:15);
        0,0,0,v(16:18);
        0,0,0,0,v(19:20);
        0,0,0,0,0,v(21)];
    I=triu(I,1)'+I;
    tau{edge_id} = 3 / trace(inv(I(1:3, 1:3)) );
    kappa{edge_id} = 3 / (2 *trace(inv(I(4:6, 4:6))));
end


measurements.edges=edges;
measurements.R=R;
measurements.t=t;
measurements.kappa=kappa;
measurements.tau=tau;

[p0, tchordal] = initialization_chordal(measurements);

end
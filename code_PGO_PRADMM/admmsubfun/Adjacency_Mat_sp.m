function A=Adjacency_Mat_sp(v,e)
num_v=size(v,1);   
edge=e(:,1:2);
num_edge = size(edge,1); 
A=sparse(edge(:,1),edge(:,2),(1:1:num_edge),num_v,num_v);
end
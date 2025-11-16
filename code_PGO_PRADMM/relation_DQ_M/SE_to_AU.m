function q_vertex = SE_to_AU(SE_vertex)
R=SE_vertex.R;
t=SE_vertex.t;
num_v=size(t,1);
q_vertex=zeros(num_v,7);
q_vertex(:,1:3)=t;
for i=1:num_v
    q_vertex(i,4:7)=SO3_to_Q(R{i});
end
end
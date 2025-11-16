function x7n=DQ_to_AQ(q)
[m,n]=size(q);
if m>n
    q=q';
end

q_r = q(1:4);
q_d = q(5:8);
t = Q_to_Vector(2*multi_Q(q_d,Q_star(q_r)));
x7n = [t,q_r];
end
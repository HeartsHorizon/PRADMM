function q=AQ_to_DQ(x7n)

[m,n]=size(x7n);
if m>n
    x7n=x7n';
end

q_r = x7n(4:7);
t = x7n(1:3);
q_d=0.5*multi_Q([0,t],q_r);
q=[q_r,q_d];
end
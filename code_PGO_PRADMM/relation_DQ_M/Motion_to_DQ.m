function q=Motion_to_DQ(x)
[m,n]=size(x);
if m>n
    x=x';
end
r=x(1:3);
t=x(4:6);

q_r=axisangle_to_Q(r);
q_d=0.5*multi_Q([0,t],q_r);
q=[q_r,q_d];
end
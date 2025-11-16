function x7n=Motion_to_AQ(x)

[m,n]=size(x);
if m>n
    x=x';
end

r = x(1:3);
t = x(4:6);
x7n = [t,axisangle_to_Q(r)];
end
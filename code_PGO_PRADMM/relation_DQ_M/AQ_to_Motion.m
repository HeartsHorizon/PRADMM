function x=AQ_to_Motion(x7n)

[m,n]=size(x7n);
if m>n
    x7n=x7n';
end

r = Q_to_axisangle(x7n(4:7));
t = x7n(1:3);
x = [r,t];
end
function q_r=axisangle_to_Q(x)

[m,n]=size(x);
if m>n
    x=x';
end


if (norm(x)-0)<1e-10
    q_r=[1,0,0,0];
else
    theta=mod(norm(x,2),2*pi);
    q_r=[cos(theta/2),x/norm(x)*sin(theta/2)];
end
end
function R=axisangle_to_SO3(x) 

[m,n]=size(x);
if m<n
    x=x';
end

x=x/norm(x)*mod(norm(x),2*pi);

n=x/norm(x);
theta=mod(norm(x),2*pi);
R=cos(theta)*eye(3)+(1-cos(theta))*(n*n')+sin(theta)*[0,-n(3),n(2);n(3),0,-n(1);-n(2),n(1),0];
end
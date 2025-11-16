function T=se_to_SE3(xi)
rho=xi(1:3);
x=xi(4:6);
R=axisangle_to_SO3(x);

theta=norm(x); 
n=x/theta; 
J=sin(theta)/theta*eye(3)+(1-sin(theta)/theta)*(n*n')+(1-cos(theta))/theta*[0,-n(3),n(2);n(3),0,-n(1);-n(2),n(1),0];

t=J*rho;
T=[R,t;zeros(1,3),1];
end
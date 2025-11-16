function xi=SE3_to_se(T)

x=SO3_to_axisangle(T(1:3,1:3));
theta=norm(x);
n=x/theta; 
J_inv=theta/2*cot(theta/2)*eye(3)+(1-theta/2*cot(theta/2))*(n*n')-theta/2*[0,-n(3),n(2);n(3),0,-n(1);-n(2),n(1),0];
rho=J_inv*T(1:3,4);
xi=[rho;x];
end
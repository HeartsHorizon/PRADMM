function T_inv=SE3_inv(T)
R=T(1:3,1:3);
t=T(1:3,4);
T_inv=[R',-R'*t;
    zeros(1,3),1];
end
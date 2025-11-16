function R=Q_to_SO3(q)
s=q(1);
v=[q(2),q(3),q(4)]';
Kv=[0,-q(4),q(3);
    q(4),0,-q(2);
    -q(3),q(2),0];
R=v*v'+s^2*eye(3)+2*s*Kv+Kv*Kv;
end
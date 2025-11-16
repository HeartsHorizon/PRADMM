function x=Q_to_axisangle(q)
[m,~]=size(q);
if m>1
    q=q';
end
if abs(q(1)^2-1)<1e-10   
    x = [0,0,0];
else
    x = 2*acos(q(1))*q(2:4)/norm(q(2:4));
end
end
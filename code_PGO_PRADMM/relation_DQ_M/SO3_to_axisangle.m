function x=SO3_to_axisangle(R)

theta=acos((trace(R)-1)/2);
[v,e]=eig(R);
if norm(e(1,1)-1)<1e-5
    n=v(:,1);
elseif  norm(e(2,2)-1)<1e-5
    n=v(:,2);
else
    n=v(:,3);
end
x=n*theta;

if norm(R-eye(3))<1e-5
    x = [0,0,0]';
end

if norm(axisangle_to_SO3(x)-R,'fro')>1e-10
    x=-x;
end
end
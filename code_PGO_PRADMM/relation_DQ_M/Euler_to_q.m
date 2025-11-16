function q=Euler_to_q(a)
a1=a(1);
a2=a(2);
a3=a(3);
q=[cos(a3/2)*cos(a2/2)*cos(a1/2)+sin(a3/2)*sin(a2/2)*sin(a1/2),
    cos(a3/2)*cos(a2/2)*sin(a1/2)-sin(a3/2)*sin(a2/2)*cos(a1/2),
    cos(a3/2)*sin(a2/2)*cos(a1/2)+sin(a3/2)*cos(a2/2)*sin(a1/2),
    sin(a3/2)*cos(a2/2)*cos(a1/2)-cos(a3/2)*sin(a2/2)*sin(a1/2)];
end
function epsilon = norm_quater(p,q)
[m_p,n_p]=size(p);
[m_q,~]=size(q);
if m_p==4
    p=p';
    n=n_p;
else
    n=m_p;
end
if m_q==4
    q=q';
end
epsilon=0;
for i=1:n
    z=min(norm(p(i,:)-q(i,:)),norm(p(i,:)+q(i,:)));
    epsilon = epsilon+z;
end
end
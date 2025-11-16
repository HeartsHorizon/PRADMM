function x_star = TRS_gen_eig(A,b)
n = size(b,1);
B = eye(n);
M0 = [-B, A;A,-b*b'];
M1 = [zeros(n),B;B,zeros(n)];
lambda = max(real(eig(M0,-M1)));
x_star = -(A+lambda*B)\b;
end
function [X1,X2]=Adjacency_Mat_row(A,i)

B=A(i,:);
X1=find(B)';
X2=full(B(X1))';


end
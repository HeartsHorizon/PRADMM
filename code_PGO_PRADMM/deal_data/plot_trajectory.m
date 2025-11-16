function plot_trajectory(X)
x1=X(:,1);
x2=X(:,2);
x3=X(:,3);
scatter3(x1,x2,x3,10,linspace(1,10,length(x1)),'filled')
end
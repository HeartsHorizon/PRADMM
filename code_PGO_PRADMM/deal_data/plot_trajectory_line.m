function plot_trajectory_line(X)
x1=X(:,1);
x2=X(:,2);
x3=X(:,3);
scatter3(x1,x2,x3,[],linspace(1,10,length(x1)),'filled')
line(x1,x2,x3,'Color','b','LineWidth',1)
end
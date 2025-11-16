function c=multi_DQ(x,y)
x_r=x(1:4);
x_d=x(5:8);
y_r=y(1:4);
y_d=y(5:8);
c_r=multi_Q(x_r,y_r);
c_d=multi_Q(x_r,y_d)+multi_Q(x_d,y_r);
c=[c_r,c_d];
end
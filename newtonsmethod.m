fun_x=@(x1,x2)(100*(x2-(x1^2))^2)+(1-x1)^2;
 %Newton'sMethod
 x_initial=[-1;1];
 err1=100;
 err2=100;
 tol=1e-06;
 %Gradient=[diff(fun_x(x1,x2),x1);diff(fun_x(x1,x2),x2)]
 count=0;
 while err1>tol||err2>tol
 % err=abs(x_new-x_initial)/x_initial;
 Hessian=[-400*x_initial(2)+1200*x_initial(1)^2+2,-400*x_initial(1);-400*x_initial(1),200]
 inv(Hessian)
 Gradient=[2*x_initial(1)-400*x_initial(1)*(-x_initial(1)^2+ x_initial(2))-2;-200*x_initial(1)^2+200*x_initial(2)]
 x_new=x_initial-inv(Hessian)*Gradient
 err1=abs((x_new(1)-x_initial(1))/x_initial(1));
 err2=abs((x_new(2)-x_initial(2))/x_initial(2));
 x_initial=x_new;
 count=count+1;
 end
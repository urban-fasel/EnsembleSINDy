function dx = lorenzForcing(t,x,u,Beta)
dx = [
Beta(1)*(x(2)-x(1))+u;
x(1)*(Beta(3)-x(3))-x(2);
x(1)*x(2)-Beta(2)*x(3);
];
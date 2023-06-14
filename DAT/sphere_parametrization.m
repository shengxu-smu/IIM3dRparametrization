N=101;
ua=-1;
ub=1;
va=0;
vb=pi;
uh=(ub-ua)/(N-1);
vh=(vb-va)/(N-1);
u=[ua:uh:ub];
v=[va:vh:vb];
for j=1:N
   for i=1:N
      x(i,j)=2*u(i)/(u(i)*u(i)+v(j)*v(j)+1);
      y(i,j)=2*v(j)/(u(i)*u(i)+v(j)*v(j)+1);
      z(i,j)=(u(i)*u(i)+v(j)*v(j)-1)/(u(i)*u(i)+v(j)*v(j)+1);
   end
end
figure(1)
mesh(x,y,z)
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

for j=1:N
   for i=1:N
      x(i,j)=u(i);
      y(i,j)=sqrt(1-u(i)*u(i))*sin(v(j));
      z(i,j)=sqrt(1-u(i)*u(i))*cos(v(j));
   end
end
figure(2)
mesh(x,y,z)
axis equal
xlabel('x')
ylabel('y')
zlabel('z')


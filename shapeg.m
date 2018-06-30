function g=shapeg(xi,dx)
g(1,:)=[-(1.-xi(2))/dx(1),-(1.-xi(1))/dx(2)];
g(2,:)=[(1.-xi(2))/dx(1),-xi(1)/dx(2)];
g(3,:)=[xi(2)/dx(1),xi(1)/dx(2)];
g(4,:)=[-xi(2)/dx(1),(1.-xi(1))/dx(2)];
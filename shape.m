function n=shape(xi)
n(1,:)=(1-xi(1))*(1-xi(2));
n(2,:)=xi(1)*(1-xi(2));
n(3,:)=xi(1)*xi(2);
n(4,:)=(1-xi(1))*xi(2);
function N=getGrad(np,lx,ly)
s=np(1); t=np(2);
N(1,:)=0.25*[-(1-t),  (1-t), (1+t), -(1+t)]*2/lx;
N(2,:)=0.25*[-(1-s), -(1+s), (1+s),  (1-s)]*2/ly;



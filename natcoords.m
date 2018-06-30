function [xpn,nep]=natcoords(xp,dx,xmin,ne)
np=size(xp,2);
for ip=1:np;
    xx=(xp(:,ip)-xmin')./dx';
    i=fix(xx);
    xpn(:,ip)=xx-i;
    nep(:,ip)=i(1)+1+i(2)*ne(1);
end
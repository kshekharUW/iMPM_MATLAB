function [dx,Nn,Ne,x,icon,numOC,ntop,nbot,nright,nleft,freeDOFs]=initgrid(nx,xmin,xmax)
dx=[xmax-xmin]./(nx-1);
x1=linspace(xmin(1),xmax(1),nx(1));
y1=linspace(xmin(2),xmax(2),nx(2));
Nn=0;
for j=1:nx(2);
    for i=1:nx(1);
        Nn=Nn+1;
        x(Nn,:)=[x1(i),y1(j)];
    end
end
% element connectivity
ne=nx-1;
Ne=0;
for j=1:ne(2);
    for i=1:ne(1);
        Ne=Ne+1;
        icon(1,Ne)=(j-1)*nx(1)+i;
        icon(2,Ne)=(j-1)*nx(1)+i+1;
        icon(3,Ne)=j*nx(1)+i+1;
        icon(4,Ne)=j*nx(1)+i;
    end
end

% number of cells a node is part of
numOC=zeros(Nn,1);
for i=1:Nn
    numOC(i)=numel(find(icon==i));
end

% determine boundary nodes
freeDOFs=[1:Nn];
ntop=[];
nbot = [];
nleft=[];
nright = [];

for i = 1:length(x)
    if mod(i,nx(1))==0
        nright(length(nright)+1) = i;
        nleft(length(nleft)+1) = i-nx(1)+1;
    end
end
nbot = linspace(nleft(1),nright(1),nx(1));
ntop = linspace(nleft(end),nright(end),nx(1));

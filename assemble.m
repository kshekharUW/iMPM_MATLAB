function [Ksystem,C,DxM,DyM,G,M,DivM,detJ] = assemble( dx,Nn,Ne,icon )
    lx=dx(1); ly=dx(2);
    detJ=lx*ly/4;

    Ksystem = zeros(Nn,Nn);     
    C=zeros(Nn,Nn); %Npt*Np on all elements
    DxM=zeros(Nn,Nn); %Npt*Np,x on all elements
    DyM=zeros(Nn,Nn); %Npt*Np,y on all elements 
    M=zeros(2*Nn,2*Nn);% Nvt*Nv on all elements 
    G = zeros(2*Nn,Nn); %Nvt*Grad(Np)on all elements
    DivM = zeros(Nn,2*Nn); %Np*Div(Nv)on all elements
    
    kelement=1/(6*lx*ly)*((lx^2)*[2,1,-1,-2;
                                  1,2,-2,-1;
                                 -1,-2,2,1;
                                 -2,-1,1,2]+...
                          (ly^2)*[2,-2,-1,1;
                                 -2,2,1,-1;
                                 -1,1,2,-2;
                                  1,-1,-2,2]);
    
           
    celement=(lx*ly/36)*[4 2 1 2;
                         2 4 2 1;
                         1 2 4 2;
                         2 1 2 4];
    celementBig=zeros(8,8);
    celementBig(1:4,1:4)=celement; 
    celementBig(5:8,5:8)=celement;
    
    
                           
    ddx=ly/12*[-2 2 1 -1;
               -2 2 1 -1;
               -1 1 2 -2;
               -1 1 2 -2];

    ddy=lx/12*[-2 -1 1 2;
               -1 -2 2 1;
               -1 -2 2 1;
               -2 -1 1 2];
    gp=[ddx;ddy]; 
    divMsmall = [ddx,ddy];
            
    for i=1:Ne
        nodes=icon(:,i);
        C(nodes,nodes)=C(nodes,nodes)+celement;    
        DxM(nodes,nodes)=DxM(nodes,nodes)+ddx;  
        DyM(nodes,nodes)=DyM(nodes,nodes)+ddy;  
        Ksystem(nodes,nodes)=Ksystem(nodes,nodes)+kelement;
        G([nodes,nodes+Nn],nodes)=G([nodes,nodes+Nn],nodes)+gp;
        M([nodes,nodes+Nn],[nodes,nodes+Nn])=M([nodes,nodes+Nn],[nodes,nodes+Nn])+celementBig;
        DivM(nodes,[nodes,nodes+Nn])=DivM(nodes,[nodes,nodes+Nn])+divMsmall;
    end 
    
    Ksystem=sparse(Ksystem);
    DxM=sparse(DxM);
    DyM=sparse(DyM);
    C=sparse(C);
	G=sparse(G);
    M=sparse(M);
    
end %end function


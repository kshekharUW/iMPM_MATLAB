function [nodePressure]=solveP(L,DivM,vstar,qn,pSupportDOF,rhoWater,dt)
%     prhs = -rhoWater/dt * DivM*vstar + qn;
%     nodePressure(pFreeDOF) = L(pFreeDOF, pFreeDOF)\...
%         (prhs(pFreeDOF)-L(pFreeDOF,pSupportDOF)*pressureBC(pSupportDOF));
%     nodePressure(pSupportDOF) = pressureBC(pSupportDOF);
    
L(pSupportDOF,pSupportDOF)=1e20;
prhs = -rhoWater/dt * DivM*vstar + qn;  
prhs(pSupportDOF)=0;
nodePressure = L\prhs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linéarisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[A,B]=Linearise(X2e,Ue,Param)

A=  [0 -sqrt(Param.alpha*Param.beta);0 -sqrt(Param.alpha*Param.beta)];

B=  [-sqrt(Param.beta/Param.alpha);sqrt(Param.beta/Param.alpha)];
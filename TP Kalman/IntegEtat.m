
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Int�gration num�rique de l'�quation d'�tat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[Xkplus]=IntegEtat(Xk,u,Param,Te)

options=odeset;
[Auxt,AuxXkplus]=ode45(@EqEtatContinue,[0,Te],Xk,options,Param,u);

Xkplus=AuxXkplus(end,:)';
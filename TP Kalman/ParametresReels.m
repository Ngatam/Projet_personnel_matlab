
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D�finition des param�tres r�els
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Param�tres du mod�le
ParamReels.beta=1;
ParamReels.alpha=4;


%�tat d'�quilibre
ParamReels.Equilibre.X2e=sqrt(ParamReels.beta/ParamReels.alpha);
ParamReels.Equilibre.Ue=ParamReels.beta/ParamReels.Equilibre.X2e;

%Appel de la fonction Linearise pour calculer les matrices A et B
[ParamReels.A,ParamReels.B]=Linearise(ParamReels.Equilibre.X2e,...
    ParamReels.Equilibre.Ue,ParamReels);
ParamReels.C=[1,0];

%etat initial
ParamReels.X0=[4 ;ParamReels.Equilibre.X2e];

%etat initial du syst�me lin�aris�
ParamReels.X0lin=ParamReels.X0-[0;ParamReels.Equilibre.X2e];

%p�riode d'�chantillonnage
ParamReels.Te=0.1;

%Covariance des bruits de mesure
ParamReels.R=1e-2;
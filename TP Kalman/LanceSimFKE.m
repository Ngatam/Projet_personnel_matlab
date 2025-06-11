

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulations et Filtre de Kalman Etendu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all %ferme toutes les figures
clear all %efface le "workspace"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PARAMETRES DU SYSTEME REEL

eval('ParametresReels'); %execute le fichier ParametresReels.m

%temps de fin de simulation (en secondes)
ParamReels.tfin=10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PARAMETRES DU FILTRE

ParamFilt=ParamReels; %recopie des paramètres réels

%Estimé initial
ParamFilt.X00=ParamReels.X0;

%Covariance de l'erreur d'estimation initiale
ParamFilt.P00=diag([1e-3,1e-3]);

%Covariance des bruits de dynamique
ParamFilt.Q=diag([1e-5,1e-5]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SIMULATION


%execution du fichier Simulink : SimNL
%sim('SimNL',ParamReels.tfin);
load DataNLin

NbMes=length(YReel); %Nombre de points de mesure
%--------------------------------------------------------
%Bruitage des mesures

%Génération d'un vecteur, de meme taille que YReel,
%composé de de bruits gaussiens d'espérance nulle
% de covariance 1
BruitsCov1=randn(NbMes,1);

%Génération d'un vecteur, de meme taille que YReel,
%composé de de bruits gaussiens d'espérance nulle
% de covariance ParamReels.R
Bruits=sqrt(ParamReels.R)*BruitsCov1;

%addition des bruits à YReel -> Y mesuré
Ymes = YReel + Bruits;
%--------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FILTRAGE

%-----------------------------------------------------------
% initialisation

%mise en forme
Xkk=ParamFilt.X00;
Pkk=ParamFilt.P00;
U=UReel;

%Stockage
Resultat.HatX(1,:)=Xkk'; %Estimé
Resultat.DiagPkk(1,:)=diag(Pkk)'; %Composantes diag de Pkk
%-----------------------------------------------------------


%-----------------------------------------------------------
%Récurrence

for k=1:NbMes-1,

      %..........................................
      %Prédiction et intégration numérique
      Xkplus=IntegEtat(Xkk,U(k,1),ParamFilt,ParamFilt.Te);
      
      %linéarisation
      ParamFilt.A=[0 -U(k,1);0 -2*ParamFilt.alpha*Xkk(2,1)+U(k,1)];
      ParamFilt.B=Xkk;

       %discrétisation
      ParamFilt.F= expm(ParamFilt.A*ParamFilt.Te);
      Pkplus=ParamFilt.F*Pkk*ParamFilt.F'+ParamFilt.Q;
      K= Pkplus*ParamFilt.C'*inv(ParamFilt.C*Pkplus*ParamFilt.C'+ParamFilt.R);
      Xkk= Xkplus+K*(Ymes(k+1,1)-ParamFilt.C*Xkplus); 
      Pkk= (eye(2)-K*ParamFilt.C)*Pkplus; 

      %Stockage
      Resultat.HatX(k+1,:)=Xkk';
      Resultat.DiagPkk(k+1,:)=diag(Pkk)';
      Resultat.K(k+1,:)=K';

end;
%-----------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TRACES DES RESULTATS
numfig=0;
Texte={'X1','X2'};

% Tracé de l'état réel et de l'estimé sur une meme figure
    for i=1:2,
        numfig=numfig+1;
        figure(numfig)
        xlabel('t')
        hold on
        a=plot(Tps,XReel(:,i),'b-');
        b=plot(Tps,Resultat.HatX(:,i),'r-');
        legend([a,b],'Reel','Estime')
        title(['Estimation de ',Texte{i}])
        grid on
    end;

%erreur d'estimation
for i=1:2,
    numfig=numfig+1;
    figure(numfig)
    xlabel('t')
    hold on
    a=plot(Tps,XReel(:,i)-Resultat.HatX(:,i),'m-');
    title(['Erreur d estimation de ',Texte{i}])
    grid on
end;

%composantes diagonales de Pkk
%for i=1:2,
%    numfig=numfig+1;
%    figure(numfig)
%    xlabel('t')
%    hold on
%    a=plot(Tps,Resultat.DiagPkk(:,i),'g-');
%    title(['Composante diagonale de Pkk, numero ',num2str(i)])
%    grid on
%end;

%composantes de K
for i=1:2,
    numfig=numfig+1;
    figure(numfig)
    xlabel('t')
    hold on
    a=plot(Tps,Resultat.K(:,i),'k-');
    title(['Composante K, numero ',num2str(i)])
    grid on
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

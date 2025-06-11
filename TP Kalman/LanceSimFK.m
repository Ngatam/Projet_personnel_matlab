%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulations et Filtre de Kalman
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Mise à zézo du workspace

close all %ferme toutes les figures
clear all %efface le "workspace"
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PARAMETRES DU SYSTEME REEL

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
%--------------------------------------------------------
%Pour Discrétiser

%Définition de la représentation d'état linéaire à temps continu
SystCont=ss(ParamFilt.A,ParamFilt.B,ParamFilt.C,0);

%discrétisation
SystDiscr=c2d(SystCont,ParamFilt.Te);

%Mise en forme pour utilisation plus pratique
ParamFilt.F=SystDiscr.a;
ParamFilt.G=SystDiscr.b;
%--------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATION

%execution du fichier Simulink Linéaire : SimLin
%sim('SimLin',ParamReels.tfin);  si pas Simulink alors télécharger les données
load DataLin;  % Donées système linéarisé



%Execution du fichier Simulink Non linaire: SimNL
%sim('SimNL',ParamReels.tfin);   % a utiliser si vous avez simulink
%load DataNLin;  % Données système non-linéaire

%-----------------------------------------------------

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
%% FILTRAGE

%-----------------------------------------------------------
% initialisation

%mise en forme
Xkk=ParamFilt.X00-[0;ParamFilt.Equilibre.X2e];
Pkk=ParamFilt.P00;
U=UReel-ParamFilt.Equilibre.Ue*ones(NbMes,1);

%Stockage
Resultat.HatX(1,:)=Xkk'; %Estimé
Resultat.DiagPkk(1,:)=diag(Pkk)'; %Composantes diag de Pkk
%-----------------------------------------------------------


%-----------------------------------------------------------
%Récurrence

for k=1:NbMes-1,
   %...........................................................

   %Prediction
    Xkplus= ParamFilt.F*Xkk + ParamFilt.G*U(k,1);
    Pkplus= ParamFilt.F*Pkk*ParamFilt.F' + ParamFilt.Q;
   %Correction
    K= Pkplus*ParamFilt.C'*inv(ParamFilt.C*Pkplus*ParamFilt.C'+ParamFilt.R);
    Xkk= Xkplus+K*(Ymes(k+1,1)-ParamFilt.C*Xkplus);
    Pkk= (eye(2)-K*ParamFilt.C)*Pkplus;
   %...........................................................

   %...........................................................

   %Stockage
    % codé le Filtre de Kalman :
    Resultat.HatX(k+1,:)=Xkk';
    Resultat.DiagPkk(k+1,:)=diag(Pkk)';
    Resultat.K(k+1,:)=K';
   %...........................................................
end;
%-----------------------------------------------------------


%-----------------------------------------------------------
%Prise en compte de la linéarisation : x2 = X2lin + X2e
% et construction de l'erreur d'estimation
% De-commentez les lignes ci-dessous lorsque vous aurez
% codé le Filtre de Kalman :

Resultat.HatX(:,2)=Resultat.HatX(:,2)+...
ParamFilt.Equilibre.X2e*ones(NbMes,1);

ErrEstim = XReel-Resultat.HatX;
%-----------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TRACES
numfig=0;
Texte={'X1','X2'};

%estimé et valeur réelle
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
    %plot(Tps,XReel(:,i)-Resultat.HatX(:,i),'m-');
      plot(Tps,ErrEstim (:,i),'m-');
    title(['Erreur d estimation de ',Texte{i}])
    grid on
end;

%composantes diagonales de Pkk
%for i=1:2,
%    numfig=numfig+1;
%    figure(numfig)
%    xlabel('t')
%    hold on
%    plot(Tps,Resultat.DiagPkk(:,i),'g-');
%    title(['Composante diagonale de Pkk, numero ',num2str(i)])
%    grid on
%end;

%composantes de K
for i=1:2,
    numfig=numfig+1;
    figure(numfig)
    xlabel('t')
    hold on
    plot(Tps,Resultat.K(:,i),'k-');
    title(['Composante K, numero ',num2str(i)])
    grid on
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

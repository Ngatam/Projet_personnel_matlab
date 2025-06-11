
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulations et  de MonteCarlo pour Filtre de Kalman Etendu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all %ferme toutes les figures
clear all %efface le "workspace"
clc; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PARAMETRES DU SYSTEME REEL

eval('ParametresReels'); %execute le fichier ParametresReels.m

%temps de fin de simulation (en secondes)
ParamReels.tfin=10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------

%Etat initial Reel
ParamReelsX0Vrai=[0.1 ;ParamReels.Equilibre.X2e];
ParamReelsP00=diag([1e-10,1e-10]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PARAMETRES DU FILTRE

ParamFilt=ParamReels; %recopie des paramètres réels

%Estimé initial
ParamFilt.X00=ParamReels.X0;

%Covariance de l'erreur d'estimation initiale
ParamFilt.P00=diag([1e-5,1e-5]);

%Covariance des bruits de dynamique
ParamFilt.Q=diag([1e-10,1e-1]);

%Matrice de sortie
ParamFilt.C=[1,0];

%Nombre de simulations de MonteCarlo
ParamFiltNs=30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Execution du fichier Simulink Non linaire: SimNL
%sim('SimNL',ParamReels.tfin);   % a utiliser si vous avez simulink
load DataNLin;  % Données système non-linéaire

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SIMULATION ET FILTRAGE

for n=1:ParamFiltNs
    disp(['Simulation de Monte Carlo numero ',num2str(n)]);
    %-------------------------------------------------------------------------
    % SIMULATION

    %Etat initial
    ParamReels.X0=ParamReelsX0Vrai+sqrtm(ParamReelsP00)*randn(2,1);
    %simulation du système Non Linéaire 
    %donc nouvelle réalisation des bruits
   
     %sim('SimNL',ParamReels.tfin);  si vous avez Matlab 
     % Sinon passage par l'intégration du sysème d'tat 
      %............................
      %mise en forme
      Xkk=ParamFilt.X00;
      Pkk=ParamFilt.P00;
      U=UReel;
      XReel(1,:)=Xkk';    
     % boucle simulation avec intération numérique discrète 
    for k=1:length(Ymes)-1,
        u=U(k,1);
        %...........................................................
        %equation d'état
        Xkplus=IntegEtat(Xkk,u,ParamFilt,ParamFilt.Te);
        %equation d'observation
        YR=ParamReels.C*Xkk;
        Xkk=Xkplus;
        ...........................................................
        %Stockage
        XReel(k+1,:)=Xkk'; % Etat réel
        YReel(k,:)=YR; % Etat réel
      
        %...........................................................
     end;
     % RMZ   
    Xkk=[];
    
    %--------------------------------------------------------
    %Bruitage des mesures
    NbMes=length(YReel); %Nombre de points de mesure
    BruitsCov1=randn(NbMes,1);
    Bruits=sqrt(ParamReels.R)*BruitsCov1;
    Ymes = YReel + Bruits;  %sortie mesurée
    %-------------------------------------------------------------------------
    %mise en forme
    Xkk=ParamFilt.X00;
    Pkk=ParamFilt.P00;
    U=UReel;
    Y=Ymes; %sortie mesurée

    %Stockage
    Resultat.HatX(1,:)=Xkk'; %Estimé
    Resultat.DiagPkk(1,:)=diag(Pkk)'; %Composantes diag de Pkk

    %-------------------------------------------------------------------------
    % FILTRAGE

    for k=1:length(Y)-1,
        u=U(k,1);
        %...........................................................
        %PREDICTION

        Xkplus=IntegEtat(Xkk,u,ParamFilt,ParamFilt.Te);

        [A,B]=Linearise(Xkk(2,1),u,ParamFilt); %Linéarisation
        ParamFilt.F=expm(ParamFilt.Te*A); %discrétisation

        Pkplus=ParamFilt.F*Pkk*ParamFilt.F'+ParamFilt.Q;
        %...........................................................

        %...........................................................
        %Correction
        y=Y(k+1,1);
        K=Pkplus*ParamFilt.C'*inv(ParamFilt.C*Pkplus*ParamFilt.C'+ParamFilt.R);
        Xkk=Xkplus+K*(y-ParamFilt.C*Xkplus);
        Pkk=Pkplus-K*ParamFilt.C*Pkplus;
        %...........................................................

        %...........................................................
        %Stockage
        Resultat.HatX(k+1,:)=Xkk'; %Estimé
        Resultat.DiagPkk(k+1,:)=diag(Pkk)';
        Resultat.K(k+1,:)=K';
        %...........................................................
    end;
    %-------------------------------------------------------------------------

    %Calcul de l'erreur d'estimation
    Resultat.Err=XReel-Resultat.HatX;

    %-------------------------------------------------------------------------
    %Stockage pratique des résultats des simulations de MonteCarlo
    ResultatsMC.ErrX1(:,n)=Resultat.Err(:,1);
    ResultatsMC.ErrX2(:,n)=Resultat.Err(:,2);

    ResultatsMC.Diag11Pkk(:,n)=Resultat.DiagPkk(:,1);
    ResultatsMC.Diag22Pkk(:,n)=Resultat.DiagPkk(:,2);

    %-------------------------------------------------------------------------


end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Traitement des résultats des simulations de Monte-Carlo

%-------------------------------------------------------------------------
%Calcul de la moyenne empirique de l'erreur d'estimation

%Pour X1
ResultatsMC.MeanX1=mean(ResultatsMC.ErrX1')';
%Pour X2
ResultatsMC.MeanX2=mean(ResultatsMC.ErrX2')';

%Mise en forme pour faciliter les tracés
ResultatsMC.MeanErr=[ResultatsMC.MeanX1,ResultatsMC.MeanX2];

%-------------------------------------------------------------------------
%Calcul du max et du min de l'erreur d'estimation

%Pour X1
ResultatsMC.MaxErrX1=max(ResultatsMC.ErrX1')';
ResultatsMC.MinErrX1=min(ResultatsMC.ErrX1')';
%Pour X2
ResultatsMC.MaxErrX2=max(ResultatsMC.ErrX2')';
ResultatsMC.MinErrX2=min(ResultatsMC.ErrX2')';


%Mise en forme pour faciliter les tracés
ResultatsMC.MaxErr=[ResultatsMC.MaxErrX1,ResultatsMC.MaxErrX2];
ResultatsMC.MinErr=[ResultatsMC.MinErrX1,ResultatsMC.MinErrX2];
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
%Moyenne des composantes diagonales des covariances fournies par le filtre

%Pour X1
ResultatsMC.Mean11Pkk=mean(ResultatsMC.Diag11Pkk')';
%Pour X2
ResultatsMC.Mean22Pkk=mean(ResultatsMC.Diag22Pkk')';

%Mise en forme pour faciliter les tracés
ResultatsMC.MeanPkk=[ResultatsMC.Mean11Pkk,ResultatsMC.Mean22Pkk];
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%Covariance empirique de l'erreur d'estimation

%Pour X1
ResultatsMC.CovErrX1=diag(cov(ResultatsMC.ErrX1'));
%Pour X2
ResultatsMC.CovErrX2=diag(cov(ResultatsMC.ErrX2'));


%Mise en forme pour faciliter les tracés
ResultatsMC.CovErr=[ResultatsMC.CovErrX1,ResultatsMC.CovErrX2];
%-------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TRACES DES RESULTATS
Texte={'X1','X2'};
numfig=0;


% Moyenne empirique de l'erreur d'estimation, son min et son max
for i=1:2,
    numfig=numfig+1;
    figure(numfig)
    a=plot(Tps,ResultatsMC.MeanErr(:,i),'b-');
    hold on
    b=plot(Tps,ResultatsMC.MinErr(:,i),'g-');
    c=plot(Tps,ResultatsMC.MaxErr(:,i),'r-');
    legend([a,b,c],'Moyenne Empirique','Min','Max')
    title(['Erreur d estimation de ',Texte{i}]);
end;

% Covariances
for i=1:2,
    numfig=numfig+1;
    figure(numfig)
    a=plot(Tps,ResultatsMC.MeanPkk(:,i),'b-');
    hold on
    b=plot(Tps,ResultatsMC.CovErr(:,i),'m-');
    legend([a,b],'Moyenne des Pkk fournies par le filtre','Covariance empirique')
    title(['Covariance de l erreur d estimation de ',Texte{i}]);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

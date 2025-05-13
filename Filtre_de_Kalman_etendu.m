% Note: 
%  x_estimateur_suivant_actuel = x_estimateur_k+1_k
%  x_estimateur_suivant_suivant = x_estimateur_k+1_k+1


%-------------- Nettoyage  ------------------------------------------------
clear all
close all


%-------------- Import  ---------------------------------------------------
eval('ParamReels')


%-------------- Fonctions du système --------------------------------------
function [f_1, f_2] =  f(x_estimateur_k_k) % Définition de la fonction f
    f_1 = x_estimateur_k_k(1)^2 + x_estimateur_k_k(2)^2;
    f_2 = x_estimateur_k_k(1);
end


function [g_1, g_2] =  g(u_k_k) % Définition de la fonction g
    g_1 = 0;
    g_2 = 0;
end


function [phi_1, phi_2] =  phi(x_estimateur_k_k) % Définition de la fonction phi
    phi_1 = x_estimateur_k_k(1);
    phi_2 = x_estimateur_k_k(2)^3;
end


%-------------- Paramètres du sytème --------------------------------------
H_k = [1; 2];
Q_k = 10^(-4) ; % Choix de Q, modélise la confiance qu'on a sur le modèle
y = [[0.01; 0], [0.0001; 0]] ; % Mesures prises sur le systèmes  


%-------------- Initialisation --------------------------------------------
x_estimateur_0_0 = ParamReels.x_estimateur_0_0 ; % On va chercher x_estimateur_0/0 dans le fichier ParamReels
P_0_0 = ParamReels.P_0_0; % On va chercher P_0/0 dans le fichier ParamReels

x_estimateur_actuel_actuel = x_estimateur_0_0 ; % Initialisation de la valeur actuelle de l'estimateur
u_actuel_actuel = 0 ; % Initialisation de la valeur de la consigne
P_actuel_actuel = P_0_0; % Initialisation de la matrice de covariance


%-------------- Boucle ----------------------------------------------------
k = 0 ; % Initinalisation de l'incrément k
k_final = 2 ; % Choix de la dernière incrémentation calculé

while (k<k_final)

    %-------------- Prédiction --------------------------------------------
    x_estimateur_suivant_actuel = f(x_estimateur_actuel_actuel) + g(u_actuel_actuel) ; % Calcul de x_estimateur_k+1/k

    F_actuel = jacobian(f, x_estimateur_k_k) ; % Calcul de la jacobienne de f
    subs(F_actuel, x_estimateur_k_k, x_estimateur_actuel_actuel) % Calcul de la jacobienne au point actuel connaissant le point actuel

    P_suivant_actuel = F_actuel*P_actuel_actuel*F_actuel' + H_k*Q*H_k' ; % P_k/k+1 = F_k*P_k_k*F_k' + H_k*Q*H_k'


    %-------------- Mesure ------------------------------------------------
    y_suivant = y(:,k+1) ; % Extraction de la mesure k+1 de la matrice y


    %-------------- Correction --------------------------------------------
    C_suivant = jacobian(phi, x_estimateur_k_k) ; % Calcul de la jacobienne de phi
    subs(F_actuel, x_estimateur_k_k, x_estimateur_suivant_actuel) ; % Calcul de la jacobienne au point suivant connaissant le point actuel
    
    K_suivant = P_suivant_actuel*C_suivant' * inv((C_suivant*P_suivant_actuel*C_suivant' + R_actuel)) ; % Calcul du gain de Kalman au point k+1
    K_suivant = (C_suivant*P_suivant_actuel*C_suivant' + R_actuel)/P_suivant_actuel*C_suivant' ;  % Calcul du gain de Kalman au point k+1

    x_estimateur_suivant_suivant = x_estimateur_suivant_actuel + K_suivant(y_suivant - phi(x_estimateur_suivant_actuel)) ; % Calcul de l'estimateur au point k+1

    P_suivant_suivant = (eye(size(K_suivant)))


end;
%#---- D�finition des param�tres ----#
m = 1.75;
R = 24;
k_1 = 1.9*10^-4;
k_2 = -6.4*10^-4;
L_0 = 0.3714;
alpha = 1000;
g = 9.81;

%#---- D�finition des valeurs au point de fonctionnement ----#
Z_0 = 1.5e-2;
I_0 = (Z_0 + k_2) * sqrt(m*g/k_1);
U_0 = R*I_0;
V_0 = -alpha*Z_0;

%#---- Lin�arisation des fonctions f_1, f_2, g_1 et g_2 ----#
f_1 = (2*k_1*I_0^2)/(m*(Z_0 + k_2)^3); % f_1 <=> a_1
f_2 = -(2*k_1*I_0)/(m*(Z_0 + k_2)^2); % f_2 <=> a_2
g_1 = -(2*k_1*I_0)/(Z_0+k_2)^2; % g_2 <=> b_2
g_2 = (2*k_1 + L_0*(Z_0 + k_2)) / (Z_0 + k_2); % g_1 <=> b_1

%#---- D�finition des fonction A et B de la fonction de transfert P = B/A----#
B = [0, 0, 0, -f_2*alpha/(1000*g_2)]; %Cr�ation de B � partir des bi
A = [1,   R/(10*g_2), (f_2*g_1 - f_1*g_2)/(100*g_2),-f_1*R/(1000*g_2)]; %Cr�ation de A � partir des ai
suspension = tf(B,A); % suspension = fonction de transfert de B/A

%#---- Placement des p�les BO  - RST de d�gr� relatif 1 ----#
poles_systeme = roots(A);
% les p�les du syst�mes sont -3.6497 -5.7645 et -3.9174

%#---- Calcul des polyn�mes R et S ----#
% On choisit nos p�les 4 qui viennent du syst�mes (dont -3.6497 qu'on choisit 2 fois); et 3 p�les rapides en -10
A_BF = poly([poles_systeme(1) poles_systeme(1) poles_systeme(2) poles_systeme(3) -10 -10 -10]);
M = [A_BF(2)-A(1), A_BF(3)-A(2), A_BF(4)-A(3), A_BF(5), A_BF(6), A_BF(7), A_BF(8)]; % Creation de la matrice pour le syst�me avec la matrice de Sylvester

%#---- Matrice de Sylvester ----#
sylv = [
   A(1) 0    0      0    0    0     0;
   A(2) A(1) 0      0    0    0     0;
   A(3) A(2) A(1)   0    0    0     0;
   A(4) A(3) A(2)   B(4) 0    0     0;
   0    A(4) A(3)   B(3) B(4) 0     0;
   0    0    A(4)   B(2) B(3) B(4)  0;
   0    0    0      B(1) B(2) B(3)  B(4);
    ];

   
%#---- Resolution de l'�quation: sylv * [S R]' = M' et cr�ation des vecteur R et S----#
Vecteur_solution = inv(sylv)* (M'); 

sigma_1 = Vecteur_solution(1);
sigma_2 = Vecteur_solution(2);
sigma_3 = Vecteur_solution(3);

r_0 = Vecteur_solution(4);
r_1 = Vecteur_solution(5);
r_2 = Vecteur_solution(6);
r_3 = Vecteur_solution(7);

S = [1, sigma_1, sigma_2, sigma_3];
R = [r_0, r_1, r_2, r_3];

%Correction - donn�e lors du TP
%R = 1E+04 * [0.2930, 2.9537, 7.7461, 2.6330];
%S = 1E+03 * [0.001, 0.0409, 0.6818, 6.0776, 0];

V = poly([-1 , -1, -1, -1]); % Polyn�me pour rendre causaux les blocs R et T


%#---- Choix du T(s) ----#
T = 26.330 * [1 30 300 1000]; % alpha a �t� calcul� � la main

%#---- Analyse de la robutesse ----#
corRS = tf(R, S);
L = suspension*corRS;
margin(L);
Sens = 1/(1+L);
Mm = 1/norm(Sens, inf);

%#---- Simulation ----#
sim('TP_Automatique_simulink')
figure
subplot(211)
plot(t,U); % t et U viennent de la simulation
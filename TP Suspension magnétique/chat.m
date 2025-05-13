%#---- Définition des paramètres ----#
m = 1.75;
R = 24;
k_1 = 1.9*10^-4;
k_2 = -6.4*10^-4;
L_0 = 0.3714;
alpha = 1000;
g = 9.81;

%#---- Définition des valeurs au point de fonctionnement ----#
Z_0 = 1.5e-2;
I_0 = (Z_0 + k_2) * sqrt(m*g/k_1);
U_0 = R*I_0;
V_0 = -alpha*Z_0;

%#---- Linéarisation des fonctions f_1, f_2, g_1 et g_2 ----#
f_1 = (2*k_1*I_0^2)/(m*(Z_0 + k_2)^3); % f_1 <=> a_1
f_2 = -(2*k_1*I_0)/(m*(Z_0 + k_2)^2); % f_2 <=> a_2
g_1 = -(2*k_1*I_0)/(Z_0+k_2)^2; % g_2 <=> b_2
g_2 = (2*k_1 + L_0*(Z_0 + k_2)) / (Z_0 + k_2); % g_1 <=> b_1

%#---- Définition des fonction A et B de la fonction de transfert P = B/A----#
B = poly([-f_2*alpha/(1000*g_2), 0, 0, 0]);
A = poly([-f_1*R/(1000*g_2), (f_2*g_1 - f_1*g_2)/(100*g_2),  R/(10*g_2), 1]);

d = find(B ~= 0, 1, 'first') - 1;  % Calcul du retard pur
B = B(d+1:end);                   % On enlève les zéros initiaux

% --- Degrés des polynômes ---
na = length(A) - 1;
nb = length(B) - 1;
nr = nb;          % Degré de R
ns = na;          % Degré de S
nm = max(na + nr, nb + ns);  % Degré du polynôme de référence

% --- Placement des pôles souhaités ---
Am = poly([-3.6497 -3.6497 -5.7645 -3.9174 -10 -10 -10]);   % Exemple : pôles souhaités
Am = Am(1:nm+1);            % On ajuste le degré si besoin

% --- Construction de la matrice de Sylvester ---
% Forme : [A 0 ... ; 0 A ...] et [B 0 ... ; 0 B ...] pour générer la convolution
SA = zeros(nm+1, nr+1);
SB = zeros(nm+1, ns+1);

for i = 1:(nr+1)
    SA(i:i+na, i) = A(:);
end

for i = 1:(ns+1)
    SB(i:i+nb, i) = B(:);
end

S = [SA SB];  % Matrice de Sylvester

% --- Résolution ---
theta = S \ Am(:);  % Résolution du système S * theta = Am

R = theta(1:nr+1).';
S_corr = theta(nr+2:end).';

% --- Calcul de T pour le suivi (optionnel) ---
% On veut B(q^-1) * T(q^-1) = Am(q^-1)
T = conv(Am, 1);           % T sans réduction
T = T(end-nb:end);         % Ajuster le degré selon B

% --- Affichage ---
disp('Polynôme R :'), disp(R);
disp('Polynôme S :'), disp(S_corr);
disp('Polynôme T :'), disp(T);

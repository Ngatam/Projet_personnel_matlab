%#---- Définition des paramètres ----#
m = 1.75
R = 24*10^3
k_1 = 1.9*10^-4 
k_2 = -6.4*10^-4
L_0 = 0.3714
alpha = 1000
g = 9.81

%#---- Définition des valeurs au point de fonctionnement ----#
Z_0 = 1.5*10^-2
I_0 = (Z_0 + k_2) * (m*g/k_1)^1/2
U_0 = R*I_0

%#---- Linéarisation des fonctions f_1, f_2, g_1 et g_2 ----#
f_1 = (*k_1*I_0^2)/(m*(Z_0 + k_2)^3) % f_1 <=> a_1
f_2 = -(2*k_1*I_0)/(m*(Z_0 + k_2) % f_2 <=> a_2
g_1 = (2*k_1 + L_0*Z_0 + L_0*k_2 - 2*k_1*Z_0) / (Z_0 + k_2) % g_1 <=> b_1
g_2 = -(2*k_1*I_0)/(Z_0+k_2)^2 % g_2 <=> b_2

%#---- Définition des fonction A et B de la fonction de transfert ----#
B = poly([-f_2*alpha/(1000*g_2), 0, 0, 0]);
%A = poly([1, R/(10*g_2), (f_2*g_1 - f_1*g_2)/(100*g_2), -f_1*R/(1000*g_2)]);
A = poly([-f_1*R/(1000*g_2), (f_2*g_1 - f_1*g_2)/(100*g_2),  R/(10*g_2), 1]);
suspension = tf(B,A);% suspension = fonction de transfert de B/A









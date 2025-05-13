%Paramètres
Ke = 0.34;
Kt = 1.26;
Kl = 4.854;
tau = 8;
ksi = 0.082
wo = 1

ac = 1.1

Q = diag([0 0 0 0 0 1 1])
R = diag([1 1])


%Q = diag([2000 19 1 0.1 0.1])
%R = diag([3 0.1])

m = 2
p = 2
n = 5

%Matrices
A = [0 1 0 0 0; 
    -wo*wo -2*ksi*wo 0 0 0;
    0 0 0 1 0;
    0 0 0 0 0;
    0 0 Kl/tau 0 -1/tau];
B = [0 0;
    Ke*wo^2 Ke*wo^2;
    0 0;
    Kt -Kt;
    0 0]
C = [ 1 0 0 0 0;
    0 0 0 0 1]

Aa = [A zeros(n,2);
    C  zeros(2,2)]
Ba= [B;
    zeros(2,2)]
AaMal= Aa + al_c*eye(n+p)
K = lqr(AaMal,Ba,Q,R)
Kp = K(:,1:n)
Ki = K(:,n+1:n+p)


sim('rebi_multivariable')



subplot(2, 2 ,1)
plot(t, ug)
xlabel('Temps (s)')
ylabel('ug')
subplot(2, 2 ,2)
plot(t, ud)
xlabel('Temps (s)')
ylabel('ud')
subplot(2, 2 ,3)
plot(t, eps)
xlabel('Temps (s)')
ylabel('eps')
subplot(2, 2 ,4)
plot(t, v)
xlabel('Temps (s)')
ylabel('v')

s=tf('s')
L1= ss(A,B,Ki*C,0)
L1=L1/s
L2=ss(A,B,Kp,0)
L = L1+L2



L1bis=ss(A,(B*Ki)/s,-C,0)
L2

Lo=  1




















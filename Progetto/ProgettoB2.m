% Progetto regolatore per sistema di sollevamento
% 
% Controlli Automatici T
% Alessandro Orciari   
% Filippo Giorgi
% Marco Merlonghi
% 12/2023
%
% PARAMETRI FISICI SISTEMA
%
% Coefficiente di attrito       = 15 
% Coefficiente elasticità disco = 50 
%
% EQUILIBRIO
%
% posizione di equilibrio = 5pi/12 rad
%
% SPECIFICHE
% 
% -) Errore a regime in risposta a un gradino w(t) = 10(t) e d(t) = 10(t) pari a 0.01
%
% -) Attenuazione di almeno 60dB per d(t)
%       con [omega_d_min, omega_d_MAX] = [0,0.75]
%
% -) Attenuazione di almeno 75dB per n(t)
%       con [omega_n_min, omega_n_MAX] = [10^5,10^7]
%
% -) S% <= 20%
% -) Ta,1 <= 0.01 s
%
% M_f >= 30 gradi

clc; close all; clear all;

%solo per visualizzione, pulsazione minima e massima
omega_plot_min = 1e-3;
omega_plot_max = 1e7;

%% Parametri

k=50;
beta=15;

% Coefficienti momento d'inerzia
J_0=8;
J_1=0.4;
J_2=0.4;
J_3=0.5;
J_4=0.7;

psi_1=-0.1;
psi_2=2;
psi_3=2.5;
psi_4=-1;

% Momento d'inerzia
J = @(x) J_0 + J_1 * cos(x + psi_1)+ J_2 * cos(2*x + psi_2)+ J_3 * cos(3*x + psi_3)+ J_4 * cos(4*x + psi_4);
J_dot = @(x) -(J_1 * sin(x + psi_1)+ J_2 * 2 * sin(2*x + psi_2)+ J_3 * 3 * sin(3*x + psi_3)+ J_4 * 4 * sin(4*x + psi_4));


%% Coppia di equilibrio
%x1 = theta, x2 = omega, u = Cm
x1_e = 5/12*pi;
x2_e = 0;

x_e = [x1_e; x2_e];
u_e = k*x1_e;

%% Linearizzazione e calcolo FdT
A_e = [0, 1; ((-k*J(x1_e) -u_e*J_dot(x1_e) + beta*x2_e*J_dot(x1_e) + k*x1_e*J_dot(x1_e))/J(x1_e)^2), ((-beta*J(x1_e)-u_e*J_dot(x1_e) + beta*x2_e*J_dot(x1_e) + k*x1_e*J_dot(x1_e))/J(x1_e)^2)];
% con le opportune semplificazioni rimane A_e = [0,1;-k/J(x1_e),-beta/J(x1_e)];
B_e = [0;1/J(x1_e)];
C_e = [1,0];
D_e = 0;

modello = ss(A_e,B_e,C_e,D_e);
GG = tf(modello);
s=tf('s');

figure(1);
margin(GG);

%% Specifiche
% errore a regime
WW = 10;
DD = 10;
e_star = 0.01;

% attenuazione disturbo sull'uscita
A_d = 60;
omega_d_min = 1e-3; % fissato al lower bound del plot perchè lo 0 è irraggiungibile
omega_d_MAX = 0.75;

% attenuazione disturbo di misura
A_n = 75;
omega_n_min = 1e5;
omega_n_MAX = 1e7;

% Sovraelongazione massima e tempo d'assestamento all'1%
S_star = 20;
T_star = 0.01;

% Margine di fase
Mf_esp = 30;

%% Regolatore statico

static_gain_GG = abs(evalfr(GG,j*0));
mu_s_error = (WW+DD)/e_star/static_gain_GG; %nel modulo uno questo è riferito alla f ad anello quindi dobbiamo dividere
                                            %per il guadagno statico della G
mu_s_d = 10^(A_d/20)/abs(evalfr(GG, j*omega_d_MAX)); %deve essere di almeno A_d decibel per evitare al zona proibita

mu_s = max(mu_s_error, mu_s_d);
    
RR_s = mu_s;

GG_e = GG*RR_s;

figure(2);
hold on;

% Calcolo specifiche S% => Margine di fase
xi_star = abs(log(S_star/100))/sqrt(pi^2 + log(S_star/100)^2);
Mf      = max(xi_star*100,Mf_esp);

% Specifiche su d
Bnd_d_x = [omega_d_min; omega_d_MAX; omega_d_MAX; omega_d_min];
Bnd_d_y = [A_d; A_d; -150; -150];
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche su n
Bnd_n_x = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];
Bnd_n_y = [-A_n; -A_n; 100; 100];
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche tempo d'assestamento (minima pulsazione critica)
omega_Ta_min = 1e-3; % lower bound per il plot
omega_Ta_max = 460/(Mf*T_star); % omega_c >= 460/(Mf*T^*) ~ 4.6
Bnd_Ta_x = [omega_Ta_min; omega_Ta_max; omega_Ta_max; omega_Ta_min];
Bnd_Ta_y = [0; 0; -150; -150];
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_mag = ["A_d"; "A_n"; "\omega_{c,min}"; "G(j\omega)"];
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(GG_e,{omega_plot_min,omega_plot_max});
grid on; zoom on;


% Specifiche sovraelongazione (margine di fase)
omega_c_min = omega_Ta_max;
omega_c_max = omega_n_min;

phi_up = Mf - 180;
phi_low = -270; % lower bound per il plot

Bnd_Mf_x = [omega_c_min; omega_c_max; omega_c_max; omega_c_min];
Bnd_Mf_y = [phi_up; phi_up; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_arg = ["G_e(j\omega)"; "M_f"];
legend(Legend_arg);

%% Regolatore dinamico

Mf_star = Mf + 5;
omega_c_star = 1050;

arg_G_omega_c_star = angle(evalfr(GG_e, j*omega_c_star));
phi_star = Mf_star - 180 - rad2deg(arg_G_omega_c_star);

m_G_omega_c_star = abs(evalfr(GG_e, j*omega_c_star));
M_star = 1/m_G_omega_c_star; % viene da M_star = 10^((-|G_e(jwc)|dB/20)/10)

% Formule di inversione
tau = (M_star - cos(phi_star*pi/180))/(omega_c_star*sin(phi_star*pi/180));
alpha_tau = (cos(phi_star*pi/180) - 1/M_star)/(omega_c_star*sin(phi_star*pi/180));
alpha = alpha_tau / tau;

if min(tau,alpha) < 0 % Check correttezza parametri
    fprintf('Errore: parametri rete anticipatrice negativi');
    disp(tau);
    disp(alpha);
    return;
end

RR_d = (1 + tau*s)/(1 + alpha*tau*s); 

RR = RR_s*RR_d;

%% Funzione ad anello

LL = RR*GG; 
figure(3);
hold on;

% Specifiche su ampiezza
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
Legend_arg = ["L(j\omega)"; "M_f"];
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(LL,{omega_plot_min,omega_plot_max});
grid on; zoom on;

% Specifiche su fase
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
legend(Legend_arg);

% STOP qui per sistema con controllore dinamico + specifiche
if 0
    return;
end

%% Check prestazioni in anello chiuso

% Funzione di sensitività complementare
FF = LL/(1+LL);

% Risposta al gradino 1(t)
figure(4);

T_simulation = 2*T_star;
[y_step,t_step] = step(FF, T_simulation);
plot(t_step,y_step,'b');
grid on, zoom on, hold on;

LV = evalfr(FF,0); %calcolare il valore finale della risposta del sistema a un ingresso costante

% vincolo sovraelongazione
patch([0,T_simulation,T_simulation,0],[LV*(1+S_star/100),LV*(1+S_star/100),LV*2,LV*2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);

% vincolo tempo di assestamento all'1%
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1-0.01),LV*(1-0.01),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1+0.01),LV*(1+0.01),LV*2,LV*2],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

ylim([0,LV*2]);

Legend_step = ["Risposta al gradino"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);

%% Check disturbo in uscita

% Funzione di sensitività
SS = 1/(1+LL);
figure(5);

% Simulazione disturbo in uscita a pulsazione 0.15
omega_d = 0.15;
tt = 0:1e-2:2e2;
dd = 10*sin(omega_d*tt) + 10*sin(omega_d*2*tt) + 10*sin(omega_d*3*tt) + 10*sin(omega_d*4*tt);
y_d = lsim(SS,dd,tt); %linear simulation e salvo l'output dentro y_d
hold on, grid on, zoom on
plot(tt,dd,'m')
plot(tt,y_d,'b')
grid on
legend('d(t)','y_d(t)')

%% Check disturbo di misura

% Funzione di sensitività complementare
FF = LL/(1+LL);
figure(6);

% Simulazione disturbo di misura a pulsazione 100000
omega_n = 1e5;
tt = 0:1e-5:2*1e-3;
nn = sin(omega_n*tt) + sin(omega_n*2*tt) + sin(omega_n*3*tt) + sin(omega_n*4*tt);
y_n = lsim(FF,nn,tt);
hold on, grid on, zoom on
plot(tt,nn,'m')
plot(tt,y_n,'b')
grid on
legend('n(t)','y_n(t')

%% Punto opzionale 1
% Vedere il file Figura_progetto.m

%% Punto opzionale 2
figure(7); %mostriamo che il nostro gradino non rispetta il tempo di assestamento
hold on, grid on, zoom on
    simulation = x_e;
    out = sim("Simulink_punto_opzionale_2.slx");
    plot(out.simout);
    hold on;
xlabel('Tempo');
ylabel('Ampiezza');
yticks(-3:0.5:3);  
yticklabels({'-3', '-2.5', '-2', '-1.5', '-1', '-0.5', '0', '0.5', '1', '1.5', '2', '2.5', '3'});

patch([0,T_simulation,T_simulation,0],[LV*(1+S_star/100),LV*(1+S_star/100),LV*2,LV*2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);

% vincolo tempo di assestamento all'1%
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1-0.01),LV*(1-0.01),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1+0.01),LV*(1+0.01),LV*2,LV*2],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);

legend('x1 = 5\pi/12', 'Vincolo sovraelongazione', 'Vincolo tempo di assestamento inferiore', 'Vincolo tempo di assestamento superiore');




figure(8); %per quali valori di theta sono rispettati i vincoli
hold on, grid on, zoom on
for x1 = 23*pi/36:pi/36:26*pi/36
    simulation = [x1; 0];
    out = sim("Simulink_punto_opzionale_2.slx");
    plot(out.simout);
    hold on;
end
%26pi/36 130
%23pi/36 115
xlabel('Tempo');
ylabel('Ampiezza');
yticks(-3:0.5:3);  
yticklabels({'-3', '-2.5', '-2', '-1.5', '-1', '-0.5', '0', '0.5', '1', '1.5', '2', '2.5', '3'});

patch([0,T_simulation,T_simulation,0],[LV*(1+S_star/100),LV*(1+S_star/100),LV*2,LV*2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);

% vincolo tempo di assestamento all'1%
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1-0.01),LV*(1-0.01),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1+0.01),LV*(1+0.01),LV*2,LV*2],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);

legend('x1 = 23\pi/36', 'x1 = 24\pi/36', 'x1 = 25\pi/36', 'x1 = 26\pi/36', 'Vincolo sovraelongazione', 'Vincolo tempo di assestamento inferiore', 'Vincolo tempo di assestamento superiore');

%% Punto opzionale 3

figure(9); %per quali valori di w(t) vengono rispettati i vincoli
hold on, grid on, zoom on
for w = 1.050:0.001:1.054
    w_sim = w;
    out = sim("Simulink_punto_opzionale_3.slx");
    plot(out.simout);
    hold on;
end

xlabel('Tempo');
ylabel('Ampiezza');
yticks(-3:0.5:3);  
yticklabels({'-3', '-2.5', '-2', '-1.5', '-1', '-0.5', '0', '0.5', '1', '1.5', '2', '2.5', '3'});
patch([0,T_simulation,T_simulation,0],[LV*(1+S_star/100),LV*(1+S_star/100),LV*2,LV*2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);

% vincolo tempo di assestamento all'1%
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1-0.01),LV*(1-0.01),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1+0.01),LV*(1+0.01),LV*2,LV*2],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);

legend('w(t) = 1.050(t)', 'w(t) = 1.051(t)', 'w(t) = 1.052(t)', 'w(t) = 1.053(t)', 'w(t) = 1.054(t)', 'Vincolo sovraelongazione', 'Vincolo tempo di assestamento inferiore', 'Vincolo tempo di assestamento superiore');


% Ipotesi sul tempo di assestamento per sistema non linearizzato con tempo
% di assestamento = 0.01 e epsilon% = 5%

figure(10);
hold on, grid on, zoom on

out = sim("Simulink_progetto_non_linearizzato.slx");
plot(out.simout);
% vincolo sovraelongazione
patch([0,T_simulation,T_simulation,0],[LV*(1+S_star/100),LV*(1+S_star/100),LV*2,LV*2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);

% vincolo tempo di assestamento all'5%
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1-0.05),LV*(1-0.05),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1+0.05),LV*(1+0.05),LV*2,LV*2],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

Legend_step = ["Risposta al gradino sistema non linearizzato"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);


% Ipotesi sul tempo di assestamento per sistema non linearizzato con tempo
% di assestamento = 0.2 e epsilon% = 1%

figure(11);
hold on, grid on, zoom on

out = sim('Simulink_progetto_non_linearizzato.slx', 'StopTime', '0.3');
plot(out.simout);


% vincolo sovraelongazione
T_star = 0.2;
T_simulation = 0.3;

patch([0,T_simulation,T_simulation,0],[LV*(1+S_star/100),LV*(1+S_star/100),LV*2,LV*2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);

% vincolo tempo di assestamento all'1%
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1-0.01),LV*(1-0.01),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1+0.01),LV*(1+0.01),LV*2,LV*2],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);


Legend_step = ["Risposta al gradino sistema non linearizzato"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);



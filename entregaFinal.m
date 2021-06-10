% Equipo 4:                     Matrículas
% Daniel Degollado Rodríguez 	A00832555
% Yael Ariel Márquez Mas		A01745310
% Kailin Wu                     A00830574
% Alberto Horacio Orozco Ramos	A00831719
% Müsel Emmanuel Tabares Pardo 	A00830710

clear all; close all; clc;

%Parámetros
m = 0.01;           % Masa (kg)
g = 9.81;          % aceleración (m/s^2)
U = 1000000;       % Momento dipolar máximo
u0 = 4*pi*10^-7;   % Permeabilidad del vacío
a = 0.08;          % Radio (m)
R = 0.00009;       % Resistencia

% Función para calcular la ecuación
f1 = @(t,z,v) v;
f2 = @(t,z,v) -g - ((9.*(U*u0).^2.*a.^4)./(4.*R.*m)).* (z.^2./(z.^2+a.^2).^5).*v;

% Datos utilizados para el método de Runge-Kutta 4to orden
t0 = 0;
tf = 4;
z0 = 20;
v0 = 0;
h = 0.0001;

% Función de Runge-Kutta 4to orden
[t,z,v] = RK4(f1,f2,t0,tf,z0,v0,h);

% Vector para almacenar los valores de la aceleración
a = f2(t,z,v);

% Gráficas de posición, velocidad y aceleración
subplot(1,3,1);
plot(t,z);
title('Posición');
xlabel('Tiempo (s)');
ylabel('Altura (m)');
hold on;
grid on;
subplot(1,3,2);
plot(t,v,'-r');
title('Velocidad');
xlabel('Tiempo (s)');
ylabel('Velocidad (m/s)');
hold on;
grid on;
subplot(1,3,3);
plot(t,a,'-g');
title('Aceleración');
xlabel('Tiempo (s)');
ylabel('Aceleración (m/s^2)');
hold on;
grid on;

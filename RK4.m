% Equipo 4:                     Matrículas
% Daniel Degollado Rodríguez 	A00832555
% Yael Ariel Márquez Mas		A01745310
% Kailin Wu                     A00830574
% Alberto Horacio Orozco Ramos	A00831719
% Müsel Emmanuel Tabares Pardo 	A00830710

% Runge - Kutta 4to orden para obtener posición y velocidad respecto al
% tiempo

function[t,z,v] = RK4(f1,f2,t0,tf,z0,v0,h)

t = [t0:h:tf];
n = length(t);
z = zeros(1,n);
v = zeros(1,n);
z(1) = z0;
v(1) = v0;

for i = 1:n-1
    k1 = h * f1(t(i), z(i), v(i));
    l1 = h * f2(t(i), z(i), v(i));
    
    k2 = h * f1(t(i)+h/2, z(i)+k1/2, v(i)+l1/2);
    l2 = h * f2(t(i)+h/2, z(i)+k1/2, v(i)+l1/2);
    
    k3 = h * f1(t(i)+h/2, z(i)+k2/2, v(i)+l2/2);
    l3 = h * f2(t(i)+h/2, z(i)+k2/2, v(i)+l2/2);
    
    k4 = h * f1(t(i)+h/2, z(i)+k3/2, v(i)+l3/2);
    l4 = h * f2(t(i)+h/2, z(i)+k3/2, v(i)+l3/2);
    
    z(i+1) = z(i) + (k1 + 2*k2 + 2*k3 + k4)/6;
    v(i+1) = v(i) + (l1 + 2*l2 + 2*l3 + l4)/6;
    
end

end

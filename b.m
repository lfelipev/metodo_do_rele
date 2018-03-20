%% Metodo do rele

d0 = 0;
d1 = 0.5;
d2 = 2.5;
e0 = 1;
e = 0.9;
nptos = 100;
Tamostra = 0.03;
nptos1 = nptos/Tamostra;

% Encontrando o Kp caso desconhecido
kont = 0;
for t=2:nptos1
    if(data(t,1)~=data(t-1,1))
        kont = kont + 1;
        ch(kont) = t;
    end
end

Tu1 = (ch(6)-ch(5)) * Tamostra;
Tu2 = (ch(5)-ch(4)) * Tamostra;

aux1 = ch(3);
aux2 = ch(5);

i = 0;
for t = aux1:aux2
    i = i + 1;
    yi(i) = data(t,2);
    ui(i) = data(t,1);
    ti(i) = i*Tamostra;
end

a1 = 0.5*([0, yi]+[yi, 0]).*([ti 0]-[0 ti]);
a1 = sum(a1(1, 2:length(yi)));
a2 = 0.5*([0, ui]+[ui, 0]).*([ti 0]-[0 ti]);
a2 = sum(a2(1, 2:length(ui)));

kp = a1/a2;

% Encontrando Au e Ad
maximo = data(aux1,2);
for t=aux1:aux2
    if data(t,2) >= maximo
        maximo = data(t,2);
    end
end

minimo = data(aux1,2);
for t = aux1:aux2
    if data(t,2) <= minimo
        minimo = data(t,2);
    end
end

Au = maximo;
Ad = minimo;

% Funcao resultante B
s = tf('s');
B = 3* exp(-0.1*s)/(2*s+1);


%% Primeiro metodo
teta = log((d1-d0)*kp-e+e0) / ((d1-d0)*kp + Ad);
x1 = (d1+d2)*kp*exp(teta)-(d1-d0)*kp+e-e0;
x2 = (d0+d2)*kp-e0-e;
tau1 = Tu1/log(x1/x2);
theta1 = tau1 * teta;

% metodo 1
s = tf('s');
G1 = kp * exp(-theta1*s)/(tau1*s+1);

step(feedback(B, 1))
hold on
step(feedback(G1, 1))
hold on
legend('Sistema real', 'Modelo obtido')

%% Segundo metodo
h = abs(d2-d1);
a = abs(Au) + abs(Ad);
Tu = Tu1 + Tu2;
ku = 4 * h / (pi*sqrt(a^2 - e^2));

x = Tu/(2*pi);
y = sqrt((ku * kp)^2 - 1);
tau_1 = x*y;

y2 = pi-atan((2*pi/Tu) * tau_1);
theta_1 = x * y2;

% metodo 2
G2 = kp * exp(-theta_1*s)/(tau_1*s+1);


step(feedback(B, 1))
hold on
step(feedback(G2, 1))
hold on
legend('Sistema real', 'Modelo obtido')

%% Segunda ordem

tau_2 = (Tu / 2*pi) * (sqrt(ku * kp - 1));
theta_2 = (Tu / 2*pi) * (pi - 2*atan((2*pi/Tu) * tau_2));

% segunda ordem
G3 = kp * exp(-theta_2*s) / (tau_2*s + 1)^2;

step(feedback(G3, 1))

step(feedback(B, 1))
hold on
step(feedback(G1, 1))
hold on
legend('Sistema real', 'Modelo obtido')
% Variables no estocásticas
var C,C_M,C_H,h_M,h_H,k_M,k_H,z_M,z_H;

% Variables exógenas
varexo e_M,e_H;


% Parámetros del modelo
parameters beta,theta,delta,A,a,rho,eta,sigma_M,sigma_H,elasticidad,correlacion;
set_param_value('beta',beta);
set_param_value('theta',theta);
set_param_value('delta',delta);
set_param_value('rho',rho);
set_param_value('eta',eta);
set_param_value('sigma_M',sigma_M);
set_param_value('sigma_H',sigma_H);
set_param_value('elasticidad',elasticidad);
set_param_value('correlacion',correlacion);
set_param_value('A',A);
set_param_value('a',a);


% Modelo
model;
C=( a*((C_M)^(elasticidad)) + (1-a)*(C_H^(elasticidad))  )^(1/elasticidad);

C_M = exp(z_M) * ((k_M(-1))^(theta)) * (h_M^(1-theta) ) + (1-delta)*(k_M(-1)+k_H(-1)) -((k_M+k_H)) ;

C_H = exp(z_H) * (k_H(-1))^(eta) * (h_H^(1-eta));

z_M = rho*z_M(-1) + e_M;

z_H = rho* z_H(-1) + e_H;


(1/C)*a*((C_M)^(elasticidad-1)) = beta * ( (1/C(+1))*(C_M(+1))^(elasticidad-1) 
* (exp(z_M(+1)) *theta * (k_M)^(theta-1) * (h_M(+1)^(1-theta) )) +(1-delta)  );
(1/C)*a*(C_M)^(elasticidad-1) = beta*((1/C(+1))*a*(C_M(+1))^(elasticidad-1) * (1-delta) + (1/C(+1))*(1-a)* (C_H(+1))^(elasticidad-1) * exp(z_H(+1)) * eta * (k_H)^(eta-1) * (h_H(+1))^(1-eta)  );
(1/C) * a * (C_M)^(elasticidad-1) * (exp(z_M) * (k_M(-1))^(theta) * (1-theta) * (h_M)^(-1*theta) ) = (A/(1-h_M-h_H));
(1/C)*(1-a)*(C_H)^(elasticidad-1) * ( exp(z_H) * (k_H(-1))^(eta) * (1-eta) * (h_H)^(-1*eta) ) = (A/(1-h_M-h_H));


end;

% Valores iniciales
initval;

k_M=10;
k_H=2;
C_M= 0.5;
C_H= 0.5;
h_M=0.33;
h_H=0.28;
C=(a*C_M^(elasticidad)+(1-a)*C_H^(elasticidad))^(1/(elasticidad));
z_M=0;
z_H=0;
e_M=0;
e_H=0;
end;

% EE
steady;

% Condiciones Blanchard y Kahn
check;

%Choques
shocks;
var e_M; stderr sigma_M*sigma_M;
var e_H; stderr sigma_H*sigma_H;
var e_M, e_H = correlacion*sigma_M*sigma_H;
end;
%Simulacion
stoch_simul(order=1,periods=200,drop=20);
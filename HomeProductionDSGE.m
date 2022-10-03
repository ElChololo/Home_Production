
beta =0.99;
theta =0.36;
delta=0.025;
rho=0.95;
eta=0.08;
sigma_M=0.007;
sigma_H=0.007;
elasticidad=0.8;
correlacion=2/3;
% Calibración de A y a
h_m = 0.33;
h_h = 0.28;
k_m = (theta/ (1/beta) - 1 + delta )^(1/(1-theta)) * h_m;
k_h = (k_m/h_m) * h_h (eta * (1-theta)) / (theta*(1-eta));
c_h = k_h^(theta) * h_h**(e-eta);
c_m = k_m**(theta) * h_m**(1-theta) - delta*(k_m+k_h);
num = (1-eta) * (k_h/h_h)^(eta);
den = (1-theta) * ( (c_m/c_h)^(elasticidad-1) ) * (k_m/h_m)^(theta) + (1-eta) * ((k_h/h_h) ^(eta));
a = num/den;
C= a*((c_m)^(elasticidad)) + (1-a) * (c_h^(elasticidad));
A=(1-h_m-h_h)* (1/C) * (1-a) * (c_h^(elasticidad-1)) * (1-eta) * ((k_h)^(eta)) * (h_h^(-1*eta));

%ejecución

dynare HomeProduction.mod


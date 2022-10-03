h_m = 0.33
h_h = 0.28
eta =0.08
teta = 0.36
beta = 0.99
delta=0.025
epsilon = 0.8
k_m = (teta/( (1/beta) - 1 + delta) )**(1/(1-teta)) * h_m
k_h = (k_m/h_m) * h_h *  ( eta * (1-teta)/(teta*(1-eta)) ) 
c_h = k_h**(eta) * h_h**(1-eta)
c_m = k_m**(teta)*h_m**(1-teta) - delta*(k_m+k_h)
num = (1-eta)* (k_h/h_h)**eta
den = (1-teta) *( (c_m/c_h)**(epsilon-1)) * ((k_m/h_m) ** teta) + (1-eta) * (k_h/h_h) ** (eta)
a = num/den
print(k_m)
print(k_h)
print(c_h)
print(c_m)
print(a)
C= a*(c_m**epsilon) + (1-a)* (c_h**epsilon)
A = (1-h_m-h_h) * 1/C  * (1-a) *(c_h**(epsilon-1)) * (1-eta) * (k_h **(eta)) *(h_h **(-1*eta))
print(A)
# Pacotes 
import numpy as np 
import matplotlib.plotly as plt
from scipy.stats import norm

# Parametros
alpha = 0.25
sigma = 0.5 
beta = 0.96
kappa = 3.5/(n-1)
eps = np.finfo(float).eps
n = 10000
m = 10000
dz = 6*sigma/(n-1)
a = 1/(alpha*beta)


k = kappa +(np.transpose(np.array(range(1,n+1)))- 1) * kappa
z = -3*sigma +(np.transpose(np.array(range(1,m+1))) - 1)* dz

phi = norm.pdf(z/0.5)
q = phi/sum(norm.pdf(phi))

y = a*k**alpha

c_new = []
for j in range(0,n):
    c_colum =[]
    value = np.exp(z)*(y)-k[j]
    for i in value:
        if i > eps :
            c_colum.append(i)
        else:
            c_colum.append(eps)
    c_new.append(c_colum)

c = np.array(c_new)
u = np.log(c)
v = np.zeros((n,1))

for i in range(0,n):
    RHS = u + beta*(np.transpose(q)*v)
    Tv= RHS.max(0)
    j = np.array([1]*n)
    g = k[j]
    error = np.linalg.norm(Tv -v, np.inf)
    v=Tv

P = np.exp(z)*(alpha*beta*a*(k**alpha))

plt.plot(P)
plt.plot(v)
plt.show()
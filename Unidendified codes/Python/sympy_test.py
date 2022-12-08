import sympy as sp
from sympy.abc import a
x,k= sp.symbols('x k')
print(sp.inverse_fourier_transform(sp.sin(k-2)/(k-2),k,x))
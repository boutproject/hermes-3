from sympy import symbols, Function, lambdify
from sympy import sin, log, diff, cos, sqrt
from sympy.matrices import Matrix
import numpy as np
# see notes for formulae
# https://bout-dev.readthedocs.io/en/stable/user_docs/coordinates.html#the-perpendicular-laplacian-in-divergence-form

x = symbols('x')
y = symbols('y')
z = symbols('z')

g11 = 1.0# + 0.1*sin(x)*sin(y)*sin(z)
g22 = 1.0# + 0.1*sin(x)*sin(y)*sin(z)
g33 = 1.0# + 0.1*sin(x)*sin(y)*sin(z)
g12 = 0.0#*cos(x)*cos(y)*cos(z)
g23 = 0.5#*cos(x)*cos(y)*cos(z)
g13 = 0.0#*cos(x)*cos(y)*cos(z)
gup = Matrix(3,3,[g11,g12,g13,g12,g22,g23,g13,g23,g33])
gdown = gup.inv()
g_11 = gdown[0,0]
g_12 = gdown[0,1]
g_13 = gdown[0,2]
g_22 = gdown[1,1]
g_23 = gdown[1,2]
g_33 = gdown[2,2]

detgup = gup.det()
#detgup = g11*g22*g33 - g11*g23*g23 - g12*g12*g33 + g12*g13*g23 - g13*g13*g22 + g13*g12*g23
J = 1/sqrt(detgup)

g11_str = str(g11)
g22_str = str(g22)
g33_str = str(g33)
g12_str = str(g12)
g13_str = str(g13)
g23_str = str(g23)



f = x**2#*sin(2*y)*sin(2*z)
a = 1.0
print("f(x,y,z) = ",f)
print("a(x,y,z) = ",a)
fstr = str(f)
astr = str(a)

dfdx = diff(f,x)
dfdy = diff(f,y)
dfdz = diff(f,z)

df1 = dfdx - (g12/g22)*dfdy
df3 = dfdz - (g23/g22)*dfdy

a_grad_perp_f_x = a*J*(g11*df1 + g13*df3)
a_grad_perp_f_y = a*J*(g12*df1 + g23*df3)
a_grad_perp_f_z = a*J*(g13*df1 + g33*df3)

g = diff(f,x)

div_a_grad_perp_f = (1/J)*(diff(a_grad_perp_f_x,x)+diff(a_grad_perp_f_y,y)+diff(a_grad_perp_f_z,z)) 

print("div_a_grad_perp_f(x,y,z) = ",div_a_grad_perp_f)
div_a_grad_perp_f_str = str(div_a_grad_perp_f)
div_a_grad_perp_f_func = lambdify((x,y,z),div_a_grad_perp_f)

xval = 0.5
yval = 0.7
zval = 0.345
print(f"div_a_grad_perp_f({xval},{yval},{zval}) = ",div_a_grad_perp_f_func(xval,yval,zval))

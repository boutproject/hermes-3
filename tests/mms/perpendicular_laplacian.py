from sympy import symbols, Function, lambdify
from sympy import sin, log, diff, cos, sqrt, pi
from sympy.matrices import Matrix
import numpy as np

x = symbols('x')
y = symbols('y')
z = symbols('z')

# Div . ( a Grad_perp f )
# see notes for formulae
# https://bout-dev.readthedocs.io/en/stable/user_docs/coordinates.html#the-perpendicular-laplacian-in-divergence-form
def div_a_grad_perp_f_symbolic(g11, g12, g13, g22, g23, g33, a, f):
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
   #print("f(x,y,z) = ",f)
   #print("a(x,y,z) = ",a)

   dfdx = diff(f,x)
   dfdy = diff(f,y)
   dfdz = diff(f,z)

   df1 = dfdx - (g_12/g_22)*dfdy
   df3 = dfdz - (g_23/g_22)*dfdy

   a_grad_perp_f_x = a*J*(g11*df1 + g13*df3)
   a_grad_perp_f_y = a*J*(g12*df1 + g23*df3)
   a_grad_perp_f_z = a*J*(g13*df1 + g33*df3)

   div_a_grad_perp_f = (1/J)*(diff(a_grad_perp_f_x,x)+diff(a_grad_perp_f_y,y)+diff(a_grad_perp_f_z,z)) 

#   print("div_a_grad_perp_f(x,y,z) = ",div_a_grad_perp_f)
#   div_a_grad_perp_f_func = lambdify((x,y,z),div_a_grad_perp_f)
#   xval = 0.5
#   yval = 0.7
#   zval = 0.345
#   print(f"div_a_grad_perp_f({xval},{yval},{zval}) = ",div_a_grad_perp_f_func(xval,yval,zval))
   return div_a_grad_perp_f
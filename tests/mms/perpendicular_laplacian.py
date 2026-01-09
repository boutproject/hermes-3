from sympy import symbols, Function, lambdify
from sympy import sin, log, diff, cos, sqrt, pi
from sympy.matrices import Matrix
import numpy as np

from boutdata.mms import identity, x, y, z

# A class very similar to boutdata.mms.Metric
# to permit construction of similar functions
# and use for Div_par and Grad_par from boutdata.mms
# A key difference is that the symbol variables
# x, y, z are the same as the global Cartesian variables x, y, z
# rather than the transformed x', y', z'.
class GeneralMetric:
   def __init__(self,
         g11=1.0,
         g12=0.0,
         g13=0.0,
         g22=1.0,
         g23=0.0,
         g33=1.0):
      # obtain g_{ij} from g^{ij}
      gup = Matrix(3,3,[g11,g12,g13,g12,g22,g23,g13,g23,g33])
      gdown = gup.inv()
      g_11 = gdown[0,0]
      g_12 = gdown[0,1]
      g_13 = gdown[0,2]
      g_22 = gdown[1,1]
      g_23 = gdown[1,2]
      g_33 = gdown[2,2]
      detgup = gup.det()
      J = 1/sqrt(detgup)
      #detgup = g11*g22*g33 - g11*g23*g23 - g12*g12*g33 + g12*g13*g23 - g13*g13*g22 + g13*g12*g23

      self.g11 = g11
      self.g12 = g12
      self.g13 = g13
      self.g22 = g22
      self.g23 = g23
      self.g33 = g33
      self.g_11 = g_11
      self.g_12 = g_12
      self.g_13 = g_13
      self.g_22 = g_22
      self.g_23 = g_23
      self.g_33 = g_33
      self.J = J
      self.B = sqrt(g_22)/J
      # for the GeneralMetric, use the cartesian x, y, z
      # global variables from boutdata as the internal variables
      self.x = x
      self.y = y
      self.z = z

# Div . ( a Grad_perp f )
# see notes for formulae
# https://bout-dev.readthedocs.io/en/stable/user_docs/coordinates.html#the-perpendicular-laplacian-in-divergence-form
def Div_a_Grad_perp_f(a, f, metric=identity):
   g11, g12, g13, g22, g23, g33 = metric.g11, metric.g12, metric.g13, metric.g22, metric.g23, metric.g33
   g_11, g_12, g_13, g_22, g_23, g_33 = metric.g_11, metric.g_12, metric.g_13, metric.g_22, metric.g_23, metric.g_33
   J = metric.J
   # use the metric internal variables for derivatives
   x, y, z = metric.x, metric.y, metric.z

   dfdx = diff(f,x)
   dfdy = diff(f,y)
   dfdz = diff(f,z)

   df1 = dfdx - (g_12/g_22)*dfdy
   df3 = dfdz - (g_23/g_22)*dfdy

   a_grad_perp_f_x = a*J*(g11*df1 + g13*df3)
   a_grad_perp_f_y = a*J*(g12*df1 + g23*df3)
   a_grad_perp_f_z = a*J*(g13*df1 + g33*df3)

   div_a_grad_perp_f = (1/J)*(diff(a_grad_perp_f_x,x)+diff(a_grad_perp_f_y,y)+diff(a_grad_perp_f_z,z))
   return div_a_grad_perp_f

def Div_par_k_Grad_par_f(a, f, metric=identity):
   g11, g12, g13, g22, g23, g33 = metric.g11, metric.g12, metric.g13, metric.g22, metric.g23, metric.g33
   g_11, g_12, g_13, g_22, g_23, g_33 = metric.g_11, metric.g_12, metric.g_13, metric.g_22, metric.g_23, metric.g_33
   J = metric.J
   # use the metric internal variables for derivatives
   x, y, z = metric.x, metric.y, metric.z
   dfdy = diff(f,y)
   k_grad_par_f = a*J*dfdy/g_22
   div_par_k_grad_par_f = (1/J)*(diff(k_grad_par_f,y))
   return div_par_k_grad_par_f

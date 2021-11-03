import numpy as np

def blob_velocity(n, method='COM', return_index=False):
  """
  Calculate blob velocity in normalized time and normalized grid spacing
  
  Input: Blob density as a 3D vector in the form  n[t,x,z] where t is time and x,z are the perpendicular spatial coordinates
  
  Keywords
  --------
  
  method    Method to use:
            'peak' -> Calculate velocity of the peak density
            'COM' -> Calculate centre of mass velocity
  return_index   return indices used to create velocity
  
  """
  from boututils import calculus as Calc

  size = n.shape

  x = np.zeros(size[0])
  z = np.zeros(size[0])
  
  if method == 'peak':
    for i in np.arange(size[0]):
      nmax,nmin = np.amax((n[i,:,:])),np.amin((n[i,:,:]))
      xpos,zpos = np.where(n[i,:,:]==nmax)
      x[i] = xpos[0]
      z[i] = zpos[0]
      
  elif method == 'COM':
    x = np.zeros(size[0])
    z = np.zeros(size[0])
    for i in np.arange(size[0]):
      data = n[i,:,:] - n[0,0,0]   #use corner cell rather than nmin
      ntot = np.sum(data[:,:])
      
      z[i] = np.sum(np.sum(data[:,:],axis=0)*(np.arange(size[2]))) / ntot
      x[i] = np.sum(np.sum(data[:,:],axis=1)*(np.arange(size[1]))) / ntot

  else:
    raise ValueError("Invalid method {}".format(method))
      
  vx = Calc.deriv(x)
  vz = Calc.deriv(z)

  if return_index:
    return vx, vz, x, z
  else:
    return vx, vz

if __name__ == "__main__":
  from boutdata.collect import collect
  import pickle
  import sys

  path=sys.argv[1] # First argument is the path to the data

  n = collect('Ne', path=path, info=False)

  vx,vy,xx,yy = blob_velocity(n[:,:,0,:], method='COM', return_index=True)

  f = open('Velocity.dat','wb')
  pickle.dump(vx,f)
  f.close()
  
  f = open('Position.dat','wb')
  pickle.dump(xx,f)
  f.close()
  
  f = open('Velocity.dat','rb')
  vx = pickle.load(f)
  f.close()
  
  try:
    import matplotlib.pyplot as plt
    plt.plot(vx)
    plt.show()
  except ImportError:
    pass

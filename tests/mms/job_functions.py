# module to provide a test function for running the MMS tests of differential operators
import os
from xbout import open_boutdataset
import numpy as np
from scipy.optimize import curve_fit

def lin_func(x,b,a):
    return b*x + a
    
def collectvar(datasets, var, mesh=0):
    return datasets[mesh][var]

def run_manufactured_solutions_test(test_input):
   # expand inputs from input dictionary
   # the number of grids generated
   ntest = test_input["ntest"]
   # the minimum number of points in each of the x, y, z grids in the test
   # the number of points in the ith test is ngrid*i
   nnbase = test_input["ngrid"]
   # list of [name, symbolic string] pairs
   differential_operator_test_list = test_input["differential_operator_list"]
   astr = test_input["a_string"]
   fstr = test_input["f_string"]
   g11_str = test_input["g11_string"]
   g22_str = test_input["g22_string"]
   g33_str = test_input["g33_string"]
   g12_str = test_input["g12_string"]
   g13_str = test_input["g13_string"]
   g23_str = test_input["g23_string"]
   base_test_dir = test_input["test_dir"]
   interactive_plots = test_input["interactive_plots"]

   # create directory 
   if not os.path.isdir(base_test_dir):
       os.system("mkdir "+base_test_dir)
   workdirs = []
   # make test for each resolution based on template file
   for i in range(0,ntest):
      workdir = f"{base_test_dir}/slab-mms-test-{i}"
      workdirs.append(workdir)
      # create directory 
      if not os.path.isdir(workdir):
         os.system("mkdir "+workdir)
      # copy template
      file = workdir+"/BOUT.inp"
      os.system(f"cp BOUT.inp.template "+file)
      # update with mesh values for test
      nn = nnbase*(i+1) # number of points in each grid
      dd = 2.0*np.pi/nn # y z grid spacing, z, y on [0, 2pi]
      ddx = 1.0/(nn-4) # x grid spacing to account for guard cells, x on [0,1] 
      # symmetricGlobalX = true ensures that x = 0 and x = 1 sits
      # halfway between the last true grid point and the first guard point.
      with open(file,"a") as file:
         mesh_string = f"""
   [mesh]
   symmetricGlobalX = true
   extrapolate_y = false
   extrapolate_x= false
   extrapolate_z= false

   nx = {nn}
   dx = {ddx}
   ny = {nn}
   dy = {dd}
   nz = {nn}
   dz = {dd}

   g11 = {g11_str}
   g22 = {g22_str}
   g33 = {g33_str}
   g12 = {g12_str}
   g23 = {g23_str}
   g13 = {g13_str}
   x_input = x
   y_input = y
   z_input = z
   a = {astr}
   f = {fstr}
   """
         n_operators = len(differential_operator_test_list)
         mesh_string += f"""
   # information about differential operators
   n_operators = {n_operators}
   """
         for i in range(0,n_operators):
               mesh_string += f"""
   differential_operator_name_{i} = {differential_operator_test_list[i][0]}          
   expected_result_{i} = {differential_operator_test_list[i][1]}
   """
         file.write(mesh_string.replace("**","^"))

      # run job on this input
      print("../.././hermes_mms_tests -d "+workdir+" > "+workdir+"/output.txt")
      os.system("../.././hermes_mms_tests -d "+workdir+" > "+workdir+"/output.txt")

   # now analyse the results of the test
   # this slice avoids including guard cells in the test
   # need guard cells in x (assume 2 here) and guard cells in y (assume 1)
   # no guard cells in z
   s = slice(2, -2), slice(1, -1), slice(None)

   # a dictionary of plot data, filled later on
   plot_data = dict()

   # open the series of "workdir/BOUT.0.nc" files, 
   # saving them in a list `datasets`
   datasets = []
   for workdir in workdirs:
      boutmeshpath = workdir+"/"+f'BOUT.0.nc'
      boutinppath = workdir+"/"+'BOUT.inp'
      datasets.append(open_boutdataset(boutmeshpath, inputfilepath=boutinppath, keep_yboundaries=False))

   # make a easy scan over the two operators, generalisation to N operators possible
   for i, values in enumerate(differential_operator_test_list):
      label = values[0]
      l2norm = []
      nylist = []
      dylist = []
      for m in range(0,ntest):
         numerical = collectvar(datasets, f"result_{i}", m)
         expected = collectvar(datasets, f"expected_result_{i}", m)
         xx = collectvar(datasets, "x_input", m)
         yy = collectvar(datasets, "y_input", m)
         zz = collectvar(datasets, "z_input", m)
         ff = collectvar(datasets, "f", m)
         aa = collectvar(datasets, "a", m)
      
         error_values = (numerical - expected)[s]
         thisl2 = np.sqrt(np.mean(error_values**2))
         l2norm.append(thisl2)
         nylist.append(yy.shape[1])
         # proxy for grid spacing
         dylist.append(1.0/yy.shape[1])
               
      # cast lists as numpy arrays for further manipulation
      l2norm = np.array(l2norm)
      #print("test error: ",l2norm)
      nylist = np.array(nylist)
      dylist = np.array(dylist)
      
      # find linear fit coefficients to test convergence rate
      # and construct fit function for plotting
      try:
         logl2norm = np.log(l2norm)
         logdylist = np.log(dylist)
         outvals = curve_fit(lin_func,logdylist,logl2norm)
         coeffs = outvals[0]
         slope = coeffs[0] # the convergence order
         offset = coeffs[1]
         logfit = slope*logdylist + offset
         fitfunc = np.exp(logfit)
      except ValueError:
         print("Infs/Nans encountered in fit, skipping")
         slope = None
         offset = None
         fitfunc = None
      
      # record results in dictionary and plot
      #label = attrs["operator"] + " : f = " + attrs["inp"]
      #label = "FV::Div_a_Grad_perp(a, f)"
      plot_data[label] = [dylist, l2norm, fitfunc, slope, offset]

   # plot the results
   try:
      import matplotlib.pyplot as plt
      ifig = 0
      for key, variable_set in plot_data.items():
         (xaxis, yaxis, fit, slope, offset) = variable_set
         plt.figure()
         plt.plot(xaxis, yaxis, "x-", label="$\\epsilon(\\mathcal{L}\\ast f)$: "+key)
         plt.plot(xaxis, yaxis[0]*(xaxis/xaxis[0])**2, "x-", label="$\\propto \\Delta^2$")
         if not fit is None:
               plt.plot(xaxis, fit, "x-", label="$e^{{{:.2f}}}\\Delta^{{{:.2f}}}$".format(offset,slope))
         plt.xlabel("$\\Delta = 1/N_y$")
         plt.title(key)
         plt.legend()
         if not fit is None:
               plt.gca().set_yscale("log")
               plt.gca().set_xscale("log")
         else:
               print("l2 error: ",yaxis)
         plt.savefig(f"{base_test_dir}/fig_{ifig}.png")
         if interactive_plots:
             plt.show()
         plt.close()
         ifig+=1
   except:
   # Plotting could fail for any number of reasons, and the actual
   # error raised may depend on, among other things, the current
   # matplotlib backend, so catch everything
      pass

   # test the convergence rates
   success = True
   output_message = ""
   for key, variable_set in plot_data.items():
      this_test_success = True
      (xaxis, yaxis, fit, slope, offset) = variable_set
      # check slope of fit ~= 2
      if not slope is None:
         if slope < 1.95:
               this_test_success = False
      else: # or permit near-zero errors, but nothing larger
         for error in yaxis:
               if error > 1.0e-10:
                  this_test_success = False
      # append test message and set global success variable
      if this_test_success:
         output_message += f"{key} convergence order {slope:.2f} => Test passed \n"
      else:
         output_message += f"{key} convergence order {slope:.2f}  => Test failed \n"
         success = False

   return success, output_message

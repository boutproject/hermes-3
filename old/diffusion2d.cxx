
#include "diffusion2d.hxx"

#include <bout/constants.hxx>
#include <bout/fv_ops.hxx>
#include "div_ops.hxx"

using bout::globals::mesh;

Diffusion2D::Diffusion2D(Solver *solver, Mesh*, Options &options) : NeutralModel(options) {
  // 2D (X-Z) diffusive model
  // Neutral gas dynamics
  solver->add(Nn, "Nn");
  solver->add(Pn, "Pn");
  
  Dnn = 0.0; // Neutral gas diffusion
  
  SAVE_REPEAT(Dnn);

  Lmax = options["Lmax"].doc("Maximum mean free path [m]").withDefault(1.0);

  // Set Laplacian inversion to null
  inv = 0;
}

Diffusion2D::~Diffusion2D() {
  if(inv) {
    delete inv;
  }
}

void Diffusion2D::update(const Field3D &Ne, const Field3D &Te, const Field3D &UNUSED(Ti), const Field3D &UNUSED(Vi)) {
  
  mesh->communicate(Nn, Pn);
  
  Nn = floor(Nn, 1e-8);
  Field3D Tn = Pn / Nn; 
  Tn = floor(Tn, 0.01/Tnorm);
  
  Field3D Nelim = floor(Ne, 1e-19); // Smaller limit for rate coefficients
  
  // Calculate atomic processes
  for(int i=0;i<mesh->LocalNx;i++)
    for(int j=0;j<mesh->LocalNy;j++)
      for(int k=0;k<mesh->LocalNz;k++) {
        // Charge exchange frequency, normalised to ion cyclotron frequency
        BoutReal sigma_cx = Nelim(i,j,k)*Nnorm*hydrogen.chargeExchange(Te(i,j,k)*Tnorm)/Fnorm;
	
        // Ionisation frequency, normalised to ion cyclotron frequency
        BoutReal sigma_iz = Nelim(i,j,k)*Nnorm*Nn(i,j,k) * hydrogen.ionisation(Te(i,j,k)*Tnorm)/Fnorm;
        
        // Neutral thermal velocity
        BoutReal vth_n = sqrt(Tn(i,j,k)); // Normalised to Cs0
	  
        // Neutral-neutral mean free path
        BoutReal a0 = PI*SQ(5.29e-11);
        BoutReal lambda_nn = 1. / (Nnorm*Nn(i,j,k)*a0); // meters
        if(lambda_nn > Lmax) {
          // Limit maximum mean free path
          lambda_nn = Lmax;
        }
	
        lambda_nn /= Lnorm; // Normalised length to Lnorm
        // Neutral-Neutral collision rate, normalised to ion cyclotron frequency
        BoutReal sigma_nn = vth_n / lambda_nn;
	
        // Total neutral collision frequency, normalised to ion cyclotron frequency
        BoutReal sigma = sigma_cx + sigma_nn + sigma_iz;
	
        // Neutral gas diffusion
        Dnn(i,j,k) = SQ(vth_n) / sigma;
	
        // Rates
        BoutReal R_rc = SQ(Nelim(i,j,k)) * hydrogen.recombination(Nelim(i,j,k)*Nnorm, Te(i,j,k)*Tnorm) * Nnorm / Fnorm; // Time normalisation
        BoutReal R_iz = Nelim(i,j,k)*Nn(i,j,k) * hydrogen.ionisation(Te(i,j,k)*Tnorm) * Nnorm / Fnorm; // Time normalisation
        BoutReal R_cx = sigma_cx * Nn(i,j,k);
        
        // Plasma sink / neutral source
        S(i,j,k) = R_rc - R_iz;
        
        // Power transfer from plasma to neutrals
	
        Qi(i,j,k) = R_cx * (3./2)*(Te(i,j,k) - Tn(i,j,k));
        
        // Power transfer due to ionisation and recombination
        Qi(i,j,k) += (3./2)*(Te(i,j,k)*R_rc - Tn(i,j,k)*R_iz);
        
        // Ion-neutral friction
        Fperp(i,j,k) = R_cx  // Charge-Exchange
          + R_rc; // Recombination
        
        // Radiated power from plasma
        // Factor of 1.09 so that recombination becomes an energy source at 5.25eV
        Rp(i,j,k) = (1.09*Te(i,j,k) - 13.6/Tnorm)*R_rc
          + (Eionize/Tnorm) * R_iz; // Ionisation energy
        
      }
  
  // Neutral density
  ddt(Nn) = 
    + S 
    + FV::Div_a_Laplace_perp(Dnn, Nn);
  
  // Neutral pressure
  ddt(Pn) = (2./3)*Qi
    + FV::Div_a_Laplace_perp(Dnn, Pn)
    ;
  
}

void Diffusion2D::precon(BoutReal, BoutReal gamma, BoutReal) {
  // Neutral gas diffusion
  // Solve (1 - gamma*Dnn*Delp2)^{-1} 
  if(!inv) {
    inv = Laplacian::create();
    // Zero value outer boundary
    
    inv->setInnerBoundaryFlags(INVERT_DC_GRAD | INVERT_AC_GRAD);
    
    inv->setCoefA(1.0);
  }
  
  inv->setCoefD(-gamma*Dnn);
  
  ddt(Nn) = inv->solve(ddt(Nn));
  
  ddt(Pn) = inv->solve(ddt(Pn));
}

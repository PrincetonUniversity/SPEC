#include <biest.hpp>
#include <sctl.hpp>

int main(int argc, char** argv) {
#ifdef SCTL_HAVE_PETSC
  PetscInitialize(&argc, &argv, NULL, NULL);
#elif defined(SCTL_HAVE_MPI)
  MPI_Init(&argc, &argv);
#endif

  { // Compute vacuum field and Taylor state
    typedef double Real;
    sctl::Comm comm = sctl::Comm::Self();
    sctl::Profile::Enable(true);

    sctl::Vector<biest::Surface<Real>> Svec(1);
    { // Initialize Svec
      long Nt = 70*28, Np = 14*28;
      Svec[0] = biest::Surface<Real>(Nt, Np, biest::SurfType::W7X_);
      // To modify the surface,
      // for (long t = 0; t < Nt; t++) {
      //   for (long p = 0; p < Np; p++) {
      //     Svec[0].Coord()[0*Nt*Np + t*Np + p] = X(theta(t),phi(p));
      //     Svec[0].Coord()[1*Nt*Np + t*Np + p] = Y(theta(t),phi(p));
      //     Svec[0].Coord()[2*Nt*Np + t*Np + p] = Z(theta(t),phi(p));
      //   }
      // }
    }

    { // Compute vacuum field
      sctl::Vector<Real> B, J;
      Real tor_flux = 1.0;
      //biest::VacuumField<Real,1,55,75>::Compute(B, J, tor_flux, 0, Svec, comm, 1e-12, 100);

      // Generate VTK visualization
      //WriteVTK("B0", Svec, B, comm);
      //B.Write("tmp-B0.data");
      //Svec[0].Coord().Write("svec0.data");
    }
    
    { // Compute Taylor state
      sctl::Vector<Real> B;
      Real lambda = 1.0, tor_flux = -1.0;
      biest::TaylorState<Real,1,40,65>::Compute(B, tor_flux, 0, lambda, Svec, comm, 1e-11, 100);

      // Generate VTK visualization
      WriteVTK("B1", Svec, B, comm);
      B.Write("tmp-B1.data");
      Svec[0].Coord().Write("svec1.data");
      }

    // Print profiling output
    sctl::Profile::print(&comm);
  }

#ifdef SCTL_HAVE_PETSC
  PetscFinalize();
#elif defined(SCTL_HAVE_MPI)
  MPI_Finalize();
#endif
  return 0;
}


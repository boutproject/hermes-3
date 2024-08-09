#include "./boundary_iterator.hxx"
#include "bout/parallel_boundary_region.hxx"

class YBoundary {
public:
  template <class T>
  void iter_regions(const T& f) {
    ASSERT1(is_init);
    for (auto& region : boundary_regions) {
      f(*region);
    }
    for (auto& region : boundary_regions_par) {
      f(*region);
    }
  }

  template <class F>
  void iter(const F& f){
    return iter_regions(f);
  }

  void init(Options& options, Mesh* mesh=nullptr){
    if (mesh == nullptr) {
      mesh = bout::globals::mesh;
    }

    bool lower_y = options["lower_y"].doc("Boundary on lower y?").withDefault<bool>(true);
    bool upper_y = options["upper_y"].doc("Boundary on upper y?").withDefault<bool>(true);
    bool outer_x = options["outer_x"].doc("Boundary on inner y?").withDefault<bool>(true);
    bool inner_x = options["inner_x"].doc("Boundary on outer y?").withDefault<bool>(false);

    if (mesh->isFci()) {
      if (outer_x) {
	for (auto& bndry : mesh->getBoundariesPar(BoundaryParType::xout)) {
	  boundary_regions_par.push_back(bndry);
	}
      }
      if (inner_x) {
	for (auto& bndry : mesh->getBoundariesPar(BoundaryParType::xin)) {
	  boundary_regions_par.push_back(bndry);
	}
      }
    } else {
      if (lower_y) {
	boundary_regions.push_back(std::make_shared<NewBoundaryRegionY>(mesh, true, mesh->iterateBndryLowerY()));
      }
      if (upper_y) {
	boundary_regions.push_back(std::make_shared<NewBoundaryRegionY>(mesh, false, mesh->iterateBndryUpperY()));
      }
    }
    is_init=true;
  }


private:
  std::vector<std::shared_ptr<BoundaryRegionPar>> boundary_regions_par;
  std::vector<std::shared_ptr<NewBoundaryRegionY>> boundary_regions;

  bool is_init{false};
};

/// Limited free gradient of log of a quantity
/// This ensures that the guard cell values remain positive
/// while also ensuring that the quantity never increases
///
///  fm  fc | fp
///         ^ boundary
///
/// exp( 2*log(fc) - log(fm) )
///
inline BoutReal limitFree(BoutReal fm, BoutReal fc) {
  if (fm < fc) {
    return fc; // Neumann rather than increasing into boundary
  }
  if (fm < 1e-10) {
    return fc; // Low / no density condition
  }
  BoutReal fp = SQ(fc) / fm;
#if CHECKLEVEL >= 2
  if (!std::isfinite(fp)) {
    throw BoutException("SheathBoundaryParallel limitFree: {}, {} -> {}", fm, fc, fp);
  }
#endif

  return fp;
}

inline BoutReal limitFree(const Field3D& f, const BoundaryRegionParIter& pnt) {
  if (pnt.valid() > 0) {
    return limitFree(pnt.yprev(f), f[pnt.ind()]);
  }
  return f[pnt.ind()];
}

inline BoutReal limitFree(const Field3D& f, const BoundaryRegionIter& pnt) {
  return limitFree(pnt.yprev(f), f[pnt.ind()]);
}


extern YBoundary yboundary;

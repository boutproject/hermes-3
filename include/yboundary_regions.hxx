#pragma once

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
    bool outer_x = options["outer_x"].doc("Boundary on outer x?").withDefault<bool>(true);
    bool inner_x = options["inner_x"].doc("Boundary on inner x?").withDefault<bool>(false);

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

extern YBoundary yboundary;

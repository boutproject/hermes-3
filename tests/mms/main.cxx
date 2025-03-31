#include "bout/field_factory.hxx"

//#include "hermes-2.hxx"
#include "bout/version.hxx"
#include "div_ops.hxx"
#include "bout/fv_ops.hxx"
#include "bout/difops.hxx"
#include "../include/div_ops.hxx"

// Class used to store function of one argument
// in operator list below
class nameandfunction1 {
public:
  std::string name;
  std::function<Field3D(const Field3D&)> func;
};

// Class used to store function of two arguments
// in operator list below
class nameandfunction2 {
public:
  std::string name;
  std::function<Field3D(const Field3D&, const Field3D&)> func;
};

class nameandfunction3 {
public:
  std::string name;
  std::function<Field3D(const Field3D&, const Field3D&, Field3D&)> func;
};

// Class used to store function of four arguments
// in operator list below. Second pair of arguments are
// internal diagnostics to the function
class nameandfunction4 {
  public:
    std::string name;
    std::function<Field3D(const Field3D&, const Field3D&, Field3D&, Field3D&)> func;
  };
  

// List of tested operators of 1 arguments
const auto differential_operators_1_arg = {
    nameandfunction1{"Div_par(f)",
                 [](const Field3D &f) {
                   return Div_par(f);
                 }},
    nameandfunction1{"Grad_par(f)",
                 [](const Field3D &f) {
                   return Grad_par(f);
                 }},
};
// List of tested operators of 2 arguments
const auto differential_operators_2_arg = {
    nameandfunction2{"FV::Div_a_Grad_perp(a, f)",
                 [](const Field3D &a, const Field3D &f) {
                   return FV::Div_a_Grad_perp(a, f);
                 }},
};
// List of tested operators of 3 arguments
const auto differential_operators_3_arg = {
  nameandfunction3{"Div_par_K_Grad_par_mod(a, f)",
               [](const Field3D &a, const Field3D &f, Field3D &flow_ylow) {
                 return Div_par_K_Grad_par_mod(a, f, flow_ylow, true);
               }},
};
// List of tested operators of 4 arguments
const auto differential_operators_4_arg = {
    nameandfunction4{"Div_a_Grad_perp_nonorthog(a, f)",
                 [](const Field3D &a, const Field3D &f, 
                    Field3D &flow_xlow, Field3D &flow_ylow) {
                   return Div_a_Grad_perp_nonorthog(a, f, flow_xlow, flow_ylow);
                 }},
    nameandfunction4{"Div_a_Grad_perp_flows(a, f)",
                  [](const Field3D &a, const Field3D &f, 
                     Field3D &flow_xlow, Field3D &flow_ylow) {
                    return Div_a_Grad_perp_flows(a, f, flow_xlow, flow_ylow);
                  }},
};

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);
  
  int n_operators = Options::root()["mesh"]["n_operators"].withDefault(1);
  
  Mesh* mesh = Mesh::create(&Options::root()["mesh"]);
  mesh->load();

  Options dump;

  // Note that sticking all the inputs in the [mesh] section is a bit of a hack, but makes
  // for a more compact input file than using the 'standard' variable-initialisation
  // routines, which would require a separate section for each variable.

  // Load and save coordinate variables
  Field3D xx{mesh}, yy{mesh}, zz{mesh};
  mesh->get(xx, "x_input", 0.0, false);
  mesh->get(yy, "y_input", 0.0, false);
  mesh->get(zz, "z_input", 0.0, false);
  dump["x_input"] = xx;
  dump["y_input"] = yy;
  dump["z_input"] = zz;

  // Get coefficient from input file
  Field3D a{mesh};
  mesh->get(a, "a", 1.0, false);
  mesh->communicate(a);
  dump["a"] = a;

  // Get test variable from input file
  Field3D f{mesh};
  mesh->get(f, "f", 0.0, false);
  mesh->communicate(f);
  dump["f"] = f;

  // diagnostic variables
  Field3D flow_xlow{mesh};
  Field3D flow_ylow{mesh};

  for (int i = 0; i < n_operators; i++){
      std::string inputname = "differential_operator_name_"+std::to_string(i);
      std::string expectedname = "expected_result_"+std::to_string(i);
      std::string outname = "result_"+std::to_string(i);
      std::string outname_flow_xlow = "result_flow_xlow_"+std::to_string(i);
      std::string outname_flow_ylow = "result_flow_ylow_"+std::to_string(i);
      std::string differential_operator_name = Options::root()["mesh"][inputname].withDefault("FV::Div_a_Grad_perp(a, f)");
      // the for loop and if statement below should be replaced
      // by a neater indexing syntax below if possible
      for (const auto& difop: differential_operators_1_arg) {
          if (difop.name.compare(differential_operator_name) == 0){
              // Get result of applying the named differential operator
              Field3D result = difop.func(f);
              dump[outname] = result;
              dump[outname].setAttributes({
                    {"operator", difop.name},
                });
              // Get expected result from input file
              Field3D expected_result{mesh};
              mesh->get(expected_result, expectedname, 0.0, false);
              dump[expectedname] = expected_result;
          }
      }
      for (const auto& difop: differential_operators_2_arg) {
          if (difop.name.compare(differential_operator_name) == 0){
              // Get result of applying the named differential operator
              Field3D result = difop.func(a, f);
              dump[outname] = result;
              dump[outname].setAttributes({
                    {"operator", difop.name},
                });
              // Get expected result from input file
              Field3D expected_result{mesh};
              mesh->get(expected_result, expectedname, 0.0, false);
              dump[expectedname] = expected_result;
          }
      }
      for (const auto& difop: differential_operators_3_arg) {
        if (difop.name.compare(differential_operator_name) == 0){
            // Get result of applying the named differential operator
            Field3D result = difop.func(a, f, flow_ylow);
            dump[outname] = result;
            dump[outname].setAttributes({
                  {"operator", difop.name},
              });
            // Get expected result from input file
            Field3D expected_result{mesh};
            mesh->get(expected_result, expectedname, 0.0, false);
            dump[expectedname] = expected_result;
            // dump diagnostics
            dump[outname_flow_ylow] = flow_ylow;
        }
      }
      for (const auto& difop: differential_operators_4_arg) {
        if (difop.name.compare(differential_operator_name) == 0){
            // Get result of applying the named differential operator
            Field3D result = difop.func(a, f, flow_xlow, flow_ylow);
            dump[outname] = result;
            dump[outname].setAttributes({
                  {"operator", difop.name},
              });
            // Get expected result from input file
            Field3D expected_result{mesh};
            mesh->get(expected_result, expectedname, 0.0, false);
            dump[expectedname] = expected_result;
            // dump diagnostics
            dump[outname_flow_xlow] = flow_xlow;
            dump[outname_flow_ylow] = flow_ylow;
        }
      }
  }
  //Field3D result = FV::Div_a_Grad_perp(a, f);
  //dump["result"] = result;

  //Field3D result_nonorthog = Div_a_Grad_perp_nonorthog(a, f);
  //dump["result_nonorthog"] = result_nonorthog;

  mesh->outputVars(dump);

  std::string outname = fmt::format(
      "{}/BOUT.{}.nc",
      Options::root()["datadir"].withDefault<std::string>("data"), BoutComm::rank());
  
  bout::OptionsIO::create(outname)->write(dump);

  BoutFinalise();

  return 0;
}

#include "Eclipse.hpp"

// Code to generate tensorial contraction code for fic-MRCI

int main()
{

  Eclipse eclipse;
  //  eclipse.generate_fock_vector();
  //  eclipse.generate_fock_diag();
  //  eclipse.generate_overlap();
  eclipse.generate_diag();
  eclipse.generate_g_e();
  eclipse.generate_e_g();
  eclipse.generate_e_e();

}

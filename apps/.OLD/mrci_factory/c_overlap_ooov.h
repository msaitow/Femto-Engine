

#ifndef C_OVERLAP_OOOV_H
#define C_OVERLAP_OOOV_H

// #include <tensor/tensor.h>                  
// #include <sci/hint/hintmo/hintmo.h>         
// #include <sci/ctnew2/ctclass_input.h>       
// #include <sci/ctnew2/ctclass_symblock.h>    
// #include <sci/ctnew2/ctclass_rdmpack.h>     
// #include <sci/ctnew2/ctclass_bareamppack.h> 
                                               
extern "C"{                                  
                                               
                                               
//  __/\\\\\\\\\\\\\\\____________________________________________________________________                                   
//   _\/\\\///////////_____________________________________________________________________                                             
//    _\/\\\_______________________________________________________/\\\_____________________                                         
//     _\/\\\\\\\\\\\__________/\\\\\\\\______/\\\\\__/\\\\\_____/\\\\\\\\\\\______/\\\\\____ 
//      _\/\\\///////_________/\\\/////\\\___/\\\///\\\\\///\\\__\////\\\////_____/\\\///\\\__               
//       _\/\\\_______________/\\\\\\\\\\\___\/\\\_\//\\\__\/\\\_____\/\\\________/\\\__\//\\\_       
//        _\/\\\______________\//\\///////____\/\\\__\/\\\__\/\\\_____\/\\\_/\\___\//\\\__/\\\__            
//         _\/\\\_______________\//\\\\\\\\\\__\/\\\__\/\\\__\/\\\_____\//\\\\\_____\///\\\\\/___    
//          _\///_________________\//////////___\///___\///___\///_______\/////________\/////_____                                   

void FC_FUNC(g_if_overlap_ooov_no0_x0, G_IF_OVERLAP_OOOV_NO0_X0)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const O2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

void FC_FUNC(g_if_overlap_ooov_no0_x1, G_IF_OVERLAP_OOOV_NO0_X1)
  (const FC_INT &sa, const FC_INT &ia, 
   const double * const T2, const double * const O2, 
   const FC_INT &nir, const FC_INT * const nsym,  const FC_INT * const psym);

      
 }     
       
       
 #endif
       
       
 
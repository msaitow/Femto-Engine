-------------------------------------------------------------------------------------
   
                :::::::::: ::::::::::   :::   ::: ::::::::::: :::::::: 
               :+:        :+:         :+:+: :+:+:    :+:    :+:    :+: 
              +:+        +:+        +:+ +:+:+ +:+   +:+    +:+    +:+  
             :#::+::#   +#++:++#   +#+  +:+  +#+   +#+    +#+    +:+   
            +#+        +#+        +#+       +#+   +#+    +#+    +#+    
           #+#        #+#        #+#       #+#   #+#    #+#    #+#     
          ###        ########## ###       ###   ###     ########       

      FEMTO :: An Integrated Toolset For The Automated Tensor Generation

-------------------------------------------------------------------------------------

# Program Overview

**Femto is a tensor generator designed to be a powerful means to develop and evaluate the state-of-the-art electronic structure ansatz.** The object-oriented programming is heavily utilized throughout the code. The main component of Femto is an expression generator, which derives the explicit working equations of the many-body ansatz, and a code generator. The generated tensor contraction codes are based on the hybrid parallelism of message-passing interface (MPI) and global-arrays (GA) toolkit and can be interfaced to *Orz* quantum chemical program package. The main object of this package is the ***libfemto*** in which the classes to represent the quantum chemical expressions in the second-quantized language are implemented.

Due to the fact that Femto is implemented by using the C++ programming language, the expression/code generation by this package is fairly fast in comparison with the other packages written in the scripting languages such as *Second-Quantized Algebra* package developed by Eric Neuscamman.

## Basic structure of the package:

- **Femto/core** - The core component for the expression generation and the implementation of the basis classes that constitute the second-quantized ansatz and the tensor contraction equations. In addition, this component includes the functions for the normal ordering of the second-quantized operators and the like-terms combining. Other utility functions for constructing the core Fock matrices and several other purposes are also provided. Especially, it is notable that the specific data structure of the *SQterm* object, enabling to store the information about the tensor structure on the basis of the *Orz* quantum chemistry program package. The main objects provided in this module are as follows:
	- **Femto::Core::SQindex** - The index object that represents the molecular orbital index such as *core*, *active* and *virtual*. It can also be specified as whether there is a summation over the associated orbital index and whether the index is an *external* index (where "external" represents that the index is the loading index of such out-of-core tensors as the electron-repulsion integral and the double excitation amplitudes in the *Orz* implementation). In the code generation process (carried out in the **Femto::reactor::SQreaktor** class), the dependence of the external indices on a stream of the binary contractions are analyzed to determine the best loop-structure.

	- **Femto::Core::SQtensor** - The tensor class that represents all the tensorial quantities including the one-, two-body integrals, the amplitudes, and the various second-quantized operators. The operators are defined as the *incommutative* tensor while the usual tensor quantity is initiated as a *commutative* one. The notation of the tensor such as *Mulliken* or *Dirac* can also be specified. 

- And so on ... (***uncompleted***)

## Description

This documentation isn't yet completed as mentioned above and as mentioned below, this is a tailored version; the capability is somewhat limited in comparison with the working version. However, it can generate (1) the first and (2) the second working versions of the DMRG-MRCI codes. The first version of the code was used for the original publication [J. Chem. Phys. 139, 044118 (2013)] while the second one was for the application to the water-oxidation reaction [Phys. Chem. Chem. Phys. doi:10.1039/c3cp55225j]. The body of Femto is a library named *libfemto.a* and C++ headerfiles located at Femto/include. The user of Femto can add new functionalities and compose the executable on the basis of *libfemto.a*. These are the actual applications:

- (1) **Femto/apps/new_mrci_factory**: This is a heart of formula/code generator used for the first publication of the DMRG-MRCI. The generated code heavily uses the Message-Passing Interface (MPI) and level3 BLAS subroutines. By interfacing to the DMRG code, the large-scale MRCI problem with up to approx. 24 active and 300 MOs is tractable.

- (2) **Femto/apps/WorldEngine**: This bunch of code generates the world of quantum chemists, i.e. it generates any tensor contraction code once the appropriate *tensor files* are put in *TensorFiles/*. The example tensor files found in *TensorFiles/* are for *FIC-MRCI*, *Canonical Transformation with Singles and Doubles (CTSD)* and *CASPT2* models. All of them are generated by a newly developed expression generator, which is far more powerful than the expression generator in *Femto* itself. The generated tensor contraction code also uses MPI parallelisms and the level 3 BLAS subroutines. Because worldEngine is based on a more effective and generalized code generation algorithm, the resultant tensor contraction code can use up to approximately 32 active and 500 virtual MOs.

## License

This is a tailored version of Femto generator for exhibition. It is released under the open-source license as long as publications that use, or is based on this code cite the appropriate references:

  [1] M. Saitow, Y.Kurashige and T. Yanai, J. Chem. Phys. 139, 044118 (2013).
  
  [2] M. Saitow, Femto :: An Integrated Toolset for the Automated Tensor Generation.

  

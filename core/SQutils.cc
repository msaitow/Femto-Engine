//
//  SQutils.cc
//  
//
//  Created by Masaaki Saitow on 12/07/01.
//  Copyright (c) 2012 Masaaki Saitow. All rights reserved.
//

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <Femto.hpp>

using namespace std;

namespace Femto {

  // *********************************************************
  // Returns n-tuple
  // *********************************************************
  IIvector makePermutations(const int n)
  {
    if(n == 1){
      IIvector outList;
      Ivector temp;
      temp.push_back(0);
      outList.push_back(temp);
      return outList;
    }
    IIvector outList;
    for(size_t i = 0;i < (size_t)n;++i){
      IIvector temp = makePermutations(n-1);
      for(size_t j = 0;j < temp.size();++j){
        for(size_t k = 0;k < temp[j].size();++k){
          if(temp[j][k] >= i) ++temp[j][k];
	} // End k
        outList.push_back(temp[j]); 
        // Looks a little bia tricky, but at least, it works ....
	IIvector::iterator it = outList.end(); --it;
        it->insert(it->begin(), 1, i);
      } // End j
    } // End i
    return outList;
  }

  // *********************************************************
  // Returns n-tuple, composed of elements of (0,1,,...,p) 
  // with redundancy
  // *********************************************************
  IIvector makeTuples1(int n, Ivector &inList)
  {
    IIvector outList;
    int p = (int)inList.size();

    if(n == 0) {
      Ivector temp;
      outList.push_back(temp);
      return outList;
    }
    if(n == 1) {
      for(size_t i = 0;i < p;++i) {
	Ivector temp;
        temp.push_back(inList[i]);
        outList.push_back(temp);
      } // End i      
      return outList;
    } // End if
    if(n == p) {
      Ivector temp;
      for(size_t i = 0;i < p;++i) temp.push_back((int)i);
      outList.push_back(temp);
      return outList;
    } // End i
    if(p < n){
      cout << "Size of p should be larger than that of n" << endl;
      abort();
    } // End if

    Ivector tempList;
    for(size_t i = 0;i < p;++i) tempList.push_back(inList[i]);
    for(size_t i = 0;i < tempList.size();++i){
      Ivector temp2;
      for(size_t t = 0;t < inList.size();++t){
        if(inList[i]!=inList[t]) temp2.push_back(inList[t]);
      }
      IIvector subList = makeTuples1(n-1, temp2);
      for(size_t j = 0;j < subList.size();++j){
        Ivector temp2;
        temp2.push_back(tempList[i]);
        outList.push_back(temp2);
        for(size_t k = 0;k < subList[j].size();++k){
	  IIvector::iterator it = outList.end(); --it;
          it->push_back(subList[j][k]);
	} // End k
      } // End j
    } // End i

    return outList;
  }

  // *********************************************************
  // Returns n-tuple, composed of elements of (0,1,,...,p)
  // without redundancy
  // *********************************************************
  IIvector makeTuples2(int n, Ivector &inList)
  {
    size_t p = inList.size();
    IIvector outList;
    if(n == 0) {
      Ivector temp;
      outList.push_back(temp);
      return outList;
    }
    if(n == (int)p) {
      Ivector temp;
      for(size_t i = 0;i < p;++i) temp.push_back(inList[i]);
      outList.push_back(temp);
      return outList;
    } // End if
    if(n == 1){
      for(size_t i = 0;i < p;++i) {
	Ivector temp;
        temp.push_back(inList[i]);
        outList.push_back(temp);
      } // End i
      return outList;
    } // End if
    if((int)p < n){
      cout << "Size of p should be larger than that of n" << endl;
      abort();
    } // End if
    
    Ivector tempList;
    for(size_t i = 0;i < p-n+1;++i) tempList.push_back(inList[i]);
    for(size_t i = 0;i < tempList.size();++i){
      Ivector temp2;
      for(size_t t = i+1;t < inList.size();++t) temp2.push_back(inList[t]);
      IIvector subList = makeTuples2(n-1, temp2);
      for(size_t j = 0;j < subList.size();++j){
        Ivector temp2;
        temp2.push_back(tempList[i]);
        outList.push_back(temp2);
        for(size_t k = 0;k < subList[j].size();++k){
	  IIvector::iterator it = outList.end(); --it;
          it->push_back(subList[j][k]);
	} // End k
      } // End j
    } // End i

    return outList;
  }

  // *********************************************************
  // Returns n-tuple, composed of elements of (0,1,,...,p) 
  // with partial redundancy
  // *********************************************************
  IIvector makeTuples3(int n, int order, Ivector &inList)
  {
    IIvector outList;
    int p = (int)inList.size();

    if(p%2 != 0){
      cout << "Function makeTupleRDM is designed only for RDM class. So, inList.size()" << endl;
      cout << "should be an even number." << endl;
      abort();
    }

    if(n == 0) return outList;
    if(n == 1) {
      for(size_t i = 0;i < p;++i) {
	Ivector temp;
        temp.push_back(inList[i]);
        outList.push_back(temp);
      } // End i      
      return outList;
    } // End if
    if(n == p) {
      Ivector temp;
      for(size_t i = 0;i < p;++i) temp.push_back((int)i);
      outList.push_back(temp);
      return outList;
    } // End i
    if(p < n){
      cout << "Size of p should be larger than that of n" << endl;
      abort();
    } // End if

    Ivector tempList;
    for(size_t i = 0;i < p;++i) tempList.push_back(inList[i]);
    for(size_t i = 0;i < tempList.size();++i){
      Ivector temp2;
      for(size_t t = 0;t < inList.size();++t){
        if(inList[i]!=inList[t] and (inList[i]>=order ? inList[i]-order : inList[i]+order)!=inList[t]) temp2.push_back(inList[t]);
      }
      IIvector subList = makeTuples3(n-1, order, temp2);
      for(size_t j = 0;j < subList.size();++j){
        Ivector temp2;
        temp2.push_back(tempList[i]);
        outList.push_back(temp2);
        for(size_t k = 0;k < subList[j].size();++k){
	  IIvector::iterator it = outList.end(); --it;
          it->push_back(subList[j][k]);
	} // End k
      } // End j
    } // End i

    return outList;
  }

  // *********************************************************
  // Returns the number of permutations between two Ivectors
  // *********************************************************
  int get_num_perms(vector<int> &ti, vector<int> &bi)
  {
    if(ti.size()!=bi.size()) abort();
    typedef pair<int, int> p_int;
    vector<p_int> x;
    for(size_t i = 0;i < ti.size();++i){
      p_int temp;
      temp.first = ti[i];
      temp.second = bi[i];
      x.push_back(temp);
    } // End i
    sort(x.begin(), x.end(), SecGreat());
    vector<int> y;
    for(size_t i = 0;i < x.size();++i) y.push_back(x[i].first);

    int n_perms = 0;
    for(size_t i = 0;i < y.size();){
      if(y[i] != (int)i){
        int t = y[i];
        y[i]  = y[t];
        y[t]  = t;
        ++n_perms; 
      } // End if
      else ++i;
    } // End i
    return n_perms;
  }

  // *********************************************************
  // Returns factorial up to n
  // *********************************************************
  int fact(const int n){
    if(n==0) return 1;
    int retval = 1;
    for(int i = 0;i < n;++i) retval *= i + 1;
    return retval;
  }

  // *********************************************************
  // Returns symmetry vector for one-body integrals
  // *********************************************************
  Symmetry h1_symm()
  {
    Symmetry h1_symm;

    Femto::Ivector h0; // {0,1}
    h0.push_back(0);
    h0.push_back(1);

    Femto::Ivector h1; // {1,0}
    h1.push_back(1);
    h1.push_back(0);
 
    h1_symm.first.push_back(h0);
    h1_symm.first.push_back(h1);
    h1_symm.second.push_back(1);
    h1_symm.second.push_back(1);

    return h1_symm;
  }

  // *********************************************************
  // Returns symmetry vector for ERI
  // *********************************************************
  Symmetry h2_symm()
  {
    Symmetry V2_symm;

    Femto::Ivector S0; // {0,1,2,3}
    S0.push_back(0);
    S0.push_back(1); 
    S0.push_back(2); 
    S0.push_back(3); 

    Femto::Ivector S1; // {2,1,0,3};
    S1.push_back(2);
    S1.push_back(1); 
    S1.push_back(0); 
    S1.push_back(3); 

    Femto::Ivector S2; // {0,3,2,1};
    S2.push_back(0);
    S2.push_back(3); 
    S2.push_back(2); 
    S2.push_back(1); 

    Femto::Ivector S3; // {1,0,3,2};
    S3.push_back(1);
    S3.push_back(0); 
    S3.push_back(3); 
    S3.push_back(2); 

    V2_symm.first.push_back(S0);
    V2_symm.first.push_back(S1);
    V2_symm.first.push_back(S2);
    V2_symm.first.push_back(S3);

    V2_symm.second.push_back(1);
    V2_symm.second.push_back(1);
    V2_symm.second.push_back(1);
    V2_symm.second.push_back(1);

    return V2_symm;
  }

  // *********************************************************
  // Returns symmetry vector for ERI in Mulliken notation
  // *********************************************************
  Symmetry h2_symmM()
  {
    Symmetry V2_symm;

    Femto::Ivector S0; // {0,1,2,3}
    S0.push_back(0);
    S0.push_back(1); 
    S0.push_back(2); 
    S0.push_back(3); 

    Femto::Ivector S1; // {1,0,2,3};
    S1.push_back(1);
    S1.push_back(0); 
    S1.push_back(2); 
    S1.push_back(3); 

    Femto::Ivector S2; // {0,1,3,2};
    S2.push_back(0);
    S2.push_back(1); 
    S2.push_back(3); 
    S2.push_back(2); 

    Femto::Ivector S3; // {2,3,0,1};
    S3.push_back(2);
    S3.push_back(3); 
    S3.push_back(0); 
    S3.push_back(1); 

    V2_symm.first.push_back(S0);
    V2_symm.first.push_back(S1);
    V2_symm.first.push_back(S2);
    V2_symm.first.push_back(S3);

    V2_symm.second.push_back(1);
    V2_symm.second.push_back(1);
    V2_symm.second.push_back(1);
    V2_symm.second.push_back(1);

    return V2_symm;
  }

  // *********************************************************
  // Returns unit symmetry for 2-index tensor
  // *********************************************************
  Symmetry u2_symm()
  {
    Symmetry u2_symm;
    Femto::Ivector S2_0;
    S2_0.push_back(0);
    S2_0.push_back(1);
    u2_symm.first.push_back(S2_0);
    u2_symm.second.push_back(1);

    return u2_symm;
  }

  // *********************************************************
  // Returns unit symmetry for 4-index tensor
  // *********************************************************
  Symmetry u4_symm()
  {
    Symmetry u4_symm;
    Femto::Ivector S2_0;
    S2_0.push_back(0);
    S2_0.push_back(1);
    S2_0.push_back(2);
    S2_0.push_back(3);
    u4_symm.first.push_back(S2_0);
    u4_symm.second.push_back(1);

    return u4_symm;
  }

  // *********************************************************
  // Returns unit symmetry for n-rank tensor
  // *********************************************************
  Symmetry uni_symm(int n)
  {
    Symmetry un_symm;
    Femto::Ivector S2_0;
    for(int i = 0;i < n;++i) S2_0.push_back(i);
    un_symm.first.push_back(S2_0);
    un_symm.second.push_back(1);

    return un_symm;
  }

  // *********************************************************
  // Returns symmetry for t2 amplitude
  // *********************************************************
  Symmetry t2_symm()
  {
    Symmetry t2_symm;
    Femto::Ivector S2_0; // {0,1,2,3}
    S2_0.push_back(0);
    S2_0.push_back(1);
    S2_0.push_back(2);
    S2_0.push_back(3);

    Femto::Ivector S2_1; // {1.0,3,2}
    S2_1.push_back(1);
    S2_1.push_back(0);
    S2_1.push_back(3);
    S2_1.push_back(2);

    t2_symm.first.push_back(S2_0);
    t2_symm.first.push_back(S2_1);
    t2_symm.second.push_back(1);
    t2_symm.second.push_back(1);

    return t2_symm;
  }

  // *********************************************************
  // Returns current date and time
  // *********************************************************
  string Femto_date()
  {
    time_t timer;
    time(&timer);
    return string(ctime(&timer));
  }

  // *********************************************************
  // Returns a logo
  // *********************************************************
  string Femto_logo(const string s)
  {
    // All these logos were generated from :
    // http://patorjk.com/software/taag/#p=display&f=Graffiti&t=Type%20Something%20
    srand((unsigned)time(NULL));
    int num = rand() % 13;
    string retval("");
    if(num == 0){
      retval += s + " ___________                __               \n" ;
      retval += s + " \\_   _____/____    _____ _/  |_  ____      \n" ;
      retval += s + "  |    __)_/ __ \\  /     \\\\   __\\/  _ \\ \n" ;
      retval += s + "  |     \\ \\  ___/ |  Y Y  \\|  | (  <_> )  \n" ;
      retval += s + "  \\___  /  \\___  >|__|_|  /|__|  \\____/   \n" ;
      retval += s + "      \\/       \\/       \\/                \n" ;
    } // Graffiti
    else if(num == 1){
      retval += s + "     ______                  __           \n" ;
      retval += s + "    / ____/___   ____ ___   / /_ ____     \n" ;
      retval += s + "   / /_   / _ \\ / __ `__ \\ / __// __ \\ \n" ;
      retval += s + "  / __/  /  __// / / / / // /_ / /_/ /    \n" ;
      retval += s + " /_/     \\___//_/ /_/ /_/ \\__/ \\____/  \n" ;
    } // Slant
    else if(num == 2){
      retval += s + " __/\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\____________________________________________________________________                                   \n" ;   
      retval += s + "  _\\/\\\\\\///////////_____________________________________________________________________                                             \n" ;    
      retval += s + "   _\\/\\\\\\_______________________________________________________/\\\\\\_____________________                                         \n" ;        
      retval += s + "    _\\/\\\\\\\\\\\\\\\\\\\\\\__________/\\\\\\\\\\\\\\\\______/\\\\\\\\\\__/\\\\\\\\\\_____/\\\\\\\\\\\\\\\\\\\\\\______/\\\\\\\\\\____ \n" ;          
      retval += s + "     _\\/\\\\\\///////_________/\\\\\\/////\\\\\\___/\\\\\\///\\\\\\\\\\///\\\\\\__\\////\\\\\\////_____/\\\\\\///\\\\\\__               \n" ;            
      retval += s + "      _\\/\\\\\\_______________/\\\\\\\\\\\\\\\\\\\\\\___\\/\\\\\\_\\//\\\\\\__\\/\\\\\\_____\\/\\\\\\________/\\\\\\__\\//\\\\\\_       \n" ;         
      retval += s + "       _\\/\\\\\\______________\\//\\\\///////____\\/\\\\\\__\\/\\\\\\__\\/\\\\\\_____\\/\\\\\\_/\\\\___\\//\\\\\\__/\\\\\\__            \n" ;            
      retval += s + "        _\\/\\\\\\_______________\\//\\\\\\\\\\\\\\\\\\\\__\\/\\\\\\__\\/\\\\\\__\\/\\\\\\_____\\//\\\\\\\\\\_____\\///\\\\\\\\\\/___    \n" ;              
      retval += s + "         _\\///_________________\\//////////___\\///___\\///___\\///_______\\/////________\\/////_____                                   \n" ;         
    } // Slant Relief
    else if(num == 3){
      retval += s + "   o__ __o__/_                            o                         \n" ;
      retval += s + "  <|    v                                <|>                        \n" ;
      retval += s + "  < >                                    < >                        \n" ;
      retval += s + "   |         o__  __o   \\o__ __o__ __o    |        o__ __o         \n" ;
      retval += s + "   o__/_    /v      |>   |     |     |>   o__/_   /v     v\\        \n" ;
      retval += s + "   |       />      //   / \\   / \\   / \\   |      />       <\\    \n" ;
      retval += s + "  <o>      \\o    o/     \\o/   \\o/   \\o/   |      \\         /   \n" ;
      retval += s + "   |        v\\  /v __o   |     |     |    o       o       o        \n" ;
      retval += s + "  / \\        <\\/> __/>  / \\   / \\   / \\   <\\__    <\\__ __/>  \n" ;
    } // Acrobatic
    else if(num == 4){
      retval += s + "  .----------------.  .----------------.  .----------------.  .----------------.  .----------------.       \n" ;
      retval += s + " | .--------------. || .--------------. || .--------------. || .--------------. || .--------------. |      \n" ;
      retval += s + " | |  _________   | || |  _________   | || | ____    ____ | || |  _________   | || |     ____     | |      \n" ;
      retval += s + " | | |_   ___  |  | || | |_   ___  |  | || ||_   \\  /   _|| || | |  _   _  |  | || |   .'    `.   | |     \n" ;
      retval += s + " | |   | |_  \\_|  | || |   | |_  \\_|  | || |  |   \\/   |  | || | |_/ | | \\_|  | || |  /  .--.  \\  | | \n" ;
      retval += s + " | |   |  _|      | || |   |  _|  _   | || |  | |\\  /| |  | || |     | |      | || |  | |    | |  | |     \n" ;
      retval += s + " | |  _| |_       | || |  _| |___/ |  | || | _| |_\\/_| |_ | || |    _| |_     | || |  \\  `--'  /  | |    \n" ;
      retval += s + " | | |_____|      | || | |_________|  | || ||_____||_____|| || |   |_____|    | || |   `.____.'   | |      \n" ;
      retval += s + " | |              | || |              | || |              | || |              | || |              | |      \n" ;
      retval += s + " | '--------------' || '--------------' || '--------------' || '--------------' || '--------------' |      \n" ;
      retval += s + "  '----------------'  '----------------'  '----------------'  '----------------'  '----------------'       \n" ;                                                            
    }
    else if(num == 5){
      retval += s + "                                                             \n" ;
      retval += s + "  _______________                                  ______    \n" ;
      retval += s + " |          |                 .'. .`. `````|`````.~      ~.  \n" ;
      retval += s + " |______    |______         .'   `   `.    |    |          | \n" ; 
      retval += s + " |          |             .'           `.  |    |          | \n" ;
      retval += s + " |          |___________.'               `.|     `.______.'  \n" ;
      retval += s + "                                                             \n" ;          
    }
    else if(num == 6){
      retval += s + "       :::::::::: ::::::::::   :::   ::: ::::::::::: :::::::: \n" ;
      retval += s + "      :+:        :+:         :+:+: :+:+:    :+:    :+:    :+: \n" ;
      retval += s + "     +:+        +:+        +:+ +:+:+ +:+   +:+    +:+    +:+  \n" ;
      retval += s + "    :#::+::#   +#++:++#   +#+  +:+  +#+   +#+    +#+    +:+   \n" ;   
      retval += s + "   +#+        +#+        +#+       +#+   +#+    +#+    +#+    \n" ;        
      retval += s + "  #+#        #+#        #+#       #+#   #+#    #+#    #+#     \n" ;         
      retval += s + " ###        ########## ###       ###   ###     ########       \n" ; 
    }
    else if(num == 7){
      retval += s + " 8 8888888888   8 8888888888            ,8.       ,8.    8888888 8888888888 ,o888888o.     \n" ; 
      retval += s + " 8 8888         8 8888                 ,888.     ,888.         8 8888    . 8888     `88.   \n" ; 
      retval += s + " 8 8888         8 8888                .`8888.   .`8888.        8 8888   ,8 8888       `8b  \n" ;  
      retval += s + " 8 8888         8 8888               ,8.`8888. ,8.`8888.       8 8888   88 8888        `8b \n" ;   
      retval += s + " 8 888888888888 8 888888888888      ,8'8.`8888,8^8.`8888.      8 8888   88 8888         88 \n" ;   
      retval += s + " 8 8888         8 8888             ,8' `8.`8888' `8.`8888.     8 8888   88 8888         88 \n" ;   
      retval += s + " 8 8888         8 8888            ,8'   `8.`88'   `8.`8888.    8 8888   88 8888        ,8P \n" ;    
      retval += s + " 8 8888         8 8888           ,8'     `8.`'     `8.`8888.   8 8888   `8 8888       ,8P  \n" ;     
      retval += s + " 8 8888         8 8888          ,8'       `8        `8.`8888.  8 8888    ` 8888     ,88'   \n" ;          
      retval += s + " 8 8888         8 888888888888 ,8'         `         `8.`8888. 8 8888       `8888888P'     \n" ;          
    }
    else if(num == 8){
      retval += s + " 8888888888                     888                  \n" ;       
      retval += s + " 888                            888                  \n" ;           
      retval += s + " 888                            888                  \n" ;             
      retval += s + " 8888888  .d88b.  88888b.d88b.  888888  .d88b.       \n" ;                   
      retval += s + " 888     d8P  Y8b 888 \"888 \"88b 888    d88\"\"88b  \n" ;               
      retval += s + " 888     88888888 888  888  888 888    888  888      \n" ;               
      retval += s + " 888     Y8b.     888  888  888 Y88b.  Y88..88P      \n" ;       
      retval += s + " 888      \"Y8888  888  888  888  \"Y888  \"Y88P\"   \n" ;         
    }
    else if(num == 9){
      retval += s + " `MMMMMMM                                         \n" ;    
      retval += s + "  MM    \\                         /              \n" ;  
      retval += s + "  MM       ____  ___  __    __   /M      _____    \n" ;       
      retval += s + "  MM   ,  6MMMMb `MM 6MMb  6MMb /MMMMM  6MMMMMb   \n" ;    
      retval += s + "  MMMMMM 6M'  `Mb MM69 `MM69 `Mb MM    6M'   `Mb  \n" ;        
      retval += s + "  MM   ` MM    MM MM'   MM'   MM MM    MM     MM  \n" ;          
      retval += s + "  MM     MMMMMMMM MM    MM    MM MM    MM     MM  \n" ;           
      retval += s + "  MM     MM       MM    MM    MM MM    MM     MM  \n" ;          
      retval += s + "  MM     YM    d9 MM    MM    MM YM.  ,YM.   ,M9  \n" ;   
      retval += s + " _MM_     YMMMM9 _MM_  _MM_  _MM_ YMMM9 YMMMMM9   \n" ;      
    }
    else if(num == 10){
      retval += s + " `7MM\"\"\"YMM                         mm               \n" ;   
      retval += s + "   MM    `7                         MM                  \n" ;           
      retval += s + "   MM   d  .gP\"Ya `7MMpMMMb.pMMMb.mmMMmm ,pW\"Wq.      \n" ;            
      retval += s + "   MM\"\"MM ,M\'   Yb  MM    MM    MM  MM  6W\'   `Wb   \n" ;            
      retval += s + "   MM   Y 8M\"\"\"\"\"\"  MM    MM    MM  MM  8M     M8 \n" ;             
      retval += s + "   MM     YM.    ,  MM    MM    MM  MM  YA.   ,A9       \n" ;                
      retval += s + " .JMML.    `Mbmmd'.JMML  JMML  JMML.`Mbmo`Ybmd9'        \n" ;                 
    }
    else if(num == 11){
      retval += s + "      #                #########      #     #   # \n" ;
      retval += s + " ########## ##########         #   #######  #   # \n" ;
      retval += s + "     #    #         #          #    # #     #   # \n" ;
      retval += s + "     #    #        #   ########     # #     #   # \n" ;
      retval += s + "    #     #     # #           #  ##########    #  \n" ;   
      retval += s + "   #   # #       #            #       #       #   \n" ;     
      retval += s + "  #     #         #    ########       #     ##    \n" ;     
    }
    else if(num == 12){
      retval += s + "     _/_/_/_/                            _/             \n" ;                
      retval += s + "    _/        _/_/    _/_/_/  _/_/    _/_/_/_/    _/_/  \n" ;                  
      retval += s + "   _/_/_/  _/_/_/_/  _/    _/    _/    _/      _/    _/ \n" ;                 
      retval += s + "  _/      _/        _/    _/    _/    _/      _/    _/  \n" ;                    
      retval += s + " _/        _/_/_/  _/    _/    _/      _/_/    _/_/     \n" ;                      
    }
    return retval;
  }

} //Femto::


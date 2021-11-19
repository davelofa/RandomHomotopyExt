 /* Copyright (c) 2021
   Davide Lofano
   lofano@math.tu-berlin.de
   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
--------------------------------------------------------------------------------
   $Project: polymake $$Id$
*/


#include "polymake/Rational.h"
#include "polymake/topaz/complex_tools.h"
#include "polymake/Matrix.h"
#include "polymake/RandomSubset.h"

namespace polymake { namespace topaz {
    
    
   
   BigObject R_H2(BigObject p_in, const Int d, const Int move, OptionSet options) {
        
       const RandomSeed seed(options["seed"]);
       UniformlyRandom<long> random_source(seed);
          
       FacetList Facets =p_in.give("FACETS");
       
       Set<Int> add_facet;
       Int mas_fac=0;
       
       const Lattice<BasicDecoration>& HD = p_in.give("HASSE_DIAGRAM");
       
       const Int dim = HD.rank()-2; 
        
       if (dim < d) { 
           cout << "wrong dimension"<<endl;
       }
       else {
           
     
       
       
       
       //Initialize list of d_faces and possible expansions
       FacetList d_Face;
       FacetList possible_Face;
       

       for (auto k = entire(HD.nodes_of_rank(d+1)); !k.at_end(); ++k) {
         const Int this_index = *k;
         d_Face.insert(HD.face(this_index));
         possible_Face.insert(HD.face(this_index));
         
       }
	    
        
        if (move!=0) {
       
        while ((!possible_Face.empty())& (mas_fac<d+1)) 
        {
            
        //choose a d_faces uniformly at random
            UniformlyRandomRanged<long> random(possible_Face.size(),random_source);
            long r_long(random.get());
                     
            FacetList::const_iterator elem=possible_Face.begin();
            for (Int elem_i=0; elem_i<r_long; ++elem_i) ++elem;
            
            
            Set< Int> this_face=*elem;
            
               //Vertices in a random order
            Array<Int> RVERTS(Facets.n_vertices(), random_permutation(Facets.n_vertices(), seed).begin());
            

            //Check one vertex at a time if I can add it
            
            for (Int i=0; i<Facets.n_vertices(); i++) {
                Set<Int> this_facet=this_face;
                this_facet+=RVERTS[i];
                if (this_facet!=this_face) 
                 {
                     bool hom=true; 
                     Int x=0;
                     BigObject subcomp=call_function("induced_subcomplex",p_in,this_facet);
                     
                     if (!subcomp.give("PURE")) hom=false ;
                    
                     
                     
                              
                     for (FacetList::subset_iterator<pm::Set<Int, pm::operations::cmp> > step=d_Face.findSubsets(this_facet); !step.at_end(); ++step){ 
                         Set<Int> V=*step;
                   
                         ++x;
                }
              

                    if ((x>0) & (x>mas_fac) & (x<d+2) & (hom==true))
                      {
                            mas_fac=x;
                            add_facet=this_facet;
                      }

                    }
                 }
                   

           possible_Face.erase(this_face);
           
        }
        }
   
            if (mas_fac!=0) {
  
                      Facets.insert(add_facet);

                    }
                    
            else {
                
                UniformlyRandomRanged<long> random(d_Face.size(),random_source);
                long r_long(random.get());
                     
                FacetList::const_iterator elem=d_Face.begin();
                for (Int elem_i=0; elem_i<r_long; ++elem_i) ++elem;
            
            
                add_facet=*elem;
                
                add_facet+=Facets.n_vertices();
                
                 Facets.insert(add_facet);  
      
        
            }

       }
           
        BigObject p_out("SimplicialComplex");
        p_out.take("INPUT_FACES") << Facets;
   return p_out;

}


UserFunction4perl("CREDIT none\n\n"
                  "# Add a simplex of dimension d without changing the homology"  
                  "# @param SimplicialComplex complex"
                  "# @param Int d"
                  "# @param Int move"
                  "# @option Int seed"
                  "# @return SimplicialComplex",
                  &R_H2,"R_H2(SimplicialComplex, Int, Int, { seed=> undef })");

} }

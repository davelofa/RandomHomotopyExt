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


#include "polymake/client.h"
#include "polymake/SparseMatrix.h"
#include "polymake/SparseVector.h"
#include "polymake/FacetList.h"
#include "polymake/Array.h"
#include "polymake/Set.h"
#include "polymake/Rational.h"
#include "polymake/topaz/complex_tools.h"
#include "polymake/Matrix.h"
#include "polymake/RandomSubset.h"

namespace polymake { namespace topaz {
    
    
   
   BigObject R_H(BigObject p_in, const Int d, Int n, OptionSet options) {
        
       const RandomSeed seed(options["seed"]);
       UniformlyRandom<long> random_source(seed);
          
       FacetList Facets =p_in.give("FACETS");
       
       const Lattice<BasicDecoration>& HD = p_in.give("HASSE_DIAGRAM");
       
       const Int dim = HD.rank()-2; 
        
       if (dim < d) { 
           cout << "wrong dimension"<<endl;
       }
       else {
           
   
       //Initialize list of Superfree_faces (d_face not contained in any d+1 face) and d_faces
       FacetList d_Face;
       FacetList Superfree_faces;

       for (auto k = entire(HD.nodes_of_rank(d+1)); !k.at_end(); ++k) {
         const Int this_index = *k;
         d_Face.insert(HD.face(this_index));
         if( HD.out_degree(this_index) == 1) {
              Set<Int> faces_this_index=HD.out_adjacent_nodes(this_index);
            const Int first_face = faces_this_index.front();
          if(HD.rank() == HD.rank(first_face)) {
                Superfree_faces.insert(HD.face(this_index));
            }
	      }
       }
	    
	  FacetList d_1_Face;  
      for (auto k = entire(HD.nodes_of_rank(d)); !k.at_end(); ++k) {
         const Int this_index = *k;
         d_1_Face.insert(HD.face(this_index));
      }

       
        while ((!Superfree_faces.empty())& (n>0)) 
        {
            
        //choose a Superfree_faces uniformly at random
            UniformlyRandomRanged<long> random(Superfree_faces.size(),random_source);
            long r_long(random.get());
 
            
            FacetList::const_iterator elem=Superfree_faces.begin();
            for (Int elem_i=0; elem_i<r_long; ++elem_i) ++elem;
            
            Set< Int> this_face=*elem;

               //Vertices in a random order
            Array<Int> RVERTS(Facets.n_vertices(), random_permutation(Facets.n_vertices(), seed).begin());

            Int VERT=-1;
            Int mas_fac=0;
            //Check one vertex at a time if I can add it, choose the one with a lot of faces
            
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
                            VERT=i;
                      }

                    }
                 }
                   

           
     
            if (VERT!=-1) {
                        // Add this_facet to the complex and collapse this_facet,this_face. Update everything.
                        n--;
                        
                        Set<Int> this_facet=this_face;
                        this_facet+=RVERTS[VERT];
                        d_Face.erase(this_face);
                        Facets.eraseSubsets(this_facet);
                        for (auto w=entire(this_face); !w.at_end(); ++w) {
                            Set<Int> f=this_face;
                            f-=*w;
                            f+=RVERTS[VERT];
                            if (d_Face.erase(f)==0) Superfree_faces.insert(f);
                            d_Face.insert(f);
                            Facets.insert(f);
                    }
                }
                        Superfree_faces.erase(this_face);
                
              }
                           

       }
           
        BigObject p_out("SimplicialComplex");
        p_out.take("INPUT_FACES") << Facets;
   return p_out;

}


UserFunction4perl("CREDIT none\n\n"
                  "# Add n simplices of dimension d without changing the homology"  
                  "# @param SimplicialComplex complex"
                  "# @param Int d"
                  "# @param Int n"
                  "# @option Int seed"
                  "# @return SimplicialComplex",
                  &R_H,"R_H(SimplicialComplex, Int, Int, { seed=> undef })");

} }

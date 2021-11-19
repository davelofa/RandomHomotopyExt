 /* Copyright (c) 2021
   Michael Joswig, Frank Lutz, Mimi Tsuruga, Davide Lofano
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
#include "polymake/graph/ShrinkingLattice.h"
#include "polymake/topaz/complex_tools.h"
#include <sys/time.h>
#include "polymake/RandomSubset.h"


namespace polymake { namespace topaz {

using graph::ShrinkingLattice;
using graph::Lattice;
using graph::lattice::BasicDecoration;


/* *@class CompareByHasseDiagram
*   @brief Compares two Ints by comparing the corresponding faces in the Hasse diagram (lex)
*/
class CompareByHasseDiagram
{
   const ShrinkingLattice<BasicDecoration>& HD_;
   const Array<Int>& relabel_;

public:
  CompareByHasseDiagram(const ShrinkingLattice<BasicDecoration>& HD,const Array<Int>& relabel) :
      HD_(HD),
      relabel_(relabel)
  { }

  Set<Int> newlabels(const Set<Int>& s) const
  {
      Set<Int> permuted_vertices;
      for (const auto it : s) {
	 const Int orig_vertex=it;
	 permuted_vertices += relabel_[orig_vertex];
      }
      return permuted_vertices;
  }

  pm::cmp_value operator() (const Int& a, const Int& b) const
  {
     return operations::cmp()( newlabels( HD_.face(a) ),newlabels( HD_.face(b) ) );
  }
};


// free_faces
//   Initialize free_face_list for new max_d. Corresponds to the list "pairs"
//   from original GAP code by Lutz. Instead of tracking both the free face F
//   (dim=max_d - 1) and the face G (dim = max_d) that F is on the boundary of,
//   we only store F. A free face is a face that is contained in exactly one
//   max_d-dimensional face.
//   Procedure: Start with fresh new free_face_list. Add all faces F of
//              dim=max_d-1 such that F is on the boundary of only facet
//              (of dim=max_d).
// @param HasseDiagram newHD: updated HD, ie, some sublattice of orig_HD
// @param Int max_d: dimension of maximal face in newHD
// @param Set<Int> free_face_list: list of free faces (dim = max_d-1) in newHD
void my_rand_free_faces(const ShrinkingLattice<BasicDecoration>& newHD,
		const Int& max_d,
		Set<Int>& free_face_list)
{

   for (auto n = entire(newHD.nodes_of_rank(max_d)); !n.at_end(); ++n) {
      const Int this_index = *n;
      if( newHD.out_degree(this_index) == 1) {
          Set<Int> faces_this_index=newHD.out_adjacent_nodes(this_index);
          const Int remove_face = faces_this_index.front();
          if(newHD.rank(this_index)+1 == newHD.rank(remove_face)) {
                free_face_list += this_index;
          }
      }
   }
}


// rand_collapse
//   Perform a collapse of specified face and update Hasse Diagram and list of
//   free faces.
//   Procedure: Find remove_face = face of remove_this. Delete nodes remove_thisbool first_removed_face(true);
//   		and remove_face from newHD. Check if boundary faces of
//		remove_face have become free, and if so add to free_face_list.
//		Also remove any faces from  free_face_list that are no longer
//              free faces.
// @param HasseDiagram newHD: ("global" param) will be updated within function
// @param Set<Int> free_face_list: also will be updated here
// @param Int remove_this: the face to begin collapse; will be of dim=max_d-1
void rand_collapse(ShrinkingLattice<BasicDecoration>& newHD, Set<Int>& free_face_list,
	      const Int& remove_this)
{
   Set<Int> faces_of_remove_this=newHD.out_adjacent_nodes(remove_this);

   if (faces_of_remove_this.size() != 1) {
      throw std::runtime_error("random_discrete_morse::collapse: collapsing a non-free face");
   }

   // node of the face of remove_this
   const Int remove_face = faces_of_remove_this.front();

   if(newHD.rank(remove_this)+1 != newHD.rank(remove_face)) {
      throw std::runtime_error("random_discrete_morse::collapse: dimensions of Hasse messed up");
   }

   // keep the nodes of boundary faces of remove_face for later use
   Set<Int> bdy_of_remove_face=newHD.in_adjacent_nodes(remove_face);

   // update free_face_list
   // remember all elements of free_face_list are of dimension max_d-1

   // first remove remove_this from free_face_list
   free_face_list-=remove_this;

   // faces that were on the boundary of remove_face are no longer free
   for (auto s = entire(bdy_of_remove_face); !s.at_end(); ++s) {
      const Int this_bdy_face = *s;
      free_face_list-=this_bdy_face;
   }

   // remove the nodes from the Hasse diagram
   newHD.delete_node(remove_this);
   newHD.delete_node(remove_face);

   // deletion of remove_face may add new free faces
   for (auto it=entire(bdy_of_remove_face); !it.at_end(); ++it) {
      const Int this_index = *it;
      if ( newHD.out_degree(this_index) == 1) {
	 free_face_list+=this_index;
      }
   }

}


FacetList do_collapses(ShrinkingLattice<BasicDecoration> newHD, const pm::SharedRandomState& random_source) {
   
    
  FacetList remaining_facets;   
       
    const Int global_d = newHD.rank()-2;   // needed for "Warning" below
    Int max_d(global_d);		   // dimension of maximum-dim face in newHD
    if ( max_d<1 ) {
        for (auto k = entire(newHD.nodes_of_rank(1)); !k.at_end(); ++k) {
         const Int this_index = *k;
//          cout<<"Facce aggiunte " <<newHD.face(this_index)<<endl;
         remaining_facets.insert(newHD.face(this_index));
               }
   
            return remaining_facets;

        }


   // find free faces of newHD
   // remember elements of free_face_list are of dim=max_d-1
   Set<Int> free_face_list;
   my_rand_free_faces(newHD,max_d,free_face_list);


   while (true) {
      if (!free_face_list.empty()) {
	 // collapse anything that can be collapsed

	 //choose remove_this uniformly at random
	 UniformlyRandomRanged<long> random(free_face_list.size(),random_source);
	 long r_long(random.get());

	 Set<Int>::const_iterator elem=free_face_list.begin();
	 for (Int elem_i=0; elem_i<r_long; ++elem_i) ++elem;
	 const Int remove_this=*elem;

	 if(!newHD.node_exists(remove_this))
	    throw std::runtime_error("random_discrete_morse::rand_discMorse::can't remove this");
        
	 rand_collapse(newHD,free_face_list,remove_this);

      }
      else {
         for (auto k = entire(newHD.nodes_of_rank(max_d+1)); !k.at_end(); ++k) {
         const Int this_index = *k;
         remaining_facets.insert(newHD.face(this_index)); 
         }
	     --max_d;

	    if(max_d>0) {

	       // reinitialize max_face_list and free_face_list
	       my_rand_free_faces(newHD,max_d,free_face_list);
	    }
	    else {
         for (auto k = entire(newHD.nodes_of_rank(1)); !k.at_end(); ++k) {
         const Int this_index = *k;
         remaining_facets.insert(newHD.face(this_index)); 
         }
           break;
        }
	 }
      }
   
return remaining_facets;

}

BigObject collapses(BigObject p_in, OptionSet options) {
    
   const Lattice<BasicDecoration>& orig_HD = p_in.give("HASSE_DIAGRAM");
    
    RandomSeed seed = options["seed"];
    UniformlyRandom<long> random_source(seed);
    
  FacetList New_Facets=do_collapses(orig_HD,random_source);
  
  
//From Hasse diagram to simplicial complex
    BigObject p_out("SimplicialComplex");
    
     
              
  
    p_out.take("INPUT_FACES") << New_Facets;
        
   return p_out;   
}

UserFunction4perl("CREDIT none\n\n"
                  "#Collapse a complex to its core without deleting critical faces"  
                  "# @param SimplicialComplex complex"
                  "# @option Int seed"
                  "# @return SimplicialComplex",
                  &collapses,"collapses(SimplicialComplex, { seed=> undef })");


} }

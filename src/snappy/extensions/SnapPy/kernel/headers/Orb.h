#ifndef _Orb_
#define _Orb_

#include "SnapPea.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

/************************************************************************/
/*                                                                      */
/*                    orb_hyperbolic_structure.c                        */
/*                                                                      */
/************************************************************************/

extern SolutionType orb_find_hyperbolic_structure(
    Triangulation *manifold,
    Boolean        manual);

/************************************************************************/
/*                                                                      */
/*                           orb_interface.c                            */
/*                                                                      */
/************************************************************************/

extern int orb_get_num_singular_edges( Triangulation *manifold);
    
extern void orb_get_singularity_info( Triangulation *manifold,
                                      int            singular_index,
                                      Real           *singular_order,
                                      Real           *inner_product);

extern void orb_set_singularity_info( Triangulation *manifold,
                                      int           singular_index,
                                      Real          singular_order);

/************************************************************************/
/*                                                                      */
/*                           orb_volume.c                               */
/*                                                                      */
/************************************************************************/

extern Real orb_volume(Triangulation *manifold);

SNAPPEA_NAMESPACE_END_SCOPE

#endif

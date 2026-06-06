#ifndef _orb_kernel_prototypes_
#define _orb_kernel_prototypes_

#include "SnapPea.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

/************************************************************************/
/*                                                                      */
/*                         orb_cusp_area.c                              */
/*                                                                      */
/************************************************************************/

extern void orb_normalize_cusps(Triangulation *manifold);

/************************************************************************/
/*                                                                      */
/*                           orb_cusps.c                                */
/*                                                                      */
/************************************************************************/

extern void orb_cusps_fill_incident_singular_edges(Triangulation * manifold);
extern Real orb_compute_orbifold_cusp_euler_characteristic(Cusp * cusp);

/************************************************************************/
/*                                                                      */
/*                    orb_hyperbolic_structure.c                        */
/*                                                                      */
/************************************************************************/

extern Real orb_minor1(GL4RMatrix matrix, int row, int col);

/************************************************************************/
/*                                                                      */
/*                    orb_identify_solution_type.c                      */
/*                                                                      */
/************************************************************************/

extern void orb_identify_solution_type(Triangulation *manifold);
extern Boolean orb_contains_flat_tetrahedra( Triangulation *manifold );
extern Boolean orb_solution_is_degenerate(Triangulation *manifold);

/************************************************************************/
/*                                                                      */
/*                            orb_tilts.c                               */
/*                                                                      */
/************************************************************************/

extern void orb_compute_tilts(Triangulation *manifold);

SNAPPEA_NAMESPACE_END_SCOPE

#endif

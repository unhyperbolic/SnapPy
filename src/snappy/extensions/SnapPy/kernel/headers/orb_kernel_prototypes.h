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
/*                          orb_diagram.c                               */
/*                                                                      */
/************************************************************************/

extern void orb_initialize_diagram(OrbDiagram *);
extern void orb_initialize_diagram_vertex(OrbDiagramVertex *vertex);
extern void orb_initialize_diagram_edge(OrbDiagramEdge * edge);
extern void orb_assign_diagram_arcs(OrbDiagram *);
extern void orb_assign_diagram_links(OrbDiagram *);
extern void orb_add_end_data_to_diagram_vertex(OrbDiagramEndData * data, OrbDiagramVertex * vertex);
extern void orb_free_diagram(OrbDiagram *);
extern OrbGraph * orb_diagram_to_graph(OrbDiagram *);
extern Triangulation * orb_triangulate_diagram_complement(OrbDiagram *, Boolean remove_finite_vertices);



/************************************************************************/
/*                                                                      */
/*                           orb_graph.c                                */
/*                                                                      */
/************************************************************************/

extern void orb_free_graph(OrbGraph *gamma);

/************************************************************************/
/*                                                                      */
/*                      orb_graph_complement.c                          */
/*                                                                      */
/************************************************************************/

Triangulation *orb_triangulate_graph_complement(OrbGraph *gamma, Boolean remove_finite_vertices);

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

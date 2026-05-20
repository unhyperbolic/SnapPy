#include "kernel.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

#define ORB_CUSP_AREA           0.3
#define ORB_CUSP_AREA_EPSILON   1e-8

static void orb_compute_cusp_areas(Triangulation *manifold);
static Real orb_compute_link_area(Tetrahedron *tet, int v);

void orb_normalize_cusps(
    Triangulation *manifold)
{
    orb_compute_cusp_areas(manifold);

    for (Cusp *cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)
        if (cusp->topology == torus_cusp || cusp->topology == Klein_cusp)
        {
            Real scalar = safe_sqrt(cusp->orb_cusp_shape->area / ORB_CUSP_AREA);

            cusp->orb_cusp_shape->inner_product[ultimate] *= scalar * scalar;

            for (EdgeClass *edge = manifold->edge_list_begin.next;
                 edge != &manifold->edge_list_end;
                 edge = edge->next)
            {
                int index = edge->incident_edge_index;
                Tetrahedron *tet = edge->incident_tet;
                Cusp *top_cusp = tet->cusp[one_vertex_at_edge[index]];
                Cusp *bottom_cusp = tet->cusp[other_vertex_at_edge[index]];

                if (cusp == top_cusp)
                    edge->orb_edge_shape->inner_product[ultimate] *= scalar;
                if (cusp == bottom_cusp)
                    edge->orb_edge_shape->inner_product[ultimate] *= scalar;
            }

            for (Tetrahedron *tet = manifold->tet_list_begin.next;
                 tet != &manifold->tet_list_end;
                 tet = tet->next)
                for (int i = 0; i < 4; i++)
                    if (tet->cusp[i] == cusp)
                        tet->orb_tet_shape->orientation_parameter[ultimate] *= scalar;
        }

    for (Tetrahedron *tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)
    {
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                if (i != j)
                    tet->orb_tet_shape->Gram_matrix[i][j]
                        = tet->edge_class[edge_between_vertices[i][j]]
                              ->orb_edge_shape->inner_product[ultimate];
                else
                    tet->orb_tet_shape->Gram_matrix[i][i]
                        = tet->cusp[i]->orb_cusp_shape->inner_product[ultimate];

        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                tet->orb_tet_shape->inverse_Gram_matrix[i][j]
                    = orb_minor1(tet->orb_tet_shape->Gram_matrix, i, j);
    }
}

static void orb_compute_cusp_areas(
    Triangulation *manifold)
{
    for (Cusp *cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)
        cusp->orb_cusp_shape->area = 0.0;

    for (Tetrahedron *tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)
        for (int v = 0; v < 4; v++)
            tet->cusp[v]->orb_cusp_shape->area += orb_compute_link_area(tet, v);
}

static Real orb_compute_link_area(
    Tetrahedron *tet,
    int          v)
{
    Real top, bottom;

    if (tet->orb_tet_shape->orientation_parameter[ultimate] < ORB_CUSP_AREA_EPSILON)
        return 0.0;

    top = -gl4R_determinant(tet->orb_tet_shape->Gram_matrix)
        * gl4R_determinant(tet->orb_tet_shape->Gram_matrix);

    bottom = 2.0;

    for (int i = 0; i < 4; i++)
        if (i != v)
            for (int j = i; j < 4; j++)
                if (j != v)
                {
                    if (i == j)
                        bottom *= tet->orb_tet_shape->inverse_Gram_matrix[i][i];
                    else
                        bottom *= sin(
                            tet->orb_tet_shape->dihedral_angle[ultimate]
                                                            [edge_between_faces[i][j]]);
                }

    return top / bottom;
}

SNAPPEA_NAMESPACE_END_SCOPE

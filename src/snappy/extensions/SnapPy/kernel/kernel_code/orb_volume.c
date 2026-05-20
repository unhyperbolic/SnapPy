#include "kernel.h"
#include "dilog.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

static Complex orb_U(Complex z, Real *angles, Boolean *ok);
static Boolean orb_flat_tet(Tetrahedron *tet);
static Real tetrahedron_volume(Real *angles, Boolean *ok);

Real orb_volume(
    Triangulation *manifold,
    Boolean       *ok)
{
    Real volume = 0.0;

    *ok = TRUE;

    for (Tetrahedron *tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)
    {
        Real angles[6];
        Boolean ok1;
        Real tet_vol;

        for (int i = 0; i < 6; i++)
            angles[i] = tet->orb_tet_shape->dihedral_angle[ultimate][i];

        tet_vol = tetrahedron_volume(angles, &ok1);

        if (!ok1 && !orb_flat_tet(tet))
        {
            *ok = FALSE;
            uFatalError("orb_volume", "orb_volume");
        }

        if (tet->orb_tet_shape->orientation_parameter[ultimate] > 0)
            volume += tet_vol;
        else
            volume -= tet_vol;
    }

    return volume;
}

static Real tetrahedron_volume(
    Real    *angles,
    Boolean *ok)
{
    static const int opposite[] = {5, 4, 3};
    GL4RMatrix G;
    Complex w1, w2, w, z1, z2, bottom;
    Real sqrt_det, real_top;

    *ok = TRUE;

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            G[i][j] = (i == j) ? (Real)1.0 : -cos(angles[edge_between_faces[i][j]]);

    real_top = 0.0;
    for (int i = 0; i < 3; i++)
        real_top -= 2.0 * sin(angles[i]) * sin(angles[opposite[i]]);

    sqrt_det = sqrt(ABS(gl4R_determinant(G)));

    z1.real = real_top;
    z1.imag = 2.0 * sqrt_det;

    z2.real = real_top;
    z2.imag = -2.0 * sqrt_det;

    bottom = Zero;

    for (int i = 0; i < 3; i++)
    {
        w1.real = cos(angles[i]);
        w1.imag = sin(angles[i]);
        w2.real = cos(angles[opposite[i]]);
        w2.imag = sin(angles[opposite[i]]);
        w = complex_mult(w1, w2);
        bottom = complex_plus(bottom, w);
    }

    for (int i = 0; i < 4; i++)
    {
        w = One;
        for (int j = 0; j < 4; j++)
            if (i != j)
            {
                w1.real = cos(angles[edge_between_faces[i][j]]);
                w1.imag = sin(angles[edge_between_faces[i][j]]);
                w = complex_mult(w, w1);
            }
        bottom = complex_plus(bottom, w);
    }

    w = One;
    for (int i = 0; i < 6; i++)
    {
        w1.real = cos(angles[i]);
        w1.imag = sin(angles[i]);
        w = complex_mult(w, w1);
    }

    bottom = complex_plus(bottom, w);

    z1 = complex_div(z1, bottom);
    z2 = complex_div(z2, bottom);

    return complex_minus(orb_U(z1, angles, ok), orb_U(z2, angles, ok)).imag / 2;
}

static Complex orb_U(
    Complex  z,
    Real    *angles,
    Boolean *ok)
{
    static const int opposite[] = {5, 4, 3};
    Complex result = complex_volume_dilog(z), w, w1, w2, dilogw;

    (void)ok;

    for (int i = 0; i < 3; i++)
    {
        w = One;

        for (int j = 0; j < 3; j++)
            if (i != j)
            {
                w1.real = cos(angles[j]);
                w1.imag = sin(angles[j]);
                w2.real = cos(angles[opposite[j]]);
                w2.imag = sin(angles[opposite[j]]);
                w = complex_mult(w, w1);
                w = complex_mult(w, w2);
            }

        w = complex_mult(w, z);
        dilogw = complex_volume_dilog(w);
        result = complex_plus(result, dilogw);
    }

    for (int i = 0; i < 4; i++)
    {
        w = MinusOne;

        for (int j = 0; j < 4; j++)
            if (i != j)
            {
                w1.real = cos(angles[edge_between_vertices[i][j]]);
                w1.imag = sin(angles[edge_between_vertices[i][j]]);
                w = complex_mult(w, w1);
            }

        w = complex_mult(w, z);
        dilogw = complex_volume_dilog(w);
        result = complex_minus(result, dilogw);
    }

    return complex_real_mult(0.5, result);
}

static Boolean orb_flat_tet(
    Tetrahedron *tet)
{
    for (int i = 0; i < 4; i++)
        for (int j = i + 1; j < 4; j++)
        {
            EdgeIndex e1 = edge_between_vertices[i][j];
            EdgeIndex e2 = edge_between_faces[i][j];
            EdgeIndex e3 = edge_between_vertices[i][one_vertex_at_edge[e2]];
            EdgeIndex e4 = edge_between_vertices[i][other_vertex_at_edge[e2]];
            EdgeIndex e5 = edge_between_vertices[j][one_vertex_at_edge[e2]];
            EdgeIndex e6 = edge_between_vertices[j][other_vertex_at_edge[e2]];

            if (ABS(tet->orb_tet_shape->dihedral_angle[ultimate][e1] - PI) < 1e-6
             && ABS(tet->orb_tet_shape->dihedral_angle[ultimate][e2] - PI) < 1e-6
             && ABS(tet->orb_tet_shape->dihedral_angle[ultimate][e3]) < 1e-6
             && ABS(tet->orb_tet_shape->dihedral_angle[ultimate][e4]) < 1e-6
             && ABS(tet->orb_tet_shape->dihedral_angle[ultimate][e5]) < 1e-6
             && ABS(tet->orb_tet_shape->dihedral_angle[ultimate][e6]) < 1e-6)
                return TRUE;
        }

    return FALSE;
}

SNAPPEA_NAMESPACE_END_SCOPE

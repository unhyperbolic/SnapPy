/**
 * orb_diagram.h
 *
 * Data structures (prefixed by OrbDiagram) to encode a knotted graph as
 * planar diagram with crossings (generalizes a knot diagram).
 *
 * Functions to convert a OrbDiagram to a OrbGraph and triangulate the complement.
 *
 * A OrbDiagram has an embedding into the plane and its edges can cross.
 * The conversion to OrbGraph turns each crossing into a four-valent vertex
 * results in a fat graph.
 *
 */

#ifndef _orb_diagram_
#define _orb_diagram_

#include "SnapPea.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

typedef struct OrbDiagramEndData OrbDiagramEndData;
typedef struct OrbDiagramEdge OrbDiagramEdge;
typedef struct OrbDiagramVertex OrbDiagramVertex;
typedef struct OrbDiagramCrossing OrbDiagramCrossing;
typedef struct OrbDiagram OrbDiagram;

typedef struct OrbGraph OrbGraph;

/*
 * Corresponds to EndType in gui/diagram_canvas.h, see
 * https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/gui/diagram_canvas.h#L27
 */
enum OrbDiagramEndType
{
    diagramBegin = 0,
    diagramEnd
};

typedef enum OrbDiagramEndType OrbDiagramEndType;

/* Corresponds to EdgeType in gui/diagram_canvas.h */
enum OrbDiagramEdgeType
{
    diagramSingular = 0,
    diagramDrilled
};

typedef enum OrbDiagramEdgeType OrbDiagramEdgeType;

/* Corresponds to EndData in gui/diagram_canvas.h */
struct OrbDiagramEndData
{
    OrbDiagramEdge    *edge;
    OrbDiagramEndType type;
    Boolean           singular;
    double            angle;
};

/* Corresponds to Vertex in gui/diagram_canvas.h */
struct OrbDiagramVertex
{
    int               x, y;
    int               connected_component;
    int               vertex_id;
    int               link_id;
    int               num_incident_end_data;
    OrbDiagramEndData **incident_end_data;
};

/* Corresponds to Edge in gui/diagram_canvas.h */
struct OrbDiagramEdge
{
    OrbDiagramVertex *vertex[2];
    int num_crossings;
    OrbDiagramCrossing **crossings;
    int arc_id;
    int link_id;
    int edge_id;
    OrbDiagramEdgeType edge_type;
};

/* Corresponds to Crossing in gui/diagram_canvas.h */
struct OrbDiagramCrossing
{
    int x, y;
    int crossing_id;
    int crossing_sign;
    OrbDiagramEdge *over, *under;
    double position_on_overstrand, position_on_understrand;
};

/* Corresponds to DiagramCanvas in gui/diagram_canvas.h */
struct OrbDiagram
{
    int num_arcs;
    int num_links;
    int num_vertices;
    OrbDiagramVertex **vertices;
    int num_edges;
    OrbDiagramEdge **edges;
    int num_crossings;
    OrbDiagramCrossing **crossings;
};

SNAPPEA_NAMESPACE_END_SCOPE

#endif

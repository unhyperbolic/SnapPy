#include "kernel.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

/* Corresponds to free_graph in gui/graph_complement.c */
extern void orb_free_graph(
    OrbGraph *gamma)
{
    int i;
    OrbGraphMeeting *theMeeting;

    if (gamma == NULL) {
	return;
    }

    if (gamma->meeting != NULL) {
	for (i = 0; i < gamma->num_meetings; i++) {
	    theMeeting = &gamma->meeting[i];

	    if (theMeeting->component != NULL)
		my_free(theMeeting->component);

	    if (theMeeting->label != NULL) my_free(theMeeting->label);

	    if (theMeeting->strand != NULL) my_free(theMeeting->strand);

	    if (theMeeting->tet != NULL) my_free(theMeeting->tet);

	    if (theMeeting->neighbor != NULL) my_free(theMeeting->neighbor);
	}

	my_free(gamma->meeting);
    }

    my_free(gamma);
}

SNAPPEA_NAMESPACE_END_SCOPE

#!/usr/bin/python3

"""

This run the final portion of MOde 2 assembly, from
creation of the AssemblyGraph2 to the end.

"""

import shasta

a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphReverseComplementVertex()
a.accessMarkerGraphEdges()
a.accessMarkerGraphReverseComplementEdge()
a.accessMarkerGraphConsensus()

a.createAssemblyGraph2()

# Missing code here.




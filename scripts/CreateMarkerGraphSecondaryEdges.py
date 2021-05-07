#!/usr/bin/python3

import shasta
import GetConfig

config = GetConfig.getConfig()

a = shasta.Assembler()

a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphReverseComplementVertex()
a.accessMarkerGraphEdges(True)

a.createMarkerGraphSecondaryEdges(
    minEdgeCoverage = int(config['MarkerGraph']['minEdgeCoverage']),
    minEdgeCoveragePerStrand = int(config['MarkerGraph']['minEdgeCoveragePerStrand']),
    neighborhoodSize = 3)



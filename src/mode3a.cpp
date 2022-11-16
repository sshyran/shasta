// Shasta.
#include "mode3a.hpp"
#include "MarkerGraph.hpp"
#include "mode3a-BubbleCleaner.hpp"
#include "mode3a-PackedMarkerGraph.hpp"
#include "Reads.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;
using namespace mode3a;

// Standard library.
#include "fstream.hpp"
#include <map>



Assembler::Assembler(
    uint64_t threadCount,
    uint64_t k, // Marker length
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    MarkerGraph& markerGraph) :
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    markers(markers),
    markerGraph(markerGraph)
{

    // This requires the marker length k to be even.
    SHASTA_ASSERT((k % 2) == 0);

    // This does not work with RLE.
    SHASTA_ASSERT(reads.representation == 0);

    // Clear all the superbubble fags in marker graph edges -
    // just in case we alreayd ran this before.
    for( MarkerGraph::Edge& edge: markerGraph.edges) {
        edge.isSuperBubbleEdge = 0;
    }

    // Create the initial PackedMarkerGraph.
    const string name0 = "Mode3a-PackedMarkerGraph-0";
    PackedMarkerGraph packedMarkerGraph0(name0, k, MappedMemoryOwner(*this), markers, markerGraph);
    packedMarkerGraph0.assembleSegmentSequences(name0);
    cout << "The initial PackedMarkerGraph has " <<
        packedMarkerGraph0.segments.size() << " segments, " <<
        packedMarkerGraph0.links.size() << " links, and " <<
        packedMarkerGraph0.segmentSequences.totalSize() <<
        " bases of assembled sequence." << endl;
    packedMarkerGraph0.writeGfa(name0);

    // Clean up the bubbles causes by errors.
    // This keeps one branch of each bubble.
    // The marker graph edges of the remaining branches are flagged as removed.
    BubbleCleaner cleaner(packedMarkerGraph0);
    cleaner.cleanup(markerGraph);

    // Create the cleaned up PackedMarkerGraph.
    const string name1 = "Mode3a-PackedMarkerGraph-1";
    PackedMarkerGraph packedMarkerGraph1(name1, k, MappedMemoryOwner(*this), markers, markerGraph);
    packedMarkerGraph1.assembleSegmentSequences(name1);
    cout << "After bubble cleanup, the PackedMarkerGraph has " <<
        packedMarkerGraph1.segments.size() << " segments, " <<
        packedMarkerGraph1.links.size() << " links, and " <<
        packedMarkerGraph1.segmentSequences.totalSize() <<
        " bases of assembled sequence." << endl;
    packedMarkerGraph1.writeGfa(name1);
}



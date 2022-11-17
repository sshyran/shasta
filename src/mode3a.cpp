// Shasta.
#include "mode3a.hpp"
#include "MarkerGraph.hpp"
#include "mode3a-AssemblyGraph.hpp"
#include "mode3a-AssemblyGraphSnapshot.hpp"
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
    packedMarkerGraph = make_shared<PackedMarkerGraph>(
        *this, "Mode3a-PackedMarkerGraph-Initial", k, markers, markerGraph);
    packedMarkerGraph->assembleSegmentSequences();
    cout << "The initial PackedMarkerGraph has " <<
        packedMarkerGraph->segments.size() << " segments, " <<
        packedMarkerGraph->links.size() << " links, and " <<
        packedMarkerGraph->segmentSequences.totalSize() <<
        " bases of assembled sequence." << endl;
    packedMarkerGraph->writeGfa();

    // Clean up the bubbles causes by errors.
    // This keeps one branch of each bubble.
    // The marker graph edges of the remaining branches are flagged as removed.
    BubbleCleaner cleaner(*packedMarkerGraph);
    cleaner.cleanup(markerGraph);
    packedMarkerGraph->remove();
    packedMarkerGraph = 0;

    // Create the final PackedMarkerGraph.
    packedMarkerGraph = make_shared<PackedMarkerGraph>(
        *this, "Mode3a-PackedMarkerGraph", k, markers, markerGraph);
    packedMarkerGraph->assembleSegmentSequences();
    cout << "After bubble cleanup, the PackedMarkerGraph has " <<
        packedMarkerGraph->segments.size() << " segments, " <<
        packedMarkerGraph->links.size() << " links, and " <<
        packedMarkerGraph->segmentSequences.totalSize() <<
        " bases of assembled sequence." << endl;
    packedMarkerGraph->writeGfa();

    // Create the AssemblyGraph.
    packedMarkerGraph->createMarkerGraphEdgeTable(threadCount);
    packedMarkerGraph->computeJourneys(threadCount);
    packedMarkerGraph->writeJourneys();
    AssemblyGraph assemblyGraph(*packedMarkerGraph);
    packedMarkerGraph->journeys.shrink_to_fit();
    packedMarkerGraph->markerGraphEdgeTable.remove();

    // Write it out.
    cout << "The initial AssemblyGraph has " <<
        num_vertices(assemblyGraph) << " segments and " <<
        num_edges(assemblyGraph) << " links." << endl;
    assemblyGraph.write("Mode3a-AssemblyGraph-Initial");
    AssemblyGraphSnapshot snapshot(assemblyGraph, "Mode3a-AssemblyGraph-Initial", *this);
}



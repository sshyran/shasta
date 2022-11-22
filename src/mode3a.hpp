#ifndef SHASTA_MODE3A_HPP
#define SHASTA_MODE3A_HPP

/*******************************************************************************

Namespace mode3a contains newer code for Mode 3 assembly.
The code in namespace mode3 will eventually become obsolete.

The top level class is Assembler.

Class PackedMarkerGraph is used to
construct a "packed" representations of the marker graph.
Each segment corresponds to a
linear sequence of marker graph edges
without any intervening incoming/outgoing edges in the marker graph.
Connectivity is generated based on marker graph connectivity.
That is, a link between segments s0 and s1 is created
if the last marker graph vertex of s0 is the same as the first
marker graph vertex of s1.
Because connectivity does not "follow the reads", the
PackedMarkerGraph is subject to fragmentation.

We use two instances of PackedMarkerGraph:
- An initial, temporary PackedMarkerGraph created from the unmodified MarkerGraph.
- A final, permanent PackedMarkerGraph created from the MarkerGraph after some of
  its edges were flagged as due to errors and marked as removed.

Class BubbleCleaner takes as input the initial, temporary PackedMarkerGraph
and flags marker graph edges of branch bubbles that are
likely to be caused by errors.

The rest of the work (mostly detangling) is done by class AssemblyGraph.
Here, each segment corresponds to a segment of the (final)
PackedMarkerGraph. There can be multiple AssemblyGraph segments
corresponding to a single PackedMarkerGraph segment.
For each of the copies we store which oriented reads are
believed to "belong" to that copy. Copies of segments
are created during the detangling process.

The AssemblyGraph is stored as a boost graph and is not persistent.
This makes it easy to manipulate, which is important for detangling.
Each segment is a vertex and each link is an edge
(contrary to what happens in the PackedMarkerGraph).

Connectivity in the AssemblyGraph is generated by "following the reads":
a link s0->s1 between segments s0 and s1 is created if an oriented
read visits segment s1 immediately after visiting segment s0.
The coverage of a link is the number of oriented reads that
have such a transition. Depending on what we are doing,
we may only be interested in links with a minimum amount of coverage.
Because of this definition, the sequence of segments visited
by an oriented read is a path in the Assembly graph.

Class AssemblyGraphSnapshot is a persistent snapshot of
the Assembly graph. It uses MemoryMapped classes so it can be accessed
later for use in the http server or the Python API.
We store at least an initial snapshot immediately after creation
and a final snapshot after detangling.

*******************************************************************************/

// Shasta.
#include "Base.hpp"
#include "MappedMemoryOwner.hpp"
#include "MemoryMappedVectorOfVectors.hpp"

// Standard library.
#include <memory>
#include "string.hpp"


namespace shasta {

    // Forward declarations of the classesd defined in this file.
    namespace mode3a {
        class Assembler;
        class PackedMarkerGraph;
        class AssemblyGraph;
    }

    // Forward declarations of Shasta classes defined elsewhere.
    class Reads;
    class CompressedMarker;
    class MarkerGraph;
}



class shasta::mode3a::Assembler :
    public MappedMemoryOwner {
public:

    Assembler(
        uint64_t threadCount,
        uint64_t k, // Marker length
        const MappedMemoryOwner&,
        const Reads&,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        MarkerGraph&);

    uint64_t k;
    const Reads& reads;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;

    // The MarkerGraph is not const because the Assembler uses the BubbleCleaner
    // to flag marker graph edges.
    MarkerGraph& markerGraph;

    std::shared_ptr<PackedMarkerGraph> packedMarkerGraph;
};



#endif

#ifndef SHASTA_MODE3A_HPP
#define SHASTA_MODE3A_HPP

/*******************************************************************************

Namespace mode3a contains newer code for Mode 3 assembly.
The code in namespace mode3 will eventually become obsolete.

The top level class is Assembler.

Class PackedMarkerGraph is used to
construct two "packed" representations of the marker graph.
Here, each segment corresponds to a path in the marker graph
Connectivity is generated based on marker graph connectivity.
That is, a link between segments s0 and s1 is created
if the last marker graph vertex of s0 is the same as the first
marker graph vertex of s1.
Because connectivity does not "follow the reads", the
mode3a::PackedMarkerGraph is subject to fragmentation.

We use two instances of PackedMarkerGraph(s).
In both versions, each segment corresponds to a
linear sequence of marker graph edges
without any intervening incoming/outgoing edges in the marker graph.
- In the first instance, all marker graph edges are used to define
the segments.
- In the second instance, only marker graph edges that were not
flagged by the BubbleCleaner are used.

Class BubbleCleaner takes as input the initial PackedMarkerGraph
and flags marker graph edges of branch bubbles that are
likely to be caused by errors.

Class Detangler is the workhorse class for detangling.
It works iteratively, and so needs to be a dynamic data structure
that can be updated easily, and for that reason it is represented
as a boost::graph::adjacency_list.
Here, each segment represents a sequence (not necessarily a path)
of mode3a::PackedMarkerGraph segments.
Connectivity is created and maintained by "following the reads".
This is constructed using as input the PackedMarkerGraph
after bubble removal.

*******************************************************************************/

// Shasta.
#include "Base.hpp"
#include "MappedMemoryOwner.hpp"
#include "MemoryMappedVectorOfVectors.hpp"

// Standard library.
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
};



#endif


#ifndef SHASTA_MODE3A_HPP
#define SHASTA_MODE3A_HPP

/*******************************************************************************

Namespace mode3a contains newer code for Mode 3 assembly.
The code in namespace mode3 will eventually become obsolete.

The top level class is mode3a::Assembler.

Class mode3a::PackedMarkerGraph is used to initially
construct a "packed" representation of the marker graph.
Here, each segment corresponds to a path in the marker graph
Connectivity is generated based on marker graph connectivity.
That is, a link between segments s0 and s1 is created
if the last marker graph vertex of s0 is the same as the first
marker graph vertex of s1.
Because connectivity does not "follow the reads", the
mode3a::PackedMarkerGraph is subject to fragmentation.

We use two mode3a::PackedMarkerGraph(s):
- An initial one in which each segment corresponds to a linear sequence of marker graph edges
  without any intervening incoming/outgoing edges.

- A cleaned up version, in which
  the bubbles dues to errors are removed.

Class mode3a::BubbleCleaner takes as input the initial mode3a::PackedMarkerGraph
and creates the cleaned up one.

Class Mode3::AssemblyGraph is the workhorse class for detangling.
It works iteratively, and so needs to be a dynamic data structure
that can be updated easily, and for that reason it is represented
as a boost::graph::adjacency_list.
Here, each segment represents a sequence (not necessarily a path)
of mode3a::PackedMarkerGraph segments (the cleaned up version).
Connectivity is created and maintained by "following the reads".

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
        class BubbleCleaner;
        class TidyMarkerGraph;
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
        const MarkerGraph&);

    uint64_t k;
    const Reads& reads;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;
    const MarkerGraph& markerGraph;
};



// In PackedMarkerGraph, each segment corresponds exactly to a
// path in the marker graph.
class shasta::mode3a::PackedMarkerGraph :
    public MappedMemoryOwner {
public:
    PackedMarkerGraph(
        const string& name,
        uint64_t k,
        const MappedMemoryOwner&,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph&);

    uint64_t k;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;
    const MarkerGraph& markerGraph;

    // For each segment, we store its path in the marker graph.
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> segments;

    // This is used for intial creation from the marker graph.
    void createSegmentsFromMarkerGraph(const string& name);

    // Get the first or last marker graph vertex of a segment.
    uint64_t getFirstSegmentVertex(uint64_t segmentId) const;
    uint64_t getLastSegmentVertex (uint64_t segmentId) const;

    // Assembled sequence of each segment.
    MemoryMapped::VectorOfVectors<Base, uint64_t> segmentSequences;
    void assembleSegmentSequences(const string& name);

    // A link between segments s0 and s1 is created if
    // the last marker graph vertex of s0 coincides with the
    // first marker graph vertex of s1.
    class Link {
    public:
        uint64_t segmentId0;
        uint64_t segmentId1;
        Link() {}
        Link(
            uint64_t segmentId0,
            uint64_t segmentId1) :
            segmentId0(segmentId0),
            segmentId1(segmentId1) {}

    };
    MemoryMapped::Vector<Link> links;

    // This creates the links once the segments are available.
    void createLinks(const string& name);

    // Given the links, create the connectivity.
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> linksBySource;
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> linksByTarget;
    void createConnectivity(const string& name);

    void writeGfa(const string& name) const;
};

#endif


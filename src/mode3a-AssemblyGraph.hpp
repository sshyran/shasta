#ifndef SHASTA_MODE3A_ASSEMBLY_GRAPH_HPP
#define SHASTA_MODE3A_ASSEMBLY_GRAPH_HPP

// See comments in mode3a.hpp.

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

namespace shasta {
    namespace mode3a {

        class AssemblyGraph;
        class AssemblyGraphVertex;
        class AssemblyGraphEdge;
        using AssemblyGraphBaseClass = boost::adjacency_list<
            boost::listS, boost::listS, boost::bidirectionalS,
            AssemblyGraphVertex, AssemblyGraphEdge>;

        class PackedMarkerGraph;
    }
}



// Each vertex represents a segment.
class shasta::mode3a::AssemblyGraphVertex {
public:

    // The PackedMarkerGraph segment that this vertex
    // corresponds to.
    uint64_t segmentId;
};



// Each edge represents a link.
class shasta::mode3a::AssemblyGraphEdge {
public:
};



class shasta::mode3a::AssemblyGraph : public AssemblyGraphBaseClass {
public:
    AssemblyGraph(const PackedMarkerGraph&);
private:

    // Each segment in the AssemblyGraph corresponds to a segment in
    // this PackedMarkerGraph.
    const PackedMarkerGraph& packedMarkerGraph;
};

#endif


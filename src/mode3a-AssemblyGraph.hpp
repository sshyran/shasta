#ifndef SHASTA_MODE3A_ASSEMBLY_GRAPH_HPP
#define SHASTA_MODE3A_ASSEMBLY_GRAPH_HPP

// See comments in mode3a.hpp.

// Shasta.
#include "ReadId.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <list>
#include "tuple.hpp"
#include "vector.hpp"

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

    AssemblyGraphVertex(uint64_t segmentId) :
        segmentId(segmentId) {}

    // The path entries that go through this vertex.
    class PathEntry {
    public:
        OrientedReadId orientedReadId;
        std::list<AssemblyGraphBaseClass::vertex_descriptor>::iterator it;

        PathEntry(
            OrientedReadId orientedReadId,
            std::list<AssemblyGraphBaseClass::vertex_descriptor>::iterator it) :
            orientedReadId(orientedReadId), it(it) {}
    };
    vector<PathEntry> pathEntries;
};



// Each edge represents a link.
class shasta::mode3a::AssemblyGraphEdge {
public:

    // A transition of an oriented read from a segment to the next.
    class Transition {
    public:
        std::list<AssemblyGraphBaseClass::vertex_descriptor>::iterator it0;
        std::list<AssemblyGraphBaseClass::vertex_descriptor>::iterator it1;
        OrientedReadId orientedReadId;

        bool operator<(const Transition& that) const
        {
            return tie(*it0, *it1, orientedReadId) < tie(*that.it0, *that.it1, that.orientedReadId);
        }
    };

    vector<Transition> transitions;

    uint64_t coverage() const
    {
        return transitions.size();
    }
};



class shasta::mode3a::AssemblyGraph : public AssemblyGraphBaseClass {
public:
    AssemblyGraph(const PackedMarkerGraph&);
    void write(const string& name) const;
private:

    // Each segment in the AssemblyGraph corresponds to a segment in
    // this PackedMarkerGraph.
    const PackedMarkerGraph& packedMarkerGraph;

    // The sequence of segments visited by each oriented read
    // is a path. Initially, we construct it from the corresponding journey
    // in the PackedMarkerGraph.
    // Indexed by OrientedReadId.getValue().
    vector< std::list<vertex_descriptor> > paths;

    void createSegmentsAndPaths();
    void createLinks();

    void writeGfa(const string& name, uint64_t minLinkCoverage) const;
    void writeLinkCoverageHistogram(const string& name) const;
};

#endif


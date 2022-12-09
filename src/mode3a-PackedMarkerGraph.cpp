// Shasta
#include "mode3a-PackedMarkerGraph.hpp"
#include "assembleMarkerGraphPath.hpp"
#include "AssembledSegment.hpp"
#include "invalid.hpp"
#include "MarkerGraph.hpp"
#include "orderPairs.hpp"
using namespace shasta;
using namespace mode3a;

// Standard library.
#include "fstream.hpp"
#include <map>

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<PackedMarkerGraph>;



PackedMarkerGraph::PackedMarkerGraph(
    const MappedMemoryOwner& mappedMemoryOwner,
    const string& name,
    uint64_t k,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph,
    bool accessExisting) :
    MappedMemoryOwner(mappedMemoryOwner),
    MultithreadedObject<PackedMarkerGraph>(*this),
    name(name),
    k(k),
    markers(markers),
    markerGraph(markerGraph)
{
    if(accessExisting) {
        accessSegments();
        accessLinks();
        accessSegmentSequences();
        accessJourneys();
    } else {
        // Initial creation from the marker graph.
        createSegmentsFromMarkerGraph();
        createLinks();
        assembleSegmentSequences();
    }
}



void PackedMarkerGraph::createSegmentsFromMarkerGraph()
{
    createNew(segments, name + "-segments");

    const bool debug = false;

    const uint64_t edgeCount = markerGraph.edges.size();
    vector<bool> wasFound(edgeCount, false);

    using MarkerGraphPath = vector<uint64_t>;
    MarkerGraphPath nextEdges;
    MarkerGraphPath previousEdges;
    MarkerGraphPath path;
    MarkerGraphPath reverseComplementedPath;

    // Main loop over all edges of the marker graph.
    // At each iteration we find a new linear path of edges.
    for(uint64_t startEdgeId=0; startEdgeId<edgeCount; startEdgeId++) {

        if(markerGraph.edges[startEdgeId].wasRemoved()) {
            continue;
        }

        // If we already found this edge, skip it.
        // It is part of a path we already found.
        if(wasFound[startEdgeId]) {
            continue;
        }

        if(debug) {
            cout << "Starting a new path at edge " << startEdgeId << endl;
        }

        // Follow the path forward.
        nextEdges.clear();
        MarkerGraph::EdgeId edgeId = startEdgeId;
        bool isCircular = false;
        while(true) {
            const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
            SHASTA_ASSERT(not edge.wasRemoved());
            const MarkerGraph::VertexId v1 = edge.target;
            if(markerGraph.outDegree(v1) != 1) {
                break;
            }
            if(markerGraph.inDegree(v1) != 1) {
                break;
            }
            edgeId = markerGraph.getFirstNonRemovedOutEdge(v1);
            if(edgeId == startEdgeId) {
                isCircular = true;
                break;
            }
            nextEdges.push_back(edgeId);
            SHASTA_ASSERT(not wasFound[edgeId]);
            if(debug) {
                cout << "Moving forward: added " << edgeId << endl;
            }
        }

        // Follow the path backward.
        previousEdges.clear();
        if(!isCircular) {
            edgeId = startEdgeId;
            while(true) {
                const MarkerGraph::Edge edge = markerGraph.edges[edgeId];
                SHASTA_ASSERT(not edge.wasRemoved());
                const MarkerGraph::VertexId v0 = edge.source;
                if(markerGraph.outDegree(v0) != 1) {
                    break;
                }
                if(markerGraph.inDegree(v0) != 1) {
                    break;
                }
                edgeId = markerGraph.getFirstNonRemovedInEdge(v0);
                previousEdges.push_back(edgeId);
                SHASTA_ASSERT(not wasFound[edgeId]);
                if(debug) {
                    cout << "Moving backward: added " << edgeId << endl;
                }
            }
        }

        // Gather the path.
        path.clear();
        copy(previousEdges.rbegin(), previousEdges.rend(), back_inserter(path));
        path.push_back(startEdgeId);
        copy(nextEdges.begin(), nextEdges.end(), back_inserter(path));

        // Mark all the edges in the path as found.
        for(const MarkerGraph::EdgeId edgeId: path) {
            if(wasFound[edgeId]) {
                cout << "Assertion failed at " << edgeId << endl;
                SHASTA_ASSERT(0);
            }
            wasFound[edgeId] = true;
        }

        // Store this path as a new segment.
        segments.appendVector(path);
    }



    // Check that all non-removed edges of the marker graph were found.
    for(uint64_t edgeId=0; edgeId<markerGraph.edges.size(); edgeId++) {
        const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
        SHASTA_ASSERT(edge.wasRemoved() or wasFound[edgeId]);
    }
}



void PackedMarkerGraph::accessSegments()
{
    accessExistingReadOnly(segments, name + "-segments");
}



void PackedMarkerGraph::accessLinks()
{
    accessExistingReadOnly(links, name + "-links");
}



// Get the first or last marker graph vertex of a PackedMarkerGraph segment.
uint64_t PackedMarkerGraph::getFirstSegmentVertex(uint64_t segmentId) const
{
    const auto segment = segments[segmentId];
    SHASTA_ASSERT(not segment.empty());
    const uint64_t firstEdgeId = segment.front();
    const MarkerGraph::Edge& edge = markerGraph.edges[firstEdgeId];
    return edge.source;
}
uint64_t PackedMarkerGraph::getLastSegmentVertex (uint64_t segmentId) const
{
    const auto segment = segments[segmentId];
    SHASTA_ASSERT(not segment.empty());
    const uint64_t lastEdgeId = segment.back();
    const MarkerGraph::Edge& edge = markerGraph.edges[lastEdgeId];
    return edge.target;
}




void PackedMarkerGraph::createLinks()
{
    // Index the segments by their first vertex.
    // This can be made faster.
    std::map<uint64_t, vector<uint64_t> > segmentMap;    // Key=first vertex, Value=segments
    for(uint64_t segmentId=0; segmentId!=segments.size(); segmentId++) {
        segmentMap[getFirstSegmentVertex(segmentId)].push_back(segmentId);
    }


    // Create the links.
    createNew(links, name + "-links");
    for(uint64_t segmentId0=0; segmentId0!=segments.size(); segmentId0++) {
        const auto& v = segmentMap[getLastSegmentVertex(segmentId0)];
        for(const uint64_t segmentId1: v) {
            links.push_back(Link(segmentId0, segmentId1));
        }
    }
}



void PackedMarkerGraph::assembleSegmentSequences()
{
    createNew(segmentSequences, name + "-sequences");

    for(uint64_t segmentId=0; segmentId!=segments.size(); segmentId++) {

        AssembledSegment assembledSegment;
        assembleMarkerGraphPath(
            0,
            k,
            markers,
            markerGraph,
            segments[segmentId],
            false,
            assembledSegment);

        segmentSequences.appendVector(
            assembledSegment.rawSequence.begin() + k/2,
            assembledSegment.rawSequence.end()   - k/2);
    }
}



void PackedMarkerGraph::accessSegmentSequences()
{
    accessExistingReadOnly(segmentSequences, name + "-sequences");
}



void PackedMarkerGraph::writeGfa() const
{
    ofstream gfa(name + ".gfa");
    ofstream csv(name + ".csv");

    // Write the headers.
    gfa << "H\tVN:Z:1.0\n";
    csv << "Segment,Sequence Length,Path Length,First marker graph vertex,Last marker graph vertex\n";

    // Write the segments.
    for(uint64_t segmentId=0; segmentId<segmentSequences.size(); segmentId++) {

        const auto sequence = segmentSequences[segmentId];
        gfa <<"S\t" << segmentId << "\t";
        copy(sequence.begin()+k/2, sequence.end()-k/2, ostream_iterator<Base>(gfa));
        gfa << "\n";

        const auto path = segments[segmentId];
        csv << segmentId << ",";
        csv << sequence.size() << ",";
        csv << path.size() << ",";
        csv << getFirstSegmentVertex(segmentId) << ",";
        csv << getLastSegmentVertex(segmentId) << "\n";
    }

    // Write the links.
    for(const Link& link: links) {
        gfa << "L\t" <<
            link.segmentId0 << "\t+\t" <<
            link.segmentId1 << "\t+\t0M\n";
    }


}


void PackedMarkerGraph::writeSegments()
{
    ofstream csv(name + "-Segments.csv");
    csv << "Segment,MarkerGraph edges\n";
    for(uint64_t segmentId=0; segmentId<segments.size(); segmentId++) {
        const auto segment = segments[segmentId];
        for(const uint64_t markerGraphEdgeId: segment) {
            csv << segmentId << "," << markerGraphEdgeId << "\n";
        }
    }

}


void PackedMarkerGraph::remove()
{
    if(segments.isOpen()) {
        segments.remove();
    }
    if(segmentSequences.isOpen()) {
        segmentSequences.remove();
    }
    if(links.isOpen) {
        links.remove();
    }
    if(linksBySource.isOpen()) {
        linksBySource.remove();
    }
    if(linksByTarget.isOpen()) {
        linksByTarget.remove();
    }
}



void PackedMarkerGraph::createMarkerGraphEdgeTable(uint64_t threadCount)
{
    // Shorthand for the marker graph edge table.
    auto& table = computeJourneysData.markerGraphEdgeTable;

    // Initialize the table.
    createNew(table, name + "-MarkerGraphEdgeTable");
    table.resize(markerGraph.edges.size());
    fill(table.begin(), table.end(), invalid<uint64_t>);

    // Fill it in.
    const uint64_t batchSize = 100;
    setupLoadBalancing(segments.size(), batchSize);
    runThreads(&PackedMarkerGraph::createMarkerGraphEdgeTableThreadFunction, threadCount);
}



void PackedMarkerGraph::createMarkerGraphEdgeTableThreadFunction(uint64_t threadId)
{
    auto& table = computeJourneysData.markerGraphEdgeTable;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all segments assigned to this batch.
        for(uint64_t segmentId=begin; segmentId!=end; ++segmentId) {
            const auto segment = segments[segmentId];

            // Loop over the marker graph path of this segment.
            for(const uint64_t markerGraphEdgeId: segment) {
                table[markerGraphEdgeId] = segmentId;
            }
        }
    }
}



#if 0
// Old code, buggy.
void PackedMarkerGraph::computeJourneys(uint64_t threadCount)
{
    journeys.clear();   // Just in case.
    const uint64_t orientedReadCount = markers.size();
    journeys.resize(orientedReadCount);

    const uint64_t batchSize = 100;
    setupLoadBalancing(orientedReadCount, batchSize);
    runThreads(&PackedMarkerGraph::computeJourneysThreadFunction, threadCount);
}



void PackedMarkerGraph::computeJourneysThreadFunction(uint64_t threadId)
{
    // Work vectors for computeJourney.
    vector<uint64_t> markerGraphVertices;
    vector<uint64_t> markerGraphEdges;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(uint64_t i=begin; i!=end; i++) {
            computeJourney(
                OrientedReadId::fromValue(ReadId(i)),
                markerGraphVertices,
                markerGraphEdges);
        }
    }

}



void PackedMarkerGraph::computeJourney(
    OrientedReadId orientedReadId,
    vector<uint64_t>& markerGraphVertices,
    vector<uint64_t>& markerGraphEdges)
{
    // The MarkerId interval for this oriented read.
    const uint64_t markerIdBegin = markers.begin(orientedReadId.getValue()) - markers.begin();
    const uint64_t markerIdEnd = markers.end(orientedReadId.getValue()) - markers.begin();

    // Get the sequence of marker graph vertices encountered by this read.
    markerGraphVertices.clear();
    for(uint64_t markerId=markerIdBegin; markerId!=markerIdEnd; markerId++) {
        const MarkerGraph::CompressedVertexId compressedVertexId = markerGraph.vertexTable[markerId];
        if(compressedVertexId != MarkerGraph::invalidCompressedVertexId) {
            markerGraphVertices.push_back(compressedVertexId);
        }
    }

    // Get the sequence of marker graph edges encountered by this read.
    markerGraphEdges.clear();
    for(uint64_t i=1; i<markerGraphVertices.size(); i++) {
        const uint64_t vertexId0 = markerGraphVertices[i-1];
        const uint64_t vertexId1 = markerGraphVertices[i];
        const MarkerGraph::Edge* edge = markerGraph.findEdge(Uint40(vertexId0), Uint40(vertexId1));
        if(edge and not edge->wasRemoved()) {
            const uint64_t edgeId = edge - &markerGraph.edges[0];
            markerGraphEdges.push_back(edgeId);
        }
    }


    // Now we can use the markerGraphEdgeTable to compute the journey.
    vector<uint64_t>& journey = journeys[orientedReadId.getValue()];
    journey.clear();
    for(const uint64_t markerGraphEdgeId: markerGraphEdges) {
        const uint64_t segmentId = markerGraphEdgeTable[markerGraphEdgeId];
        if(segmentId != invalid<uint64_t>) {
            if(journey.empty() or segmentId != journey.back()) {
                journey.push_back(segmentId);
            }
        }
    }
}
#endif



void PackedMarkerGraph::computeJourneys(uint64_t threadCount)
{
    auto& data = computeJourneysData;

    // The marker graph edge table gives the segmentId that each
    // marker graph edge is on, or invalid<uint64_t> if none.
    createMarkerGraphEdgeTable(threadCount);

    // Compute the journey pairs for each oriented read id.
    const uint64_t orientedReadCount = markers.size();
    createNew(data.journeyPairs, name + "-JourneyPairs");
    data.journeyPairs.beginPass1(orientedReadCount);
    const uint64_t batchSize = 100;
    setupLoadBalancing(segments.size(), batchSize);
    runThreads(&PackedMarkerGraph::computeJourneysPass1ThreadFunction, threadCount);
    data.journeyPairs.beginPass2();
    setupLoadBalancing(segments.size(), batchSize);
    runThreads(&PackedMarkerGraph::computeJourneysPass2ThreadFunction, threadCount);
    data.journeyPairs.endPass2(true, true);

    // We no longer need the marker graph edge table.
    data.markerGraphEdgeTable.remove();

    // For each oriented read, sort the journey pairs and use them to compute the journeys.
    // The journeys are temporarily stored in data.journeys.
    data.journeys.resize(orientedReadCount);
    setupLoadBalancing(orientedReadCount, batchSize);
    runThreads(&PackedMarkerGraph::computeJourneysPass3ThreadFunction, threadCount);

    // We no longer need the journey pairs.
    data.journeyPairs.remove();

    // Copy the journeys to their permanent location in mapped memory.
    createNew(journeys, name + "-Journeys");
    for(const auto& journey: data.journeys) {
        journeys.appendVector(journey);
    }

    // We no longer need the temporary copy of the journeys.
    data.journeys.clear();
    data.journeys.shrink_to_fit();
}



void PackedMarkerGraph::accessJourneys()
{
    accessExistingReadOnly(journeys, name + "-Journeys");
}



void PackedMarkerGraph::computeJourneysPass1ThreadFunction(uint64_t threadId)
{
    computeJourneysPass12ThreadFunction(1);
}



void PackedMarkerGraph::computeJourneysPass2ThreadFunction(uint64_t threadId)
{
    computeJourneysPass12ThreadFunction(2);
}



// This is used to compute data.journeyPairs for each oriented read.
void PackedMarkerGraph::computeJourneysPass12ThreadFunction(uint64_t pass)
{
    auto& data = computeJourneysData;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over segments assigned to this batch.
        for(uint64_t segmentId=begin; segmentId!=end; segmentId++) {
            const auto segment = segments[segmentId];

            // Loop over marker graph edges of this segment.
            for(uint64_t edgeId: segment) {

                // Loop over marker intervals of this edge.
                for(const MarkerInterval& markerInterval: markerGraph.edgeMarkerIntervals[edgeId]) {

                    if(pass == 1) {
                        data.journeyPairs.incrementCountMultithreaded(
                            markerInterval.orientedReadId.getValue());
                    } else {
                        data.journeyPairs.storeMultithreaded(
                            markerInterval.orientedReadId.getValue(),
                            make_pair(markerInterval.ordinals[0], segmentId));
                    }
                }
            }
        }
    }
}



// In pass 3 we use the journey pairs to compute journeys for each oriented read
// and we store them temporarily in data.journeys.
void PackedMarkerGraph::computeJourneysPass3ThreadFunction(uint64_t threadId)
{
    auto& data = computeJourneysData;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over oriented reads in this batch.
        for(uint64_t i=begin; i!=end; ++i) {

            // Sort by marker ordinal the journey pairs for each read.
            auto journeyPairs = data.journeyPairs[i];
            sort(journeyPairs.begin(), journeyPairs.end(), OrderPairsByFirstOnly<uint32_t, uint64_t>());

            vector<uint64_t>& journey = data.journeys[i];
            for(const auto& journeyPair: journeyPairs) {
                const uint64_t segmentId = journeyPair.second;
                if(journey.empty() or segmentId != journey.back()) {
                    journey.push_back(segmentId);
                }
            }
        }
    }
}



void PackedMarkerGraph::writeJourneys() const
{
    ofstream csv(name + "-Journeys.csv");
    for(uint64_t i=0; i<journeys.size(); i++) {
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(i));
        csv << orientedReadId << ",";
        const auto journey = journeys[i];
        for(const uint64_t segmentId: journey) {
            csv << segmentId << ",";
        }
        csv << "\n";
    }
}



// Compute average marker graph edge coverage for a segment.
double PackedMarkerGraph::averageMarkerGraphEdgeCoverage(uint64_t segmentId) const
{
    const auto path = segments[segmentId];
    uint64_t sum = 0;
    for(const uint64_t edgeId: path) {
        sum += markerGraph.edgeMarkerIntervals.size(edgeId);
    }
    return double(sum) / double(path.size());
}

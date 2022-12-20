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



void PackedMarkerGraph::computeJourneys(uint64_t threadCount)
{
    auto& data = computeJourneysData;

    // Compute the marker graph journey of each oriented read.
    const uint64_t orientedReadCount = markers.size();
    createNew(data.markerGraphJourneys, name + "-MarkerGraphJourneys");
    createNew(data.journeyPairs, name + "-JourneyPairs");
    data.journeyPairs.beginPass1(orientedReadCount);
    data.markerGraphJourneys.beginPass1(orientedReadCount);
    const uint64_t batchSize = 100;
    setupLoadBalancing(segments.size(), batchSize);
    runThreads(&PackedMarkerGraph::computeJourneysPass1ThreadFunction, threadCount);
    data.journeyPairs.beginPass2();
    data.markerGraphJourneys.beginPass2();
    setupLoadBalancing(segments.size(), batchSize);
    runThreads(&PackedMarkerGraph::computeJourneysPass2ThreadFunction, threadCount);
    data.journeyPairs.endPass2(true, true);
    data.markerGraphJourneys.endPass2(true, true);

    // For each oriented read, sort the journey pairs and use them to compute the journeys.
    // The journeys are temporarily stored in data.journeys.
    data.journeys.resize(orientedReadCount);
    data.detailedJourneys.resize(orientedReadCount);
    setupLoadBalancing(orientedReadCount, batchSize);
    runThreads(&PackedMarkerGraph::computeJourneysPass3ThreadFunction, threadCount);

    // We no longer need the marker graph journeys.
    data.journeyPairs.remove();
    data.markerGraphJourneys.remove();

    // Copy the journeys to their permanent location in mapped memory.
    createNew(journeys, name + "-Journeys");
    for(const auto& journey: data.journeys) {
        journeys.appendVector(journey);
    }
    createNew(detailedJourneys, name + "-DetailedJourneys");
    for(const auto& detailedJourney: data.detailedJourneys) {
        detailedJourneys.appendVector(detailedJourney);
    }

    // Check the the journeys and the detailed journeys agree.
    SHASTA_ASSERT(journeys.size() == detailedJourneys.size());
    for(uint64_t i=0; i<journeys.size(); i++) {
        const auto journey = journeys[i];
        const auto detailedJourney = detailedJourneys[i];
        SHASTA_ASSERT(journey.size() == detailedJourney.size());
        for(uint64_t j=0; j<journey.size(); j++) {
            SHASTA_ASSERT(journey[j] == detailedJourney[j].segmentId);
        }
    }

    // We no longer need the temporary copy of the journeys.
    data.journeys.clear();
    data.journeys.shrink_to_fit();
}



void PackedMarkerGraph::accessJourneys()
{
    accessExistingReadOnly(journeys, name + "-Journeys");
    accessExistingReadOnly(detailedJourneys, name + "-DetailedJourneys");
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

            // Loop over marker graph edges of this segment to update the journey pairs.
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

            // Loop over marker graph edges of this segment to update the marker graph journeys.
            for(uint64_t position=0; position<segment.size(); position++) {
                const uint64_t edgeId = segment[position];

                // Loop over marker intervals of this edge.
                for(const MarkerInterval& markerInterval: markerGraph.edgeMarkerIntervals[edgeId]) {

                    if(pass == 1) {
                        data.markerGraphJourneys.incrementCountMultithreaded(
                            markerInterval.orientedReadId.getValue());
                    } else {
                        data.markerGraphJourneys.storeMultithreaded(
                            markerInterval.orientedReadId.getValue(),
                            {segmentId, position, edgeId, markerInterval.ordinals[0], markerInterval.ordinals[1]});
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

            // Sort by marker ordinal the journey pairs for each oriented read.
            auto journeyPairs = data.journeyPairs[i];
            sort(journeyPairs.begin(), journeyPairs.end(), OrderPairsByFirstOnly<uint32_t, uint64_t>());

            // Sort by marker ordinal the marker graph journeys for each read.
            auto markerGraphJourney = data.markerGraphJourneys[i];
            sort(markerGraphJourney.begin(), markerGraphJourney.end());

            // Check that the journey pairs and the marker graph journeys agree.
            SHASTA_ASSERT(journeyPairs.size() == markerGraphJourney.size());
            for(uint64_t i=0; i<journeyPairs.size(); i++) {
                SHASTA_ASSERT(journeyPairs[i].first == markerGraphJourney[i].ordinal0);
                SHASTA_ASSERT(journeyPairs[i].second == markerGraphJourney[i].segmentId);
            }

            // Use the journeyPairs to compute the journey for this oriented read.
            vector<uint64_t>& journey = data.journeys[i];
            for(const auto& journeyPair: journeyPairs) {
                const uint64_t segmentId = journeyPair.second;
                if(journey.empty() or segmentId != journey.back()) {
                    journey.push_back(segmentId);
                }
            }



            // Use the marker graph journey to compute the detailed journey for this read.
            vector<JourneyStep>& detailedJourney = data.detailedJourneys[i];
            JourneyStep journeyStep;
            for(uint64_t j=0; j<markerGraphJourney.size(); j++) {
                const auto& markerGraphJourneyStep =  markerGraphJourney[j];

                // If this is the first step in a segment, initialize the JourneyStep.
                if( j==0 or
                    markerGraphJourneyStep.segmentId !=journeyStep.segmentId) {
                    journeyStep.segmentId = markerGraphJourneyStep.segmentId;
                    journeyStep.positions[0] = markerGraphJourneyStep.positionInSegment;
                    journeyStep.ordinals[0] = markerGraphJourneyStep.ordinal0;
                }

                // If this is the last step in a segment, finalize the JourneyStep and store it.
                if(j==markerGraphJourney.size()-1 or
                    markerGraphJourney[j+1].segmentId !=journeyStep.segmentId) {
                    journeyStep.positions[1] = markerGraphJourneyStep.positionInSegment + 1;
                    journeyStep.ordinals[1] = markerGraphJourneyStep.ordinal1;
                    detailedJourney.push_back(journeyStep);
                }
            }

            // Check that the journey and the detailedJourney agree.
            SHASTA_ASSERT(journey.size() == detailedJourney.size());
            for(uint64_t j=0; j<journey.size(); j++) {
                SHASTA_ASSERT(journey[j] == detailedJourney[j].segmentId);
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

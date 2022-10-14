#include "mode3-Detangler.hpp"
#include "mode3.hpp"
using namespace shasta;
using namespace mode3;



Detangler::Detangler(const AssemblyGraph& assemblyGraph)
{
    createJourneys(assemblyGraph);
    createInitialClusters();
    cout << "The initial Detangler has " << clusters.size() << " clusters." << endl;
}




// To create the journeys, simply extract the segmentIds from the assemblyGraphJourneys.
void Detangler::createJourneys(const AssemblyGraph& assemblyGraph)
{
    const uint64_t journeyCount = assemblyGraph.assemblyGraphJourneys.size();

    journeys.clear();
    journeys.resize(journeyCount);
    for(uint64_t i=0; i<journeyCount; i++) {
        const span<const AssemblyGraphJourneyEntry> assemblyGraphJourney = assemblyGraph.assemblyGraphJourneys[i];
        Journey& journey = journeys[i];

        for(const AssemblyGraphJourneyEntry& assemblyGraphJourneyEntry: assemblyGraphJourney) {
            journey.push_back(assemblyGraphJourneyEntry.segmentId);
        }
    }
}



void Detangler::createInitialClusters()
{
    const ReadId readCount = ReadId(journeys.size() / 2);
    clusterTable.resize(journeys.size());

    // Loop over all oriented reads.
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);

            // Get the Journey for this oriented read.
            const Journey& journey = journeys[orientedReadId.getValue()];

            // Loop over entries in this Journey.
            JourneyEntryId journeyEntryId;
            journeyEntryId.orientedReadId = orientedReadId;
            clusterTable[orientedReadId.getValue()].resize(journey.size());
            for(uint64_t position=0; position<journey.size(); position++) {
                journeyEntryId.position = position;
                const uint64_t segmentId = journey[position];

                // Locate the Cluster corresponding to this segment,
                // creating it if necessary.
                auto it = clusters.find(segmentId);
                if(it == clusters.end()) {
                    tie(it, ignore) = clusters.insert(make_pair(segmentId, vector<Cluster>(1)));
                }
                SHASTA_ASSERT(it->second.size() == 1);

                // Add this journey entry to the cluster.
                Cluster& cluster = it->second.front();
                cluster.journeyEntryIds.push_back(journeyEntryId);

                clusterTable[orientedReadId.getValue()][position] = &cluster;
            }
        }
    }
}

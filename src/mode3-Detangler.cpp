#include "mode3-Detangler.hpp"
#include "mode3.hpp"
using namespace shasta;
using namespace mode3;



Detangler::Detangler(const AssemblyGraph& assemblyGraph)
{
    createJourneys(assemblyGraph);
    createInitialClusters();
    cout << "The initial Detangler has " << clusters.size() << " clusters." << endl;

    for(const auto& p: clusters) {
        for(const Cluster& cluster: p.second) {
            simpleDetangle(&cluster);
        }
    }
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
            journey.push_back(Step(assemblyGraphJourneyEntry.segmentId));
        }
    }
}



// Initially, we create a Cluster for each segmentId.
void Detangler::createInitialClusters()
{

    // Loop over all oriented reads.
    const ReadId readCount = ReadId(journeys.size() / 2);
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);

            // Get the Journey for this oriented read.
            Journey& journey = journeys[orientedReadId.getValue()];

            // Loop over Step(s) in this Journey.
            StepInfo stepInfo;
            stepInfo.orientedReadId = orientedReadId;
            for(uint64_t position=0; position<journey.size(); position++) {
                stepInfo.position = position;
                Step& step = journey[position];
                const uint64_t segmentId = step.segmentId;

                // Locate the Cluster corresponding to this segment,
                // creating it if necessary.
                ClusterContainer::iterator it = clusters.find(segmentId);
                if(it == clusters.end()) {
                    tie(it, ignore) = clusters.insert(make_pair(segmentId, std::list<Cluster>()));
                    it->second.push_back(Cluster(segmentId));
                }
                std::list<Cluster>& segmentClusters = it->second;

                // Sanity check: this segmentId must correspond to exactly one Cluster.
                SHASTA_ASSERT(segmentClusters.size() == 1);
                Cluster& cluster = segmentClusters.front();

                // Add this Step to the Cluster.
                cluster.steps.push_back(stepInfo);
                step.cluster = &cluster;
            }
        }
    }
}



// Find the next/previous cluster for each of the steps in a given cluster.
// The output vector has size equal to the number of steps in this cluster,
// and the corresponding OrientedReadId(s) are the same
// as the ones in the steps vector for the given cluster.
// Some of the pointers returned can be zero. This can happen if this
// cluster is the first or last cluster in the journey of an oriented read.
void Detangler::findNextClusters(
    const Cluster* cluster0,
    vector<const Cluster*>& nextClusters
    ) const
{
    nextClusters.clear();

    // Loop over the steps of this cluster.
    for(const StepInfo& stepInfo: cluster0->steps) {
        const OrientedReadId orientedReadId = stepInfo.orientedReadId;
        const uint64_t position = stepInfo.position;

        // Get journey for this oriented read.
        const Journey& journey = journeys[orientedReadId.getValue()];

        // Locate the cluster at the next position in the journey.
        // There is none if we are at the end of the journey.
        const Cluster* cluster1 = 0;
        const uint64_t nextPosition = position + 1;
        if(nextPosition < journey.size()) {
            cluster1 = journey[nextPosition].cluster;
        }

        // Store it in the output vector.
        nextClusters.push_back(cluster1);
    };

}
void Detangler::findPreviousClusters(
    const Cluster* cluster0,
    vector<const Cluster*>& previousClusters
    ) const
{
    previousClusters.clear();

    // Loop over the steps of this cluster.
    for(const StepInfo& stepInfo: cluster0->steps) {
        const OrientedReadId orientedReadId = stepInfo.orientedReadId;
        const uint64_t position = stepInfo.position;

        // Get the journey for this oriented read.
        const Journey& journey = journeys[orientedReadId.getValue()];

        // Locate the cluster at the previous position in the journey.
        // There is none if we are at the end of the journey.
        const Cluster* cluster1 = 0;
        if(position > 0) {
            const uint64_t previousPosition = position - 1;
            cluster1 = journey[previousPosition].cluster;
        }

        // Store it in the output vector.
        previousClusters.push_back(cluster1);
    };

}



// Simple, classical detangling of a single cluster.
void Detangler::simpleDetangle(const Cluster* cluster)
{
    // ****** EXPOSE WHEN CODE STABILIZES
    const uint64_t minLinkCoverage = 6;

    std::map< pair<const Cluster*, const Cluster*>, uint64_t> tangleMatrix;

    std::map<const Cluster*, uint64_t> previousCount;
    std::map<const Cluster*, uint64_t> nextCount;

    // Loop over the Step(s) in this Cluster.
    for(const StepInfo& stepInfo: cluster->steps) {
        const OrientedReadId orientedReadId = stepInfo.orientedReadId;
        const uint64_t position = stepInfo.position;
        const Journey& journey = journeys[orientedReadId.getValue()];

        // Increment counts.
        if(position != 0) {
            const Cluster* previousCluster = journey[position - 1].cluster;
            auto it = previousCount.find(previousCluster);
            if(it == previousCount.end()) {
                previousCount.insert(make_pair(previousCluster, 1));
            } else {
                ++(it->second);
            }
        }
        if(position != journey.size()-1) {
            const Cluster* nextCluster = journey[position + 1].cluster;
            auto it = nextCount.find(nextCluster);
            if(it == nextCount.end()) {
                nextCount.insert(make_pair(nextCluster, 1));
            } else {
                ++(it->second);
            }
        }

        // If at the beginning or end of the journey, we cannot use it.
        if(position == 0) {
            continue;
        }
        if(position == journey.size()-1) {
            continue;
        }

        // Access the previous and previous Step.
        const Step& previousStep = journey[position - 1];
        const Step& nextStep = journey[position + 1];

        // Get the previous and next Cluster.
        const Cluster* previousCluster = previousStep.cluster;
        const Cluster* nextCluster = nextStep.cluster;

        // Increment the tangle matrix for this pair.
        const auto p = make_pair(previousCluster, nextCluster);
        auto it = tangleMatrix.find(p);
        if(it == tangleMatrix.end()) {
            tangleMatrix.insert(make_pair(p, 1));
        } else {
            ++(it->second);
        }

    };


    if(previousCount.size() < 2 or nextCount.size() < 2) {
        return;
    }

#if 0
    // Write out the complete detangle matrix.
    {

        cout << "Coverage for links from previous clusters:\n";
        for(const auto& p: previousCount) {
            cout << p.first->stringId() << " " << p.second << "\n";
        }
        cout << "Coverage for links to next clusters:\n";
        for(const auto& p: nextCount) {
            cout << p.first->stringId() << " " << p.second << "\n";
        }


        cout << "Tangle matrix:\n";
        for(const auto& previous: previousCount) {
            const Cluster* previousCluster = previous.first;
            for(const auto& next: nextCount) {
                const Cluster* nextCluster = next.first;
                cout << previousCluster->stringId()  << " ";
                cout << nextCluster->stringId()   << " ";
                auto it = tangleMatrix.find(make_pair(previousCluster, nextCluster));
                if(it == tangleMatrix.end()) {
                    cout << "0";
                } else {
                    cout << it->second;
                }
                cout << "\n";
            }
        }
    }
#endif


    // For detangling, ignore incoming/outgoing links
    // with coverage less than minLinkCoverage.
    for(auto it=previousCount.begin(); /* Later */ ; /* Later */) {

        // Save an incremented iterator for later.
        // We need to do this because it may be invalidated by
        // the next erase.
        auto itNext = it;
        ++itNext;

        // If low coverage, erase.
        if(it->second < minLinkCoverage) {
            previousCount.erase(it);
        }

        // Increment the iterator and check if done.
        it = itNext;
        if(it == previousCount.end()) {
            break;
        }
    }
    for(auto it=nextCount.begin(); /* Later */ ; /* Later */) {

        // Save an incremented iterator for later.
        // We need to do this because it may be invalidated by
        // the next erase.
        auto itNext = it;
        ++itNext;

        // If low coverage, erase.
        if(it->second < minLinkCoverage) {
            nextCount.erase(it);
        }

        // Increment the iterator and check if done.
        it = itNext;
        if(it == nextCount.end()) {
            break;
        }
    }


    // Only detangle if we still have at least two incoming and
    // two outgoing links.
    if(previousCount.size() < 2 or nextCount.size() < 2) {
        return;
    }

    // For now, only detangle the 2 by 2 case.
    if(not(previousCount.size() == 2 and nextCount.size() ==2)) {
        return;
    }


    if(true) {
        cout << "Detangling Cluster " << cluster->stringId()  << "\n";

        cout << "Coverage for links from previous clusters:\n";
        for(const auto& p: previousCount) {
            cout << p.first->stringId() << " " << p.second << "\n";
        }
        cout << "Coverage for links to next clusters:\n";
        for(const auto& p: nextCount) {
            cout << p.first->stringId() << " " << p.second << "\n";
        }

        cout << "Tangle matrix:\n";
        for(const auto& previous: previousCount) {
            const Cluster* previousCluster = previous.first;
            for(const auto& next: nextCount) {
                const Cluster* nextCluster = next.first;
                cout << previousCluster->stringId()  << " ";
                cout << nextCluster->stringId()   << " ";
                auto it = tangleMatrix.find(make_pair(previousCluster, nextCluster));
                if(it == tangleMatrix.end()) {
                    cout << "0";
                } else {
                    cout << it->second;
                }
                cout << "\n";
            }
        }
    }

}

#include "mode3-Detangler.hpp"
#include "deduplicate.hpp"
#include "mode3.hpp"
using namespace shasta;
using namespace mode3;



Detangler::Detangler(const AssemblyGraph& assemblyGraph)
{
    createJourneys(assemblyGraph);
    createInitialClusters();
    cout << "The initial Detangler has " << clusters.size() << " clusters." << endl;

    uint64_t count = 0;
    for(auto& p: clusters) {
        for(Cluster& cluster: p.second) {
            if(simpleDetangle(&cluster)) {
                ++count;
            }
        }
    }
    cout << "Detangled " << count << " clusters out of " << clusters.size() << endl;
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
                    it->second.push_back(Cluster(segmentId, 0));
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
bool Detangler::simpleDetangle(Cluster* cluster0)
{
    // ****** EXPOSE WHEN CODE STABILIZES
    const uint64_t minLinkCoverage = 6;
    const uint64_t maxDiscordantCount = 2;
    const uint64_t minConcordantCount = 8;

    const bool debug = true;

    // Find the previous clusters for each of the steps in this cluster.
    vector<const Cluster*> previousClusters;
    findPreviousClusters(cluster0, previousClusters);
    SHASTA_ASSERT(previousClusters.size() == cluster0->steps.size());


    // Find the next clusters for each of the steps in this cluster.
    vector<const Cluster*> nextClusters;
    findNextClusters(cluster0, nextClusters);
    SHASTA_ASSERT(nextClusters.size() == cluster0->steps.size());

    // Count the distinct previous clusters.
    // They are stored sorted.
    vector<const Cluster*> distinctPreviousClusters = previousClusters;
    vector<uint64_t > distinctPreviousClustersCoverage;
    deduplicateAndCount(distinctPreviousClusters, distinctPreviousClustersCoverage);
    SHASTA_ASSERT(distinctPreviousClusters.size() == distinctPreviousClustersCoverage.size());

    // If less than two, do nothing.
    if(distinctPreviousClusters.size() < 2) {
        return false;
    }

    // Count the distinct previous clusters.
    // They are stored sorted.
    vector<const Cluster*> distinctNextClusters = nextClusters;
    vector<uint64_t > distinctNextClustersCoverage;
    deduplicateAndCount(distinctNextClusters, distinctNextClustersCoverage);
    SHASTA_ASSERT(distinctNextClusters.size() == distinctNextClustersCoverage.size());

    // If less than two, do nothing.
    if(distinctPreviousClusters.size() < 2) {
        return false;
    }

    // Only keep the previous clusters that have sufficient coverage and are not null.
    vector< pair<const Cluster*, uint64_t> > previousWithCoverage;
    for(uint64_t i=0; i<distinctPreviousClusters.size(); i++) {
        const Cluster* cluster1 = distinctPreviousClusters[i];
        if(cluster1) {
            const uint64_t coverage = distinctPreviousClustersCoverage[i];
            if(coverage >= minLinkCoverage) {
                previousWithCoverage.push_back(make_pair(cluster1, coverage));
            }
        }
    }

    // Only keep the next clusters that have sufficient coverage and are not null.
    vector< pair<const Cluster*, uint64_t> > nextWithCoverage;
    for(uint64_t i=0; i<distinctNextClusters.size(); i++) {
        const Cluster* cluster1 = distinctNextClusters[i];
        if(cluster1) {
            const uint64_t coverage = distinctNextClustersCoverage[i];
            if(coverage >= minLinkCoverage) {
                nextWithCoverage.push_back(make_pair(cluster1, coverage));
            }
        }
    }

    // Compute the tangle matrix.
    // tangleMatrix[i][j] contains the number of oriented reads
    // that come from the i-th previous cluster and go to the j-th previous cluster.
    vector< vector<uint64_t> > tangleMatrix(previousWithCoverage.size(), vector<uint64_t>(nextWithCoverage.size(), 0));
    for(uint64_t i=0; i<previousWithCoverage.size(); i++) {
        const Cluster* previousCluster = previousWithCoverage[i].first;
        for(uint64_t j=0; j<nextWithCoverage.size(); j++) {
            const Cluster* nextCluster = nextWithCoverage[j].first;
            for(uint64_t k=0; k<previousClusters.size(); k++) {
                if((previousClusters[k] == previousCluster) and (nextClusters[k] == nextCluster)) {
                    ++tangleMatrix[i][j];
                }
            }
        }
    }

    // For now, only handle the 2 by 2 case.
    if(not(previousWithCoverage.size() == 2 and nextWithCoverage.size() == 2)) {
        return false;
    }

    // Compute the sum of diagonal and off-diagonal terms.
    const uint64_t diagonalSum = tangleMatrix[0][0] + tangleMatrix[1][1];
    const uint64_t offDiagonalSum = tangleMatrix[0][1] + tangleMatrix[1][0];

    // Check if the criteria for detangle are satisfied.
    const uint64_t concordantCount = max(diagonalSum, offDiagonalSum);
    const uint64_t discordantCount = min(diagonalSum, offDiagonalSum);
    if(concordantCount < minConcordantCount or discordantCount > maxDiscordantCount) {
        return false;
    }

    if(debug) {
        cout << "Detangling " << cluster0->stringId() << "\n";
        cout << "Previous:\n";
        for(const auto& p: previousWithCoverage) {
            cout << p.first->stringId() << " " << p.second << "\n";
        }
        cout << "Next:\n";
        for(const auto& p: nextWithCoverage) {
            cout << p.first->stringId() << " " << p.second << "\n";
        }
        cout << "Tangle matrix:\n";
        for(uint64_t i=0; i<previousWithCoverage.size(); i++) {
            const Cluster* previousCluster = previousWithCoverage[i].first;
            for(uint64_t j=0; j<nextWithCoverage.size(); j++) {
                const Cluster* nextCluster = nextWithCoverage[j].first;
                cout << previousCluster->stringId() << " ";
                cout << nextCluster->stringId() << " ";
                cout << tangleMatrix[i][j] << "\n";
            }
        }
        cout << "Diagonal " << diagonalSum << "\n";
        cout << "Off-diagonal " << offDiagonalSum << "\n";

    }



    // If getting here, we can detangle this cluster.
    // This generates two new clusters for this segment.
    const bool inPhase = diagonalSum > offDiagonalSum;

    // The new steps for cluster0.
    vector<StepInfo> newSteps0;

    // Create the two new clusters.
    const uint64_t segmentId = cluster0->segmentId;
    std::list<Cluster>& segmentClusters = clusters[segmentId];
    segmentClusters.push_back(Cluster(segmentId, segmentClusters.size()));
    Cluster& cluster1 = segmentClusters.back();
    segmentClusters.push_back(Cluster(segmentId, segmentClusters.size()));
    Cluster& cluster2 = segmentClusters.back();

    // Do the detangling. The steps that correspond to the dominant portion of the
    // tangle matrix are moved to the new clusters.
    for(uint64_t k=0; k<previousClusters.size(); k++) {
        const StepInfo& step = cluster0->steps[k];
        const OrientedReadId orientedReadId = step.orientedReadId;
        Journey& journey = journeys[orientedReadId.getValue()];
        const uint64_t position = step.position;
        const Cluster* previousCluster = previousClusters[k];
        const Cluster* nextCluster = nextClusters[k];
        if(inPhase) {
            if(previousCluster == previousWithCoverage[0].first and nextCluster == nextWithCoverage[0].first) {
                // Add it to the steps of cluster1.
                cluster1.steps.push_back(StepInfo(orientedReadId, position));
                journey[position].cluster = &cluster1;
            } else if(previousCluster == previousWithCoverage[1].first and nextCluster == nextWithCoverage[1].first) {
                // Add it to the steps of cluster2.
                cluster2.steps.push_back(StepInfo(orientedReadId, position));
                journey[position].cluster = &cluster2;
            } else {
                // Leave it in cluster0.
                newSteps0.push_back(StepInfo(orientedReadId, position));
            }
        } else {
            if(previousCluster == previousWithCoverage[0].first and nextCluster == nextWithCoverage[1].first) {
                // Add it to the steps of cluster1.
                cluster1.steps.push_back(StepInfo(orientedReadId, position));
                journey[position].cluster = &cluster1;
            } else if(previousCluster == previousWithCoverage[1].first and nextCluster == nextWithCoverage[0].first) {
                // Add it to the steps of cluster2.
                cluster2.steps.push_back(StepInfo(orientedReadId, position));
                journey[position].cluster = &cluster2;
            } else {
                // Leave it in cluster0.
                newSteps0.push_back(StepInfo(orientedReadId, position));
            }
        }
   }




    // Update the steps of the cluster we just detangled.
    cluster0->steps.swap(newSteps0);

    return true;
}


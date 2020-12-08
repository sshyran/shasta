#include "Align5.hpp"
#include "orderPairs.hpp"
#include "PngImage.hpp"
#include "SHASTA_ASSERT.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace Align5;

#include "fstream.hpp"
#include "tuple.hpp"



void shasta::align5(
    const span<const CompressedMarker>& markers0,
    const span<const CompressedMarker>& markers1,
    const Options& options,
    MemoryMapped::VectorOfVectors<Align5::MatrixEntry, uint64_t>& matrix, // Used as work area.
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    bool debug)
{
    switch(options.m) {
    case 1:
        align5<1>(markers0, markers1, options, matrix, alignment, alignmentInfo, debug);
        return;
    case 2:
        align5<2>(markers0, markers1, options, matrix, alignment, alignmentInfo, debug);
        return;
    case 3:
        align5<3>(markers0, markers1, options, matrix, alignment, alignmentInfo, debug);
        return;
    case 4:
        align5<4>(markers0, markers1, options, matrix, alignment, alignmentInfo, debug);
        return;
    default:
        SHASTA_ASSERT(0);
    }
}



// Version templated on m, the number of markers that define
// a "feature" used in the alignment.
template<uint64_t m> void shasta::align5(
    const span<const CompressedMarker>& markers0,
    const span<const CompressedMarker>& markers1,
    const Options& options,
    MemoryMapped::VectorOfVectors<Align5::MatrixEntry, uint64_t>& matrix, // Used as work area.
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    bool debug)
{
    Align5::Aligner<m> graph(markers0, markers1,
        options, matrix, alignment, alignmentInfo,
        debug);
}



template<uint64_t m> shasta::Align5::Aligner<m>::Aligner(
    const MarkerSequence& markerSequence0,
    const MarkerSequence& markerSequence1,
    const Options& options,
    MemoryMapped::VectorOfVectors<MatrixEntry, uint64_t>&, // Currently unused
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    bool debug) :
    nx(uint32_t(markerSequence0.size())),
    ny(uint32_t(markerSequence1.size())),
    deltaX(int32_t(options.deltaX)),
    deltaY(int32_t(options.deltaY))
{
    // ***************************** EXPOSE WHEN CODE STABILIZES
    const uint32_t minEntryCountPerCell = 10;

    if(debug) {
        cout << timestamp << "Align5 begins." << endl;
        cout << timestamp << "Input sequences have " <<
            nx << " and " << ny << " markers." << endl;
    }

    // Check that we are in the templated version consistent with
    // the options.
    SHASTA_ASSERT(options.m == m);

    // Create features sorted by feature.
    if(debug) {
        cout << timestamp << "Creating sorted features." << endl;
    }
    sortFeatures(markerSequence0, sortedFeatures0);
    sortFeatures(markerSequence1, sortedFeatures1);

    // Create the sparse representation of the alignment matrix.
    if(debug) {
        cout << timestamp << "Creating the alignment matrix." << endl;
    }
    createAlignmentMatrix();
    if(debug) {
        cout << timestamp << "Writing the alignment matrix." << endl;
        writeAlignmentMatrix("Align5-AlignmentMatrix.csv");
    }

    // Gather well populated cells.
    if(debug) {
        cout << timestamp << "Gathering well populated cells." << endl;
    }
    createCells(minEntryCountPerCell);
    if(debug) {
        cout << timestamp << "Writing cells." << endl;
        writeCellsCsv("Align5-Cells.csv");
        writeCellsPng("Align5-Cells.png");
    }

    if(debug) {
        cout << timestamp << "Writing alignment matrix in feature space." << endl;
        writeAlignmentMatrixInFeatureSpace("Align5-AlignmentMatrixInFeatureSpace.png");
    }


    // Create markers sorted by KmerId.
    if(debug) {
        cout << timestamp << "Creating sorted markers." << endl;
    }
    sortMarkers(markerSequence0, sortedMarkers0);
    sortMarkers(markerSequence1, sortedMarkers1);

    if(debug) {
        cout << timestamp << "Writing alignment matrix in marker space." << endl;
        writeAlignmentMatrixInMarkerSpace("Align5-AlignmentMatrixInMarkerSpace.png");
    }

    if(debug) {
        cout << timestamp << "Align5 ends." << endl;
    }
}



template<uint64_t m> void shasta::Align5::Aligner<m>::sortFeatures(
    const MarkerSequence& markerSequence,
    vector< pair<Feature, uint32_t> >& sortedFeatures)
{
    sortedFeatures.clear();

    // If the sequence is too short, it has no features.
    if(markerSequence.size() < m) {
        return;
    }
    sortedFeatures.reserve(markerSequence.size() - (m-1));

    // Start with the first m kmerIds.
    Feature feature;
    uint64_t nextOrdinal = 0;
    for(uint64_t i=0; i<m; i++) {
        feature[i] = markerSequence[nextOrdinal++].kmerId;
    }


    // Add the features.
    for(uint32_t ordinal=0; ; ordinal++) {
        sortedFeatures.push_back(make_pair(feature, ordinal));

        // Check if done.
        if(nextOrdinal == markerSequence.size()) {
            break;
        }

        // Shift by one.
        for(uint64_t i=1; i<m; i++) {
            feature[i-1] = feature[i];
        }
        feature.back() = markerSequence[nextOrdinal++].kmerId;
    }
    SHASTA_ASSERT(sortedFeatures.size() == markerSequence.size() - (m-1));

    // Sort by feature.
    sort(sortedFeatures.begin(), sortedFeatures.end(),
        OrderPairsByFirstOnly<Feature, uint32_t>());
}



template<uint64_t m> void shasta::Align5::Aligner<m>::sortMarkers(
    const MarkerSequence& markerSequence,
    vector< pair<KmerId, uint32_t> >& sortedMarkers)
{
    sortedMarkers.resize(markerSequence.size());
    for(uint32_t i=0; i<markerSequence.size(); i++) {
        sortedMarkers[i] = make_pair(markerSequence[i].kmerId, i);
    }
    sort(sortedMarkers.begin(), sortedMarkers.end(),
        OrderPairsByFirstOnly<KmerId, uint32_t>());
}



template<uint64_t m> void shasta::Align5::Aligner<m>::
    writeAlignmentMatrixInMarkerSpace(const string& fileName) const
{
    PngImage image(nx, ny);
    writeCheckerboard(image);
    image.writeGrid(   10,  15,  15,  15);      // Grey
    image.writeGrid(   50,  30,  30,  30);      // Grey
    image.writeGrid(  100,  90,  90,  90);      // Grey
    image.writeGrid(  500, 160, 160, 160);      // Grey
    image.writeGrid( 1000, 255, 255, 255);      // White
    image.writeGrid( 5000, 255, 120, 255);      // Purple
    image.writeGrid(10000, 255, 255,  60);      // Yellow
    image.writeGrid(50000, 255, 255, 120);      // Yellow



    // Joint loop over the markers, looking for common k-mer ids.
    uint64_t n = 0;
    auto begin0 = sortedMarkers0.begin();
    auto begin1 = sortedMarkers1.begin();
    auto end0 = sortedMarkers0.end();
    auto end1 = sortedMarkers1.end();

    auto it0 = begin0;
    auto it1 = begin1;
    while(it0!=end0 && it1!=end1) {
        if(it0->first < it1->first) {
            ++it0;
        } else if(it1->first < it0->first) {
            ++it1;
        } else {

            // We found a common k-mer id.
            const KmerId kmerId = it0->first;


            // This k-mer could appear more than once in each of the oriented reads,
            // so we need to find the streak of this k-mer in kmers0 and kmers1.
            auto it0Begin = it0;
            auto it1Begin = it1;
            auto it0End = it0Begin;
            auto it1End = it1Begin;
            while(it0End!=end0 && it0End->first == kmerId) {
                ++it0End;
            }
            while(it1End!=end1 && it1End->first == kmerId) {
                ++it1End;
            }

            // Loop over pairs in the streaks.
            for(auto jt0=it0Begin; jt0!=it0End; ++jt0) {
                for(auto jt1=it1Begin; jt1!=it1End; ++jt1) {
                    const uint32_t x = jt0->second;
                    const uint32_t y = jt1->second;
                    image.setPixel(x, y, 255, 0, 0);
                    ++n;
                }
            }

            // Continue the joint loop over k-mers.
            it0 = it0End;
            it1 = it1End;
        }

    }


    image.write(fileName);
    cout << "The alignment matrix in marker space contains " <<
        n << " entries." << endl;
}



template<uint64_t m> void shasta::Align5::Aligner<m>::
    writeAlignmentMatrixInFeatureSpace(const string& fileName) const
{
    PngImage image(nx, ny);
    writeCheckerboard(image);
    image.writeGrid(   10,  15,  15,  15);      // Grey
    image.writeGrid(   50,  30,  30,  30);      // Grey
    image.writeGrid(  100,  90,  90,  90);      // Grey
    image.writeGrid(  500, 160, 160, 160);      // Grey
    image.writeGrid( 1000, 255, 255, 255);      // White
    image.writeGrid( 5000, 255, 120, 255);      // Purple
    image.writeGrid(10000, 255, 255,  60);      // Yellow
    image.writeGrid(50000, 255, 255, 120);      // Yellow



    // Joint loop over the features, looking for common k-mer ids.
    uint64_t n = 0;
    auto begin0 = sortedFeatures0.begin();
    auto begin1 = sortedFeatures1.begin();
    auto end0 = sortedFeatures0.end();
    auto end1 = sortedFeatures1.end();

    auto it0 = begin0;
    auto it1 = begin1;
    while(it0!=end0 && it1!=end1) {
        if(it0->first < it1->first) {
            ++it0;
        } else if(it1->first < it0->first) {
            ++it1;
        } else {

            // We found a common k-mer id.
            const Feature feature = it0->first;


            // This k-mer could appear more than once in each of the oriented reads,
            // so we need to find the streak of this k-mer in kmers0 and kmers1.
            auto it0Begin = it0;
            auto it1Begin = it1;
            auto it0End = it0Begin;
            auto it1End = it1Begin;
            while(it0End!=end0 && it0End->first == feature) {
                ++it0End;
            }
            while(it1End!=end1 && it1End->first == feature) {
                ++it1End;
            }

            // Loop over pairs in the streaks.
            for(auto jt0=it0Begin; jt0!=it0End; ++jt0) {
                for(auto jt1=it1Begin; jt1!=it1End; ++jt1) {
                    const uint32_t x = jt0->second;
                    const uint32_t y = jt1->second;
                    image.setPixel(x, y, 255, 0, 0);
                    ++n;
                }
            }

            // Continue the joint loop over features.
            it0 = it0End;
            it1 = it1End;
        }

    }


    image.write(fileName);
    cout << "The alignment matrix in feature space contains " <<
        n << " entries." << endl;
}



template<uint64_t m> void
    shasta::Align5::Aligner<m>::writeCheckerboard(PngImage& image) const
{
    Coordinates xy;
    uint32_t& x = xy.first;
    uint32_t& y = xy.second;

    Coordinates iXY;
    uint32_t& iX = iXY.first;
    uint32_t& iY = iXY.second;

    for(y=0; y<ny; y++) {
        for(x=0; x<nx; x++) {
            iXY = getCellIndexesFromxy(xy);
            if(((iX + iY) %2) == 0) {
                image.setPixel(x, y, 0, 48, 0);
            }
        }
    }
}



// Return (X,Y) given (x,y).
template<uint64_t m> shasta::Align5::Coordinates
    shasta::Align5::Aligner<m>::getXY(Coordinates xy) const
{
    return Coordinates(
        xy.first + xy.second,
        nx + xy.second - xy.first - 1
        );
}



// Return (iX,iY) given (X,Y).
template<uint64_t m> shasta::Align5::Coordinates
    shasta::Align5::Aligner<m>::getCellIndexesFromXY(Coordinates XY) const
{
    return Coordinates(
        XY.first  / deltaX,
        XY.second / deltaY
        );
}


// Return (iX,iY) given (x,y).
template<uint64_t m> shasta::Align5::Coordinates
    shasta::Align5::Aligner<m>::getCellIndexesFromxy(Coordinates xy) const
{
    const Coordinates XY = getXY(xy);
    return getCellIndexesFromXY(XY);
}



template<uint64_t m> void shasta::Align5::Aligner<m>::createAlignmentMatrix()
{
    alignmentMatrix.clear();

    // Joint loop over the sorted features, looking for common features.
    auto begin0 = sortedFeatures0.begin();
    auto begin1 = sortedFeatures1.begin();
    auto end0 = sortedFeatures0.end();
    auto end1 = sortedFeatures1.end();

    auto it0 = begin0;
    auto it1 = begin1;
    while(it0!=end0 && it1!=end1) {
        if(it0->first < it1->first) {
            ++it0;
        } else if(it1->first < it0->first) {
            ++it1;
        } else {

            // We found a common feature.
            const Feature feature = it0->first;


            // This feature could appear more than once in each of the sequences,
            // so we need to find the streak of this feature.
            auto it0Begin = it0;
            auto it1Begin = it1;
            auto it0End = it0Begin;
            auto it1End = it1Begin;
            while(it0End!=end0 && it0End->first == feature) {
                ++it0End;
            }
            while(it1End!=end1 && it1End->first == feature) {
                ++it1End;
            }

            // Loop over pairs in the streaks.
            for(auto jt0=it0Begin; jt0!=it0End; ++jt0) {
                for(auto jt1=it1Begin; jt1!=it1End; ++jt1) {
                    const uint32_t x = jt0->second;
                    const uint32_t y = jt1->second;
                    const Coordinates iXY = getCellIndexesFromxy(Coordinates(x, y));
                    const uint32_t iX = iXY.first;
                    const uint32_t iY = iXY.second;
                    if(alignmentMatrix.size() <= iY) {
                        alignmentMatrix.resize(iY+1);
                    }
                    alignmentMatrix[iY].push_back(make_pair(iX, Coordinates(x, y)));
                }
            }

            // Continue the joint loop over features.
            it0 = it0End;
            it1 = it1End;
        }

    }

    for(auto& v: alignmentMatrix) {
        sort(v.begin(), v.end(), OrderPairsByFirstOnly<uint32_t, Coordinates>());
    }
}



template<uint64_t m> void shasta::Align5::Aligner<m>::writeAlignmentMatrix(const string& fileName) const
{
    uint64_t entryCount = 0;
    ofstream csv(fileName);
    csv << "iX,iY,X,Y,x,y\n";
    for(uint32_t iY=0; iY<alignmentMatrix.size(); iY++) {
        for(const auto& v: alignmentMatrix[iY]) {
            const uint32_t iX = v.first;
            const Coordinates& xy = v.second;
            const Coordinates XY = getXY(xy);
            csv << iX << ",";
            csv << iY << ",";
            csv << XY.first << ",";
            csv << XY.second << ",";
            csv << xy.first << ",";
            csv << xy.second << "\n";
            SHASTA_ASSERT(Coordinates(iX, iY) == getCellIndexesFromxy(xy));
            SHASTA_ASSERT(Coordinates(iX, iY) == getCellIndexesFromXY(XY));
            ++entryCount;
        }
    }
    cout << "The alignment matrix contains " << entryCount << " entries." << endl;
}



template<uint64_t m> void shasta::Align5::Aligner<m>::createCells(
    uint32_t minEntryCountPerCell)
{
    // Start with nothing.
    cells.clear();
    cells.resize(alignmentMatrix.size());

    // Loop over iY values.
    for(uint32_t iY=0; iY<alignmentMatrix.size(); iY++) {

        // Access the vectors for this value of iY.
        const vector< pair<uint32_t, Coordinates> >& iYAlignmentMatrix = alignmentMatrix[iY];
        vector< pair<uint32_t, Cell> >& iYCells = cells[iY];

        // Each vector in the alignment matrix is sorted by iX, so we can scan it,
        // creating a new cell each time we encounter a new value of iX
        // (but only if the cell is sufficiently populated).
        for(auto it=iYAlignmentMatrix.begin(); it!=iYAlignmentMatrix.end(); /* Increment later */) {

            // When getting here, we are starting a new cell.
            const auto& entry = *it;
            const uint32_t iX = entry.first;

            // Now count the the alignment matrix entries in this cell.
            const auto it0 = it;
            while(true) {
                ++it;
                if(it == iYAlignmentMatrix.end()) {
                    break;
                }
                if(it->first != iX) {
                    break;
                }
            }

            // If too few entries in this cell, don't store it.
            if(it - it0 < minEntryCountPerCell) {
                continue;
            }

            // Store this cell.
            iYCells.push_back(make_pair(iX, Cell()));
        }

    }
}



template<uint64_t m> void shasta::Align5::Aligner<m>::writeCellsCsv(
    const string& fileName) const
{
    uint64_t cellCount = 0;
    ofstream csv(fileName);
    csv << "iX,iY,minX,maxX,minY,maxY,sizeX,sizeY\n";
    for(uint32_t iY=0; iY<cells.size(); iY++) {
        for(const auto& p: cells[iY]) {
            const uint32_t iX = p.first;
            csv << iX << ",";
            csv << iY << "\n";
            ++cellCount;
        }
    }
    cout << "There are " << cellCount << " well populated cells." << endl;
}



template<uint64_t m> void shasta::Align5::Aligner<m>::writeCellsPng(
    const string& fileName) const
{
    // The size of the square in XY space.
    const uint32_t sizeXY = nx + ny -1;

    // Number of markers per pixel - to control the size of the picture.
    const uint32_t markersPerPixel = 5;

    // Create the image.
    const uint32_t imageSize = sizeXY / markersPerPixel;
    PngImage image(imageSize, imageSize);

    // Write the checkerboard in XY space.
    for(uint32_t j=0; j<imageSize; j++) {
        const uint32_t Y = j * markersPerPixel;
        for(uint32_t i=0; i<imageSize; i++) {
            const uint32_t X = i * markersPerPixel;
            const Coordinates iXY = getCellIndexesFromXY(Coordinates(X, Y));
            const uint32_t iX = iXY.first;
            const uint32_t iY = iXY.second;
            if(((iX+iY) %2) == 0) {
                image.setPixel(i, j, 0, 50, 0);
            }
        }
    }

    // Write the cells.
    for(uint32_t iY=0; iY<cells.size(); iY++) {
        for(const auto& p: cells[iY]) {
            const uint32_t iX = p.first;
            SHASTA_ASSERT(iX < sizeXY);
            SHASTA_ASSERT(iY < sizeXY);
            const uint32_t iMin = ( iX    * deltaX)  / markersPerPixel;
            const uint32_t iMax = ((iX+1) * deltaX)  / markersPerPixel;
            const uint32_t jMin = ( iY    * deltaY)  / markersPerPixel;
            const uint32_t jMax = ((iY+1) * deltaY)  / markersPerPixel;
            for(uint32_t j=jMin; j<jMax; j++) {
                for(uint32_t i=iMin; i<iMax; i++) {
                    SHASTA_ASSERT(i < imageSize);
                    SHASTA_ASSERT(j < imageSize);
                    image.setPixel(i, j, 255, 255, 255);
                }
            }
         }
    }

    // Write the image.
    image.write(fileName);
}


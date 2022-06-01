#include <assert.h>

#include "hpc.h"

/**
 * @brief Get offset in global array of slice i
 *
 * @param globalSize
 * @param nofSlices
 * @param i
 * @return index Offset of slice i
 */
index getSliceOffset(index globalSize, index nofSlices, index i) {
    index q = globalSize / nofSlices;
    index r = globalSize % nofSlices;
    return i < r ? i * (q + 1) : r * (q + 1) + (i - r) * q;
}

/**
 * @brief Get size of slice i
 *
 * @param globalSize
 * @param nofSlices
 * @param i
 * @return index Size of slice i
 */
index getSliceSize(index globalSize, index nofSlices, index i) {
    index q = globalSize / nofSlices;
    index r = globalSize % nofSlices;
    return i < r ? q + 1 : q;
}

/**
 * @brief Get slice index of global index j
 *
 * @param globalSize
 * @param nofSlices
 * @param j
 * @return index Slice index of global index j
 */
index getSliceIndex(index globalSize, index nofSlices, index j) {
    index q = globalSize / nofSlices;
    index r = globalSize % nofSlices;
    index i = j < r * (q + 1) ? j / (q + 1) : r + (j - r * (q + 1)) / q;
    assert(i >= 0);
    assert(i <= nofSlices);
    return i;
}

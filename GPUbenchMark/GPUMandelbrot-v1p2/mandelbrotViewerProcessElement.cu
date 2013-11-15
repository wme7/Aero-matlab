/** @file processMandelbrotElement.cu
 *
 * Copyright 2010 The Mathworks, Inc.
 * $Revision: 1$
 * $Date: 2010-11-08$
 */

/** Work out which piece of the global array this thread should operate on */ 
__device__ size_t calculateGlobalIndex() {
    // Which block are we?
    size_t const globalBlockIndex = blockIdx.x + blockIdx.y * gridDim.x;
    // Which thread are we within the block?
    size_t const localThreadIdx = threadIdx.x + blockDim.x * threadIdx.y;
    // How big is each block?
    size_t const threadsPerBlock = blockDim.x*blockDim.y;
    // Which thread are we overall?
    return localThreadIdx + globalBlockIndex*threadsPerBlock;

}

/** The actual Mandelbrot algorithm for a single location */ 
__device__ double doIterations( double const realPart0, 
                                      double const imagPart0, 
                                      double const escapeRadius, 
                                      unsigned int const maxIters ) {
    // Initialise: z = z0
    double const escapeRadius2 = escapeRadius*escapeRadius;
    double realPart = realPart0;
    double imagPart = imagPart0;
    unsigned int count = 0;

    // Loop until escape
    while ( ( count <= maxIters )
            && ((realPart*realPart + imagPart*imagPart) <= escapeRadius2) ) {
        ++count;
        // Update: z = z*z + z0;
        double const oldRealPart = realPart;
        realPart = realPart*realPart - imagPart*imagPart + realPart0;
        imagPart = 2.0*oldRealPart*imagPart + imagPart0;
    }

    // Correct final position for smooth shading
    double const absZ2 = ( realPart*realPart + imagPart*imagPart );
    if (absZ2<escapeRadius2) {
        return double(count) + 1.0 - log( log( escapeRadius2 ) / 2.0 ) / log(2.0);
    } else {
        return double(count) + 1.0 - log( log( absZ2 ) / 2.0 ) / log(2.0);
    }
}


/** Main entry point.
 * Works out where the current thread should read/write to global memory
 * and calls doIterations to do the actual work.
 */
__global__ void processMandelbrotElement( 
                      double * out, 
                      const double * x, 
                      const double * y,
                      const double escapeRadius, 
                      const unsigned int maxIters, 
                      const unsigned int numel ) {
    // Work out which thread we are
    size_t const globalThreadIdx = calculateGlobalIndex();

    // If we're off the end, return now
    if (globalThreadIdx >= numel) {
        return;
    }
    
    // Get our X and Y coords
    double const realPart0 = x[globalThreadIdx];
    double const imagPart0 = y[globalThreadIdx];

    // Run the itearations on this location
    double const count = doIterations( realPart0, imagPart0, escapeRadius, maxIters );
    out[globalThreadIdx] = log( count );
}

/** @internal
 ** @file     vl_siftdescriptor.c
 ** @brief    vl_siftdescriptor - MEX definition
 ** @author   Andrea Vedaldi
 **/

/*
Copyright (C) 2007-12 Andrea Vedaldi and Brian Fulkerson.
All rights reserved.

This file is part of the VLFeat library and is made available under
the terms of the BSD license (see the COPYING file).
*/

#include <mexutils.h>
#include <vl/mathop.h>
#include <vl/covdet.h>

#include <math.h>
#include <assert.h>

/* option codes */
enum {
  opt_magnif,
  opt_float_descriptors,
  opt_norm_thresh,
  opt_verbose
} ;

/* options */
vlmxOption  options [] = {
  {"Magnif",           1,   opt_magnif            },
  {"Verbose",          0,   opt_verbose           },
  {"FloatDescriptors", 0,   opt_float_descriptors },
  {"NormThresh",       1,   opt_norm_thresh       },
  {0,              0,   0                         }
} ;

/** ------------------------------------------------------------------
 ** @internal
 ** @brief Transpose descriptor
 **
 ** @param dst destination buffer.
 ** @param src source buffer.
 **
 ** The function writes to @a dst the transpose of the SIFT descriptor
 ** @a src. The transpose is defined as the descriptor that one
 ** obtains from computing the normal descriptor on the transposed
 ** image.
 **/

VL_INLINE void
transpose_descriptor (float* dst, float* src)
{
  int const BO = 8 ;  /* number of orientation bins */
  int const BP = 4 ;  /* number of spatial bins     */
  int i, j, t ;

  for (j = 0 ; j < BP ; ++j) {
    int jp = BP - 1 - j ;
    for (i = 0 ; i < BP ; ++i) {
      int o  = BO * i + BP*BO * j  ;
      int op = BO * i + BP*BO * jp ;
      dst [op] = src[o] ;
      for (t = 1 ; t < BO ; ++t)
        dst [BO - t + op] = src [t + o] ;
    }
  }
}

/** ------------------------------------------------------------------
 ** @brief MEX entry point
 **/

void
mexFunction(int nout, mxArray *out[],
            int nin, const mxArray *in[])
{
  enum {IN_GRAD=0,IN_FRAMES,IN_END} ;
  enum {OUT_DESCRIPTORS} ;

  int                verbose = 0 ;
  int                opt ;
  int                next = IN_END ;
  mxArray const     *optarg ;

  mxArray           *grad_array ;
  float             *grad ;
  int                M, N ;

  vl_bool            floatDescriptors = 0 ;
  double             magnif = 3 ;
  double             window_size = (double)VL_SIFT_NBP / 2;
  double             norm_thresh = 0 ;
  double            *ikeys = 0 ;
  int                nikeys = 0 ;

  int i,j ;

  VL_USE_MATLAB_ENV ;

  /* -----------------------------------------------------------------
   *                                               Check the arguments
   * -------------------------------------------------------------- */

  if (nin < 2) {
    mexErrMsgTxt("Two arguments required.") ;
  } else if (nout > 1) {
    mexErrMsgTxt("Too many output arguments.");
  }

  if (mxGetNumberOfDimensions (in[IN_GRAD])    != 3              ||
      mxGetClassID            (in[IN_GRAD])    != mxSINGLE_CLASS ||
      mxGetDimensions         (in[IN_GRAD])[0] != 2              ) {
    mexErrMsgTxt("GRAD must be a 2xMxN matrix of class SINGLE.") ;
  }

  if (!vlmxIsMatrix(in[IN_FRAMES], 4, -1)) {
    mexErrMsgTxt("FRAMES must be a 4xN matrix.") ;
  }
  nikeys = mxGetN (in[IN_FRAMES]) ;
  ikeys  = mxGetPr(in[IN_FRAMES]) ;

  while ((opt = vlmxNextOption (in, nin, options, &next, &optarg)) >= 0) {
    switch (opt) {

      case opt_verbose :
        ++ verbose ;
        break ;

      case opt_magnif :
        if (!vlmxIsPlainScalar(optarg) || (magnif = *mxGetPr(optarg)) < 0) {
          mexErrMsgTxt("MAGNIF must be a non-negative scalar.") ;
        }
        break ;

      case opt_float_descriptors :
        floatDescriptors = 1 ;
        break ;

      case opt_norm_thresh :
        if (!vlmxIsPlainScalar(optarg) || (norm_thresh = *mxGetPr(optarg)) < 0) {
          mexErrMsgTxt("NORMTHRESH must be a non-negative scalar.") ;
        }
        break ;

      default :
        abort() ;
    }
  }

  grad_array = mxDuplicateArray(in[IN_GRAD]) ;
  grad = (float*) mxGetData (grad_array) ;
  M    = mxGetDimensions(in[IN_GRAD])[1] ;
  N    = mxGetDimensions(in[IN_GRAD])[2] ;

  /* transpose angles */
  for (i = 1 ; i < 2*M*N ; i+=2) {
    grad [i] = VL_PI/2 - grad [i] ;
  }

  /* -----------------------------------------------------------------
   *                                                            Do job
   * -------------------------------------------------------------- */
  {
    void * descr = 0 ;

    if (verbose) {
      mexPrintf("vl_siftdescriptor: filter settings:\n") ;
      mexPrintf("vl_siftdescriptor:   magnif                = %g\n",
                magnif);
      mexPrintf("vl_siftdescriptor:   num of frames         = %d\n",
                nikeys) ;
      mexPrintf("vl_siftdescriptor:   float descriptor      = %d\n",
                floatDescriptors) ;
      mexPrintf("vl_siftdescriptor:   norm thresh           = %g\n",
                norm_thresh);
    }

    {
      mwSize dims [2] ;
      dims [0] = 128 ;
      dims [1] = nikeys ;
      out[OUT_DESCRIPTORS]= mxCreateNumericArray
        (2, dims,
         floatDescriptors ? mxSINGLE_CLASS : mxUINT8_CLASS,
         mxREAL) ;
      descr = mxGetData(out[OUT_DESCRIPTORS]) ;
    }

    /* ...............................................................
     *                                             Process each octave
     * ............................................................ */
    for (i = 0 ; i < nikeys ; ++i) {
      float  buf [128], rbuf [128] ;

      double y  = *ikeys++ - 1 ;
      double x  = *ikeys++ - 1 ;
      double s  = *ikeys++ ;
      double th = VL_PI / 2 - *ikeys++ ;

      vl_sift_calc_descriptor (grad, M, N, buf, x, y, s, th, magnif,
                               window_size, norm_thresh, 0) ;

      transpose_descriptor (rbuf, buf) ;

      if (! floatDescriptors) {
        vl_uint8 * descr_ = descr ;
        for (j = 0 ; j < 128 ; ++j) {
          float x = 512.0F * rbuf [j] ;
          x = (x < 255.0F) ? x : 255.0F ;
          *descr_++ = (vl_uint8) (x) ;
        }
        descr = descr_ ;
      } else {
        float * descr_ = descr ;
        for (j = 0 ; j < 128 ; ++j) {
          *descr_++ = 512.0F * rbuf [j] ;
        }
        descr = descr_ ;
      }
    }
    /* cleanup */
    mxDestroyArray (grad_array) ;
  } /* job done */
}

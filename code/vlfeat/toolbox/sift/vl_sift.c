/** @internal
 ** @file     sift.c
 ** @author   Andrea Vedaldi
 ** @brief    Scale Invariant Feature Transform (SIFT) - MEX
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
  opt_octaves = 0,
  opt_levels,
  opt_first_octave,
  opt_frames,
  opt_edge_thresh,
  opt_peak_thresh,
  opt_norm_thresh,
  opt_magnif,
  opt_window_size,
  opt_orientations,
  opt_float_descriptors,
  opt_verbose
} ;

/* options */
vlmxOption  options [] = {
  {"Octaves",          1,   opt_octaves           },
  {"Levels",           1,   opt_levels            },
  {"FirstOctave",      1,   opt_first_octave      },
  {"Frames",           1,   opt_frames            },
  {"PeakThresh",       1,   opt_peak_thresh       },
  {"EdgeThresh",       1,   opt_edge_thresh       },
  {"NormThresh",       1,   opt_norm_thresh       },
  {"Magnif",           1,   opt_magnif            },
  {"WindowSize",       1,   opt_window_size       },
  {"Orientations",     0,   opt_orientations      },
  {"FloatDescriptors", 0,   opt_float_descriptors },
  {"Verbose",          0,   opt_verbose           },
  {0,                  0,   0                     }
} ;

/** ------------------------------------------------------------------
 ** @internal
 ** @brief Transpose desriptor
 **
 ** @param dst destination buffer.
 ** @param src source buffer.
 **
 ** The function writes to @a dst the transpose of the SIFT descriptor
 ** @a src. The tranpsose is defined as the descriptor that one
 ** obtains from computing the normal descriptor on the transposed
 ** image.
 **/

VL_INLINE void
transpose_descriptor (float *dst, float const *src)
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
 ** @internal
 ** @brief Ordering of tuples by increasing scale
 **
 ** @param a tuple.
 ** @param b tuble.
 **
 ** @return @c a[2] < b[2]
 **/

static int
korder (void const* a, void const* b) {
  double x = ((double*) a) [2] - ((double*) b) [2] ;
  if (x < 0) return -1 ;
  if (x > 0) return +1 ;
  return 0 ;
}

/** ------------------------------------------------------------------
 ** @internal
 ** @brief Check for sorted keypoints
 **
 ** @param keys keypoint list to check
 ** @param nkeys size of the list.
 **
 ** @return 1 if the keypoints are storted.
 **/

vl_bool
check_sorted (double const * keys, vl_size nkeys)
{
  vl_uindex k ;
  for (k = 0 ; k + 1 < nkeys ; ++ k) {
    if (korder(keys, keys + 4) > 0) {
      return VL_FALSE ;
    }
    keys += 4 ;
  }
  return VL_TRUE ;
}


/** ------------------------------------------------------------------
 ** @brief MEX entry point
 **/

void
mexFunction(int nout, mxArray *out[],
            int nin, const mxArray *in[])
{
  enum {IN_I=0,IN_END} ;
  enum {OUT_FRAMES=0, OUT_DESCRIPTORS, OUT_GSS, OUT_RESP, OUT_GRAD} ;

  int                verbose = 0 ;
  int                opt ;
  int                next = IN_END ;
  mxArray const     *optarg ;

  float const *data ;
  int                M, N ;

  int                O     = - 1 ;
  int                S     =   3 ;
  int                o_min =   0 ;

  double             edge_thresh = -1 ;
  double             peak_thresh = -1 ;
  double             norm_thresh = -1 ;
  double             magnif      = -1 ;
  double             window_size = -1 ;

  mxArray           *ikeys_array = 0 ;
  double            *ikeys = 0 ;
  int                nikeys = -1 ;
  vl_bool            force_orientations = 0 ;
  vl_bool            floatDescriptors = 0 ;

  VL_USE_MATLAB_ENV ;

  /* -----------------------------------------------------------------
   *                                               Check the arguments
   * -------------------------------------------------------------- */

  if (nin < 1) {
    mexErrMsgTxt("One argument required.") ;
  } else if (nout > 5) {
    mexErrMsgTxt("Too many output arguments.");
  }

  if (mxGetNumberOfDimensions (in[IN_I]) != 2              ||
      mxGetClassID            (in[IN_I]) != mxSINGLE_CLASS  ) {
    mexErrMsgTxt("I must be a matrix of class SINGLE") ;
  }

  data = (float*) mxGetData (in[IN_I]) ;
  M    = mxGetM (in[IN_I]) ;
  N    = mxGetN (in[IN_I]) ;

  while ((opt = vlmxNextOption (in, nin, options, &next, &optarg)) >= 0) {
    switch (opt) {

    case opt_verbose :
      ++ verbose ;
      break ;

    case opt_octaves :
      if (!vlmxIsPlainScalar(optarg) || (O = (int) *mxGetPr(optarg)) < 0) {
        mexErrMsgTxt("'Octaves' must be a positive integer.") ;
      }
      break ;

    case opt_levels :
      if (! vlmxIsPlainScalar(optarg) || (S = (int) *mxGetPr(optarg)) < 1) {
        mexErrMsgTxt("'Levels' must be a positive integer.") ;
      }
      break ;

    case opt_first_octave :
      if (!vlmxIsPlainScalar(optarg)) {
        mexErrMsgTxt("'FirstOctave' must be an integer") ;
      }
      o_min = (int) *mxGetPr(optarg) ;
      break ;

    case opt_edge_thresh :
      if (!vlmxIsPlainScalar(optarg) || (edge_thresh = *mxGetPr(optarg)) < 1) {
        mexErrMsgTxt("'EdgeThresh' must be not smaller than 1.") ;
      }
      break ;

    case opt_peak_thresh :
      if (!vlmxIsPlainScalar(optarg) || (peak_thresh = *mxGetPr(optarg)) < 0) {
        mexErrMsgTxt("'PeakThresh' must be a non-negative real.") ;
      }
      break ;

    case opt_norm_thresh :
      if (!vlmxIsPlainScalar(optarg) || (norm_thresh = *mxGetPr(optarg)) < 0) {
        mexErrMsgTxt("'NormThresh' must be a non-negative real.") ;
      }
      break ;

    case opt_magnif :
      if (!vlmxIsPlainScalar(optarg) || (magnif = *mxGetPr(optarg)) < 0) {
        mexErrMsgTxt("'Magnif' must be a non-negative real.") ;
      }
      break ;

    case opt_window_size :
      if (!vlmxIsPlainScalar(optarg) || (window_size = *mxGetPr(optarg)) < 0) {
        mexErrMsgTxt("'WindowSize' must be a non-negative real.") ;
      }
      break ;

    case opt_frames :
      if (!vlmxIsMatrix(optarg, 4, -1)) {
        mexErrMsgTxt("'Frames' must be a 4 x N matrix.x") ;
      }
      ikeys_array = mxDuplicateArray (optarg) ;
      nikeys      = mxGetN (optarg) ;
      ikeys       = mxGetPr (ikeys_array) ;
      if (! check_sorted (ikeys, nikeys)) {
        qsort (ikeys, nikeys, 4 * sizeof(double), korder) ;
      }
      break ;

    case opt_orientations :
      force_orientations = 1 ;
      break ;

    case opt_float_descriptors :
      floatDescriptors = 1 ;
      break ;

    default :
      abort() ;
    }
  }

  /* -----------------------------------------------------------------
   *                                                            Do job
   * -------------------------------------------------------------- */
  {
    VlCovDet          *det ;
    VlScaleSpace        *gss;
    VlScaleSpace        *grad;
    VlScaleSpace const  *resp;
    double              *frames = 0 ;
    void                *descr  = 0 ;
    float const         *descrs = 0;
    int                  nframes = 0, reserved = 0, i,j,q ;
    vl_bool              calc_descs = nout > 1;
    vl_bool              calc_orient = VL_TRUE;
    VlFrameOrientedDisc const *keys  = 0 ;
    int nkeys = 0 ;

    /* create a filter to process the image */
    det = vl_covdet_new_disc_detector (M, N, VL_COVDET_METHOD_DOG, O, o_min, S, calc_orient, calc_descs) ;
    gss = vl_covdet_get_gss(det);
    grad = vl_covdet_get_gradient(det);
    resp = vl_covdet_get_responses(det);

    if (peak_thresh >= 0) vl_covdet_set_peak_thresh (det, peak_thresh) ;
    if (edge_thresh >= 0) vl_covdet_set_edge_thresh (det, edge_thresh) ;
    if (norm_thresh >= 0) vl_covdet_set_norm_thresh (det, norm_thresh) ;
    if (magnif      >= 0) vl_covdet_set_magnif      (det, magnif) ;
    if (window_size >= 0) vl_covdet_set_window_size (det, window_size) ;

    if (verbose) {
      mexPrintf("vl_sift: filter settings:\n") ;
      mexPrintf("vl_sift:   octaves      (O)      = %d\n",
                vl_scalespace_get_octaves_num   (gss)) ;
      mexPrintf("vl_sift:   levels       (S)      = %d\n",
                vl_scalespace_get_levels_num    (gss)) ;
      mexPrintf("vl_sift:   first octave (o_min)  = %d\n",
                vl_scalespace_get_octave_min    (gss)) ;
      mexPrintf("vl_sift:   edge thresh           = %g\n",
                vl_covdet_get_edge_thresh      (det)) ;
      mexPrintf("vl_sift:   peak thresh           = %g\n",
                vl_covdet_get_peak_thresh      (det)) ;
      mexPrintf("vl_sift:   norm thresh           = %g\n",
                vl_covdet_get_norm_thresh     (det)) ;
      mexPrintf("vl_sift:   magnif                = %g\n",
                vl_covdet_get_magnif          (det)) ;
      mexPrintf("vl_sift:   window size           = %g\n",
                vl_covdet_get_window_size     (det)) ;
      mexPrintf("vl_sift:   float descriptor      = %d\n",
                floatDescriptors) ;

      mexPrintf((nikeys >= 0) ?
                "vl_sift: will source frames? yes (%d read)\n" :
                "vl_sift: will source frames? no\n", nikeys) ;
      mexPrintf("vl_sift: will force orientations? %s\n",
                force_orientations ? "yes" : "no") ;
    }

    /* ...............................................................
     *                                             Process each octave
     * ............................................................ */

    if (nikeys < 0) {
      vl_covdet_detect(det, data);

      keys  = vl_covdet_get_oriented_discs (det) ;
      nkeys = vl_covdet_get_frames_num  (det) ;
      descrs = vl_covdet_get_descriptors   (det);

      if (verbose > 1) {
        printf ("vl_sift: detected %d (oriented) keypoints\n", nkeys) ;
      }
    } else {
      vl_scalespace_init(gss, data);
      vl_scalespace_apply(gss, grad, vl_imgradient_polar_f_callback, 0);
      nkeys = nikeys ;
    }

    /* For each keypoint ........................................ */
    for (i = 0; i < nkeys ; ++i) {
      VlFrameOrientedDisc const *k = 0 ;
      VlFrameOrientedDisc        tk ;
      VlScaleSpaceFrame    ik ;
      double angles[4];
      int nangles = 0;

      /* Obtain keypoint orientations ........................... */
      if (nikeys >= 0) {
        vl_scalespace_frame_init (resp, &ik,
                               ikeys [4 * i + 1] - 1.,
                               ikeys [4 * i + 0] - 1.,
                               ikeys [4 * i + 2]) ;

        tk.x = ik.x;
        tk.y = ik.y;
        tk.sigma = ik.sigma;
        k = &tk;

        /* optionally compute orientations too */
        if (force_orientations) {
          nangles = vl_covdet_calc_ssframe_orientations
                    (grad, angles, &ik) ;
        } else {
          angles [0] = VL_PI / 2 - ikeys [4 * i + 3] ;
          nangles    = 1 ;
        }
      } else {
        k = keys + i;
        angles[0] = k->angle;
        nangles = 1;
      }

      /* For each orientation ................................... */
      for (q = 0 ; q < nangles ; ++q) {
        float  rbuf [128] ; /* Buffer for transposed descriptor */

        /* make enough room for all these keypoints and more */
        if (reserved < nframes + 1) {
          reserved += 2 * nkeys ;
          frames = mxRealloc (frames, 4 * sizeof(double) * reserved) ;
          if (nout > 1) {
            if (! floatDescriptors) {
              descr  = mxRealloc (descr,  128 * sizeof(vl_uint8) * reserved) ;
            } else {
              descr  = mxRealloc (descr,  128 * sizeof(float) * reserved) ;
            }
          }
        }

        /* compute descriptor (if necessary) */
        if (nout > 1) {
          if (nikeys >= 0){
            float buf [128] ;
            vl_covdet_calc_ssframe_descriptor (det, buf, &ik, angles [q]) ;
            transpose_descriptor (rbuf, buf) ;
          } else {
            transpose_descriptor (rbuf, descrs + i * 128) ;
          }
        }

        /* Save back with MATLAB conventions. Notice tha the input
         * image was the transpose of the actual image. */
        frames [4 * nframes + 0] = k -> y + 1. ;
        frames [4 * nframes + 1] = k -> x + 1. ;
        frames [4 * nframes + 2] = k -> sigma ;
        frames [4 * nframes + 3] = VL_PI / 2 - angles [q] ;

        if (nout > 1) {
          if (! floatDescriptors) {
            for (j = 0 ; j < 128 ; ++j) {
              float x = 512.0F * rbuf [j] ;
              x = (x < 255.0F) ? x : 255.0F ;
              ((vl_uint8*)descr) [128 * nframes + j] = (vl_uint8) x ;
            }
          } else {
            for (j = 0 ; j < 128 ; ++j) {
              float x = 512.0F * rbuf [j] ;
              ((float*)descr) [128 * nframes + j] = x ;
            }
          }
        }

        ++ nframes ;
      } /* next orientation */
    } /* next keypoint */

    if (verbose) {
      mexPrintf ("vl_sift: found %d keypoints\n", nframes) ;
    }

    /* ...............................................................
     *                                                       Save back
     * ............................................................ */

    {
      mwSize dims [2] ;

      /* create an empty array */
      dims [0] = 0 ;
      dims [1] = 0 ;
      out[OUT_FRAMES] = mxCreateNumericArray
        (2, dims, mxDOUBLE_CLASS, mxREAL) ;

      /* set array content to be the frames buffer */
      dims [0] = 4 ;
      dims [1] = nframes ;
      mxSetPr         (out[OUT_FRAMES], frames) ;
      mxSetDimensions (out[OUT_FRAMES], dims, 2) ;

      if (nout > 1) {

        /* create an empty array */
        dims [0] = 0 ;
        dims [1] = 0 ;
        out[OUT_DESCRIPTORS]= mxCreateNumericArray
          (2, dims,
           floatDescriptors ? mxSINGLE_CLASS : mxUINT8_CLASS,
           mxREAL) ;

        /* set array content to be the descriptors buffer */
        dims [0] = 128 ;
        dims [1] = nframes ;
        mxSetData       (out[OUT_DESCRIPTORS], descr) ;
        mxSetDimensions (out[OUT_DESCRIPTORS], dims, 2) ;
        if(nout > 2) {
          out[OUT_GSS] = vlmxCreateArrayFromScaleSpace(gss);
          if (nout > 3) {
            out[OUT_RESP] = vlmxCreateArrayFromScaleSpace(resp);
            if (nout > 4) {
              out[OUT_GRAD] = vlmxCreateArrayFromScaleSpace(grad);
            }
          }
        }
      }
    }

    /* cleanup */
    vl_covdet_delete (det) ;

    if (ikeys_array)
      mxDestroyArray(ikeys_array) ;

  } /* end: do job */
}

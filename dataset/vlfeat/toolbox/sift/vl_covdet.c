/** @internal @file covdet.c
 ** @author Andrea Vedaldi
 ** @brief Scale Invariant Feature Transform (SIFT) - MEX
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
  opt_method = 0,
  opt_affineadaptation,
  opt_orientation,
  opt_octaves,
  opt_levels,
  opt_first_octave,
  opt_frames,
  opt_edge_thresh,
  opt_peak_thresh,
  opt_norm_thresh,
  opt_magnif,
  opt_window_size,
  opt_aff_win_size,
  opt_aff_max_iter,
  opt_aff_conv_thr,
  opt_patch_size,
  opt_float_descriptors,
  opt_verbose
} ;

/* options */
vlmxOption  options [] = {
  {"Method",           1,   opt_method            },
  {"AffineAdaptation", 1,   opt_affineadaptation  },
  {"Orientation",      1,   opt_orientation       },
  {"NumOctaves",       1,   opt_octaves           },
  {"NumLevels",        1,   opt_levels            },
  {"FirstOctave",      1,   opt_first_octave      },
  {"Frames",           1,   opt_frames            },
  {"PeakThreshold",    1,   opt_peak_thresh       },
  {"EdgeThreshold",    1,   opt_edge_thresh       },
  {"NormThreshold",    1,   opt_norm_thresh       },
  {"Magnif",           1,   opt_magnif            },
  {"WindowSize",       1,   opt_window_size       },
  {"AffineWindowSize", 1,   opt_aff_win_size      },
  {"AffineMaxNumIter", 1,   opt_aff_max_iter      },
  {"AffConvThr",       1,   opt_aff_conv_thr      },
  {"SIFTPatchSize",    1,   opt_patch_size        },
  {"FloatDescriptors", 0,   opt_float_descriptors },
  {"Verbose",          0,   opt_verbose           },
  {0,                  0,   0                     }
} ;

/** ------------------------------------------------------------------
 ** @internal @brief Transpose desriptor
 ** @param dst destination buffer.
 ** @param src source buffer.
 **
 ** The function writes to @a dst the transpose of the SIFT descriptor
 ** @a src. The tranpsose is defined as the descriptor that one
 ** obtains from computing the normal descriptor on the transposed
 ** image.
 **/

static void
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

/** @internal @brief Ordering of tuples by increasing scale
 ** @param a tuple.
 ** @param b tuble.
 ** @return @c a[2] < b[2]
 **/

static int
korder (void const* a, void const* b) {
  double x = ((double*) a) [2] - ((double*) b) [2] ;
  if (x < 0) return -1 ;
  if (x > 0) return +1 ;
  return 0 ;
}

/** @internal @brief Check for sorted keypoints
 ** @param keys keypoint list to check
 ** @param nkeys size of the list.
 ** @return 1 if the keypoints are storted.
 **/

static vl_bool
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
  enum {IN_I = 0, IN_END} ;
  enum {OUT_FRAMES=0, OUT_DESCRIPTORS, OUT_GSS, OUT_RESP, OUT_END} ;

  int verbose = 0 ;
  int opt ;
  int next = IN_END ;
  mxArray const *optarg ;
  VlEnumerator *pair ;

  float const *image ;
  vl_size width, height ;

  VlCovDetMethod method = VL_COVDET_METHOD_DOG;
  vl_bool doAffineAdaptation = VL_FALSE;
  vl_bool doOrientation = VL_TRUE;
  VlFrameType frameType ;
  vl_size frameSize ;
  vl_size frameNumValues ;

  int numOctaves = - 1 ;
  int numLevels = 3 ;
  int firstOctave = 0 ;
  int patchSize = 41 ;

  double edgeThreshold = -1 ;
  double peakThreshold = -1 ;
  double normThreshold = -1 ;
  double magnif = -1 ;
  double window_size = -1 ;

  int affineWindowSize = 19 ;
  int affineAdaptationMaxNumIterations = -1 ;
  double aff_conv_thr  = -1 ;

  void *inputFrames = NULL ;
  vl_index numInputFrames = -1 ; /* signed */
  VlFrameType inputFramesType = 0;
  vl_size inputFramesSize;
  vl_bool floatDescriptors = VL_FALSE ;

  VL_USE_MATLAB_ENV ;

  /* -----------------------------------------------------------------
   *                                               Check the arguments
   * -------------------------------------------------------------- */

  if (nin < IN_END) {
    vlmxError(vlmxErrNotEnoughInputArguments, 0) ;
  } else if (nout > OUT_END) {
    vlmxError(vlmxErrTooManyOutputArguments, 0) ;
  }

  if (mxGetNumberOfDimensions(IN(I)) != 2 || mxGetClassID(IN(I)) != mxSINGLE_CLASS) {
    vlmxError(vlmxErrInvalidArgument, "I must be a matrix of class SINGLE.") ;
  }

  image = (float*) mxGetData(IN(I)) ;
  height = mxGetM(IN(I)) ;
  width = mxGetN(IN(I)) ;

  while ((opt = vlmxNextOption (in, nin, options, &next, &optarg)) >= 0) {
    switch (opt) {

    case opt_verbose :
      ++ verbose ;
      break ;

    case opt_method:
      pair = vlmxDecodeEnumeration(optarg, vlCovdetMethods, VL_TRUE) ;
      if (pair == NULL) {
        vlmxError(vlmxErrInvalidArgument, "METHOD is not a supported detection method.") ;
      }
      method = (VlCovDetMethod)pair->value ;
      break;

    case opt_affineadaptation:
      if (!mxIsLogicalScalar(optarg)) {
        vlmxError(vlmxErrInvalidArgument, "AFFINEADAPTATION must be a logical scalar value.") ;
      } else {
        doAffineAdaptation = *mxGetLogicals(optarg) ;
      }
      break ;

    case opt_orientation:
      if (!mxIsLogicalScalar(optarg)) {
        vlmxError(vlmxErrInvalidArgument, "ORIENTATION must be a logical scalar value.") ;
      } else {
        doOrientation = *mxGetLogicals(optarg);
      }
      break ;

    case opt_octaves :
      if (!vlmxIsPlainScalar(optarg) || (numOctaves = (int) *mxGetPr(optarg)) < 0) {
        vlmxError(vlmxErrInvalidArgument, "OCTAVES must be a positive integer.") ;
      }
      break ;

    case opt_levels :
      if (!vlmxIsPlainScalar(optarg) || (numLevels = (int) *mxGetPr(optarg)) < 1) {
        vlmxError(vlmxErrInvalidArgument, "LEVELS must be a positive integer.") ;
      }
      break ;

    case opt_first_octave :
      if (!vlmxIsPlainScalar(optarg)) {
        vlmxError(vlmxErrInvalidArgument, "FIRSTOCTAVE must be an integer") ;
      }
      firstOctave = (int) *mxGetPr(optarg) ;
      break ;

    case opt_edge_thresh :
      if (!vlmxIsPlainScalar(optarg) || (edgeThreshold = *mxGetPr(optarg)) < 1) {
        vlmxError(vlmxErrInvalidArgument, "EDGETHRESHOLD must be a real not smaller than 1.") ;
      }
      break ;

    case opt_peak_thresh :
      if (!vlmxIsPlainScalar(optarg) || (peakThreshold = *mxGetPr(optarg)) < 0) {
        vlmxError(vlmxErrInvalidArgument, "PEAKTHRESHOLD must be a non-negative real.") ;
      }
      break ;

    case opt_norm_thresh :
      if (!vlmxIsPlainScalar(optarg) || (normThreshold = *mxGetPr(optarg)) < 0) {
        vlmxError(vlmxErrInvalidArgument, "NORMTHRESHOLD must be a non-negative real.") ;
      }
      break ;

    case opt_magnif :
      if (!vlmxIsPlainScalar(optarg) || (magnif = *mxGetPr(optarg)) < 0) {
        vlmxError(vlmxErrInvalidArgument, "MAGNIF must be a non-negative real.") ;
      }
      break ;

    case opt_window_size :
      if (!vlmxIsPlainScalar(optarg) || (window_size = *mxGetPr(optarg)) < 0) {
        vlmxError(vlmxErrInvalidArgument, "WINDOWSIZE must be a non-negative real.") ;
      }
      break ;

    case opt_aff_win_size :
      if (! vlmxIsPlainScalar(optarg) ||
          (affineWindowSize = (int) *mxGetPr(optarg)) < 1) {
        vlmxError(vlmxErrInvalidArgument, "AFFWINSIZE must be a positive integer.") ;
      }
      break ;

    case opt_aff_max_iter :
      if (! vlmxIsPlainScalar(optarg) ||
          (affineAdaptationMaxNumIterations = (int) *mxGetPr(optarg)) < 1) {
        vlmxError(vlmxErrInvalidArgument, "AFFMAXITER must be a positive integer.") ;
      }
      break ;

    case opt_aff_conv_thr :
      if (!vlmxIsPlainScalar(optarg) ||
          (aff_conv_thr = *mxGetPr(optarg)) < 0) {
        vlmxError(vlmxErrInvalidArgument, "AFFCONVTHR must be a non-negative real.") ;
      }
      break ;

    case opt_patch_size :
      if (! vlmxIsPlainScalar(optarg)
          || (affineAdaptationMaxNumIterations = (int) *mxGetPr(optarg)) < 1) {
        vlmxError(vlmxErrInvalidArgument, "SIFTPATCHSIZE must be a positive integer.") ;
      }
      break ;

    case opt_frames : {
        int mifrms, i;
        double const *mx_frms_ptr;
        if (!vlmxIsMatrix(optarg, -1, -1)) {
          vlmxError(vlmxErrInvalidArgument, "FRAMES must be a plain matrix.") ;
        }

        numInputFrames       = mxGetN  (optarg) ;
        mifrms       = mxGetM  (optarg) ;
        mx_frms_ptr = mxGetPr (optarg) ;
        switch (mifrms) {
        case 3: {
          VlFrameDisc *frm;
          inputFramesType = VL_FRAMETYPE_DISC;
          inputFramesSize = vl_frame_size(inputFramesType);
          inputFrames = vl_malloc(inputFramesSize * numInputFrames);
          frm = (VlFrameDisc*)inputFrames;
          for (i = 0; i < numInputFrames; ++i) {
            int ifrms_i = i * mifrms;
            frm->x = mx_frms_ptr[ifrms_i + 1] - 1.;
            frm->y = mx_frms_ptr[ifrms_i + 0] - 1.;
            frm->sigma = mx_frms_ptr[ifrms_i + 2];
            frm++;
          }
        } break;
        case 4: {
          VlFrameOrientedDisc *frm;
          inputFramesType = VL_FRAMETYPE_ORIENTED_DISC;
          inputFramesSize = vl_frame_size(inputFramesType);
          inputFrames = vl_malloc(inputFramesSize * numInputFrames);
          frm = (VlFrameOrientedDisc*)inputFrames;
          for (i = 0; i < numInputFrames; ++i) {
            int ifrms_i = i * mifrms;
            frm->x = mx_frms_ptr[ifrms_i + 1] - 1.;
            frm->y = mx_frms_ptr[ifrms_i + 0] - 1.;
            frm->sigma = mx_frms_ptr[ifrms_i + 2];
            frm->angle = VL_PI / 2. - mx_frms_ptr[ifrms_i + 3];
            frm++;
          }
        } break;
        case 5: {
          VlFrameEllipse *frm;
          inputFramesType = VL_FRAMETYPE_ELLIPSE;
          inputFramesSize = vl_frame_size(inputFramesType);
          inputFrames = vl_malloc(inputFramesSize * numInputFrames);
          frm = (VlFrameEllipse*)inputFrames;
          for (i = 0; i < numInputFrames; ++i) {
            int ifrms_i = i * mifrms;
            frm->x = mx_frms_ptr[ifrms_i + 1] - 1.;
            frm->y = mx_frms_ptr[ifrms_i + 0] - 1.;
            frm->e11 = mx_frms_ptr[ifrms_i + 4];
            frm->e12 = mx_frms_ptr[ifrms_i + 3];
            frm->e22 = mx_frms_ptr[ifrms_i + 2];
            frm++;
          }
        } break;
        case 6: {
          VlFrameOrientedEllipse *frm;
          inputFramesType = VL_FRAMETYPE_ORIENTED_ELLIPSE;
          inputFramesSize = vl_frame_size(inputFramesType);
          inputFrames = vl_malloc(inputFramesSize * numInputFrames);
          frm = (VlFrameOrientedEllipse*)inputFrames;
          for (i = 0; i < numInputFrames; ++i) {
            int ifrms_i = i * mifrms;
            frm->x = mx_frms_ptr[ifrms_i + 1] - 1.;
            frm->y = mx_frms_ptr[ifrms_i + 0] - 1.;
            frm->a11 = mx_frms_ptr[ifrms_i + 5];
            frm->a12 = mx_frms_ptr[ifrms_i + 4];
            frm->a21 = mx_frms_ptr[ifrms_i + 3];
            frm->a22 = mx_frms_ptr[ifrms_i + 2];
            frm++;
          }
        } break;
        default:
          vlmxError(vlmxErrInvalidArgument, "FRAMES must have at least 3 rows and at most 6.") ;
          break;
        }
      } break ;

    case opt_float_descriptors :
      floatDescriptors = 1 ;
      break ;

    default :
      abort() ;
    }
  }

  frameType = vl_frame_get_type(doAffineAdaptation, doOrientation);
  frameSize = vl_frame_size(frameType );
  frameNumValues = frameSize  / sizeof(float);

  /* -----------------------------------------------------------------
   *                                                            Do job
   * -------------------------------------------------------------- */
  {
    VlCovDet          *det ;
    VlAffineShapeEstimator         *aff_det = 0 ;
    VlAffinePatchNormalizer        *aff_norm = 0 ;
    VlScaleSpace        *gss;
    VlScaleSpace const  *resp;
    void const          *frames = 0 ;
    float const         *descrs = 0 ;
    double              *mx_frames = 0 ;
    void                *mx_descrs = 0 ;
    int                  nframes = 0, i, j;
    vl_bool              calc_descs = nout > 1;
    vl_size              desc_sz = 0;
    double time ;

    /* create a filter to process the image */
    if (doAffineAdaptation){
      det = vl_covdet_new_ellipse_detector (height, width, method,
                                            numOctaves, firstOctave, numLevels,
                                            affineWindowSize, patchSize,
                                            doOrientation, calc_descs) ;
    } else {
      det = vl_covdet_new_disc_detector (height, width, method,
                                         numOctaves, firstOctave, numLevels,
                                         doOrientation, calc_descs) ;
    }
    gss = vl_covdet_get_gss(det);
    resp = vl_covdet_get_responses(det);
    if (doAffineAdaptation) {
      aff_det = vl_covdet_get_affdet(det);
      aff_norm = vl_covdet_get_affnorm(det);
    }
    if (calc_descs){
      desc_sz = vl_covdet_get_descriptor_size(det);
    }

    if (peakThreshold >= 0) vl_covdet_set_peak_thresh (det, peakThreshold) ;
    if (edgeThreshold >= 0) vl_covdet_set_edge_thresh (det, edgeThreshold) ;
    if (calc_descs) {
      if (normThreshold >= 0)  vl_covdet_set_norm_thresh (det, normThreshold);
      if (magnif >= 0) vl_covdet_set_magnif (det, magnif) ;
      if (!doAffineAdaptation) {
        /* Do not set these parameters for ellipse frames detection
         * - they are adjusted by covdet in order to calculate descriptors
         * from patch of constant size.
         */
        if (window_size >= 0)  vl_covdet_set_window_size (det, window_size);
      }
    }
    if (doAffineAdaptation) {
      if (affineAdaptationMaxNumIterations >= 0) {
        vl_affineshapeestimator_set_max_iter (aff_det, affineAdaptationMaxNumIterations) ;
      }
      if (aff_conv_thr >= 0) {
        vl_affineshapeestimator_set_conv_thresh (aff_det, aff_conv_thr);
      }
    }
	
    if (verbose) {
      mexPrintf("vl_covdet: method = %s\n", vl_enumeration_get_by_value(vlCovdetMethods, method)->name) ;
      mexPrintf("vl_covdet: num octaves = %d\n", vl_scalespace_get_octaves_num(gss)) ;
      mexPrintf("vl_covdet: num levels per octave = %d\n", vl_scalespace_get_levels_num(gss)) ;
      mexPrintf("vl_covdet: first octave = %d\n", vl_scalespace_get_octave_min(gss)) ;
      mexPrintf("vl_covdet: edge threshold = %g\n", vl_covdet_get_edge_thresh(det)) ;
      mexPrintf("vl_covdet: peak threshold = %g\n", vl_covdet_get_peak_thresh(det)) ;
      if (calc_descs) {
        mexPrintf("vl_covdet: norm threshold = %g\n", vl_covdet_get_norm_thresh(det)) ;
        mexPrintf("vl_covdet: magnification factor = %g\n", vl_covdet_get_magnif(det)) ;
        if (!doAffineAdaptation) {
          mexPrintf("vl_covdet: window size = %g\n", vl_covdet_get_window_size(det)) ;
        }
      }
      if (doAffineAdaptation) {
        mexPrintf("vl_covdet: affine adaptation is turned on\n") ;
        mexPrintf("vl_covdet: affine window size = %d\n", vl_affineshapeestimator_get_window_size(aff_det)) ;
        mexPrintf("vl_covdet: max num of affine adaptation iterations = %d\n", vl_affineshapeestimator_get_max_iter(aff_det)) ;
        mexPrintf("vl_covdet: conv_thr = %g\n", vl_affineshapeestimator_get_conv_thresh(aff_det)) ;
        mexPrintf("vl_covdet: patch_size = %d\n", vl_affinepatchnormalizer_get_patch_size(aff_norm)) ;
      } else {
        mexPrintf("vl_covdet: affine adaptation is turned off\n") ;
      }
      mexPrintf("vl_covdet: float descriptor = %s\n", VL_YESNO(floatDescriptors)) ;
      mexPrintf((numInputFrames >= 0) ?
                "vl_covdet: will source frames? yes (%d read)\n" :
                "vl_covdet: will source frames? no\n", numInputFrames) ;
    }

    if (numInputFrames < 0) {
      /* Detect frames ........................................ */
      time = vl_get_cpu_time() ;
      vl_covdet_detect(det, image);
      time = vl_get_cpu_time() - time ;

      frames  = vl_covdet_get_frames_storage (det) ;
      nframes = vl_covdet_get_frames_num (det) ;
      descrs = vl_covdet_get_descriptors (det);

      if (verbose) {
        VL_PRINTF("vl_covdet: detected %d frames of type '%s' in %.2f seconds\n",
                  nframes, vlFrameNames[frameType], time) ;
      }
    } else {
      vl_covdet_convert_frames(det, image, inputFrames, numInputFrames, inputFramesType);

      frames  = vl_covdet_get_frames_storage (det) ;
      nframes = vl_covdet_get_frames_num (det) ;
      descrs = vl_covdet_get_descriptors (det);
    }

    /* make enough room for all the keypoints */
    mx_frames = mxMalloc(frameNumValues  * sizeof(double) * nframes) ;
    if (calc_descs) {
      if (! floatDescriptors) {
        mx_descrs  = mxMalloc(desc_sz * sizeof(vl_uint8) * nframes) ;
      } else {
        mx_descrs  = mxMalloc(desc_sz * sizeof(float) * nframes) ;
      }
    }

    /* For each keypoint ........................................ */
    for (i = 0; i < nframes ; ++i) {
      float  rbuf [128] ; /* Buffer for transposed descriptor */
      int frm_i = i * frameNumValues ;
      int desc_i = i * desc_sz;

      /* process descriptor (if necessary) */
      if (nout > 1) {
        transpose_descriptor (rbuf, descrs + desc_i) ;
      }

      /* Save back with MATLAB conventions. Notice tha the input
       * image was the transpose of the actual image. */

      switch (frameType ) {
      case VL_FRAMETYPE_DISC: {
        VlFrameDisc const *frm = (VlFrameDisc const *)frames;
        mx_frames [frm_i + 0] = frm -> y + 1. ;
        mx_frames [frm_i + 1] = frm -> x + 1. ;
        mx_frames [frm_i + 2] = frm -> sigma ;
        frames = (void const *)(frm + 1);
      } break;
      case VL_FRAMETYPE_ORIENTED_DISC:{
        VlFrameOrientedDisc const *frm = (VlFrameOrientedDisc const *)frames;
        mx_frames [frm_i + 0] = frm -> y + 1. ;
        mx_frames [frm_i + 1] = frm -> x + 1. ;
        mx_frames [frm_i + 2] = frm -> sigma ;
        mx_frames [frm_i + 3] = VL_PI / 2. - frm->angle;
        frames = (void const *)(frm + 1);
      } break;
      case VL_FRAMETYPE_ELLIPSE:{
        VlFrameEllipse const *frm = (VlFrameEllipse const *)frames;
        mx_frames [frm_i + 0] = frm -> y + 1. ;
        mx_frames [frm_i + 1] = frm -> x + 1. ;
        mx_frames [frm_i + 2] = frm -> e22 ;
        mx_frames [frm_i + 3] = frm -> e12 ;
        mx_frames [frm_i + 4] = frm -> e11 ;
        frames = (void const *)(frm + 1);
      } break;
      case VL_FRAMETYPE_ORIENTED_ELLIPSE:{
        VlFrameOrientedEllipse const *frm =
            (VlFrameOrientedEllipse const *)frames;
        mx_frames [frm_i + 0] = frm -> y + 1. ;
        mx_frames [frm_i + 1] = frm -> x + 1. ;
        mx_frames [frm_i + 2] = frm -> a22 ;
        mx_frames [frm_i + 3] = frm -> a21 ;
        mx_frames [frm_i + 4] = frm -> a12 ;
        mx_frames [frm_i + 5] = frm -> a11 ;
        frames = (void const *)(frm + 1);
      } break;
      default:
        VL_ASSERT(0,"Invalid frame type.");
        break;
      }

      if (calc_descs) {
        if (! floatDescriptors) {
          for (j = 0 ; j < 128 ; ++j) {
            float x = 512.0F * rbuf [j] ;
            x = (x < 255.0F) ? x : 255.0F ;
            ((vl_uint8*)mx_descrs) [desc_i + j] = (vl_uint8) x ;
          }
        } else {
          for (j = 0 ; j < 128 ; ++j) {
            float x = 512.0F * rbuf [j] ;
            ((float*)mx_descrs) [desc_i + j] = x ;
          }
        }
      }
    } /* next keypoint */

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
      dims [0] = frameNumValues  ;
      dims [1] = nframes ;
      mxSetPr         (out[OUT_FRAMES], mx_frames) ;
      mxSetDimensions (out[OUT_FRAMES], dims, 2) ;

      if (calc_descs) {

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
        mxSetData       (out[OUT_DESCRIPTORS], mx_descrs) ;
        mxSetDimensions (out[OUT_DESCRIPTORS], dims, 2) ;
        if(nout > 2) {
          out[OUT_GSS] = vlmxCreateArrayFromScaleSpace(gss);
          if (nout > 3) {
            out[OUT_RESP] = vlmxCreateArrayFromScaleSpace(resp);
          }
        }
      }
    }

    /* cleanup */
    vl_covdet_delete (det) ;

    if (inputFrames)
      vl_free(inputFrames);

  } /* end: job done */
}

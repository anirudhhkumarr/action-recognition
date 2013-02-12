/** @file covdet.h
 ** @brief Covariant feature detectors (@ref covdet)
 ** @author Karel Lenc
 ** @author Andrea Vedaldi
 **/

/*
Copyright (C) 2007-12 Andrea Vedaldi and Brian Fulkerson.
All rights reserved.

This file is part of the VLFeat library and is made available under
the terms of the BSD license (see the COPYING file).
*/

#ifndef VL_COVDET_H
#define VL_COVDET_H

#include <stdio.h>
#include "generic.h"
#include "stringop.h"
#include "imopv.h"
#include "scalespace.h"

/** Maximum number of angles per unoriented frame */
#define VL_COVDET_MAX_ANGLES 4
/** Bins per orientation */
#define VL_SIFT_NBO 8
/** Bins per spatial position */
#define VL_SIFT_NBP 4

/** ------------------------------------------------------------------
 ** @brief Image response functions
 **
 ** Response functions appliable to the input image used for frames
 ** detection.
 **/

typedef enum _VlCovDetMethod {
  VL_COVDET_METHOD_DOG = 1,          /**< Difference of Gaussians */
  VL_COVDET_METHOD_HESSIAN,          /**< Determinant of hessian matrix */
  VL_COVDET_METHOD_NUM
} VlCovDetMethod;

/** ------------------------------------------------------------------
 ** @brief Mapping between strings and ::VlCovDetMethod values
 **/

VL_EXPORT VlEnumerator vlCovdetMethods [VL_COVDET_METHOD_NUM];

/* ---------------------------------------------------------------- */
/** @name Feature frames
 ** @{ */

/** @brief Types of feature frames */
typedef enum _VlFrameType {
  VL_FRAMETYPE_DISC = 1,            /**< Translation and scale covariant frame */
  VL_FRAMETYPE_ORIENTED_DISC,       /**< Similarity covariant frame */
  VL_FRAMETYPE_ELLIPSE,             /**< Affine covariant frame up to rot. */
  VL_FRAMETYPE_ORIENTED_ELLIPSE,    /**< Affine covariant frame */
  VL_FRAMETYPE_NUM
} VlFrameType ;

/** @brief Names of the frame types */
VL_EXPORT const char* vlFrameNames [VL_FRAMETYPE_NUM] ;

/** @brief Mapping between string values and VlFrameType values */
VL_EXPORT VlEnumerator vlFrameTypes [VL_FRAMETYPE_NUM] ;

/** @brief Disc feature frame */
typedef struct _VlFrameDisc
{
  float x ;     /**< center x-coordinate */
  float y ;     /**< center y-coordinate */
  float sigma ; /**< radius or scale */
} VlFrameDisc ;

/** @brief Oriented disc feature frame
 ** An upright frame has @c angle equal to zero.
 **/
typedef struct _VlFrameOrientedDisc {
  float x ;     /**< center x-coordinate */
  float y ;     /**< center y-coordinate */
  float sigma ; /**< radius or scale */
  float angle ; /**< rotation angle (rad) */
} VlFrameOrientedDisc ;

/** @brief Ellipse feature frame */
typedef struct _VlFrameEllipse {
  float x ;     /**< center x-coordinate */
  float y ;     /**< center y-coordinate */
  float e11 ;   /**< */
  float e12 ;
  float e22 ;
} VlFrameEllipse ;

/** @brief Oriented ellipse feature frame
 ** The affine transformation transforms the ellipse shape into
 ** a circular region. */
typedef struct _VlFrameOrientedEllipse {
  float x ;     /**< center x-coordinate */
  float y ;     /**< center y-coordinate */
  float a11 ;   /**< */
  float a12 ;
  float a21 ;
  float a22 ;
} VlFrameOrientedEllipse;

/** @brief Get the size of a frame structure
 ** @param frameType identifier of the type of frame.
 ** @return size of the corresponding frame structure in bytes.
 **/

VL_INLINE vl_size
vl_frame_size (VlFrameType frameType) {
  switch (frameType) {
  case VL_FRAMETYPE_DISC: return sizeof(VlFrameDisc);
  case VL_FRAMETYPE_ORIENTED_DISC: return sizeof(VlFrameOrientedDisc);
  case VL_FRAMETYPE_ELLIPSE: return sizeof(VlFrameEllipse);
  case VL_FRAMETYPE_ORIENTED_ELLIPSE: return sizeof(VlFrameOrientedEllipse);
  default:
    VL_ASSERT(0, "Invalid Frame type.\n");
    break;
  }
  return 0;
}

VL_INLINE VlFrameType
vl_frame_get_type (vl_bool affineAdaptation, vl_bool orientation)
{
  if (affineAdaptation) {
    if (orientation) {
      return VL_FRAMETYPE_ORIENTED_ELLIPSE;
    } else {
      return VL_FRAMETYPE_ELLIPSE;
    }
  } else {
    if (orientation) {
      return VL_FRAMETYPE_ORIENTED_DISC;
    } else {
      return VL_FRAMETYPE_DISC;
    }
  }
}

/* ---------------------------------------------------------------- */
/*                                                  SIFT Descriptor */
/* ---------------------------------------------------------------- */

VL_EXPORT int
vl_calc_frame_orientations(float* gradient, int imageWidth, int imageHeight,
                           double angles[VL_COVDET_MAX_ANGLES],
                           float frameX, float frameY,
                           float frameSigma, int frameOctave);

VL_EXPORT void
vl_sift_calc_descriptor(float const* gradient,
                        int imageWidth, int imageHeight,
                        float *descriptor,
                        double x, double y, double sigma, double angle0,
                        double magnif, double windowSize,
                        double normThreshold, int borderSize) ;

/* ---------------------------------------------------------------- */
/*                                       Affine shape normalisation */
/* ---------------------------------------------------------------- */

/** @brief Affine shape normalisation
 **
 ** Normalises the affine neighbourhood of a keypoint into isotropic shape.
 **
 ** @sa @ref affine
 **/
typedef struct _VlAffinePatchNormalizer
{
  vl_size width ;                /**< Original image width */
  vl_size height ;               /**< Original image height */
  float* workspace ;             /**< Workspace for patch region scaling */
  vl_size wss ;                  /**< Size of the workspace */
  VlImSmoothFilt* smoothFilter ; /**< smooth filter downscaling */
  float* patch ;                 /**< Patch with normalised neighbourhood */
  vl_size patchSize ;            /**< Size of the square patch */
} VlAffinePatchNormalizer ;


/** @name Create and destroy
 ** @{
 **/
VL_EXPORT VlAffinePatchNormalizer *
vl_affinepatchnormalizer_new (vl_size imageWidth, vl_size imageHeight, vl_size patchSize) ;

VL_EXPORT void
vl_affinepatchnormalizer_delete (VlAffinePatchNormalizer *self) ;
/** @} */

VL_EXPORT int
vl_affinepatchnormalizer_normalise (VlAffinePatchNormalizer *self, float const *image,
                                    VlAffineShapeEstimatorFrame const *frame, double magnif);

VL_INLINE vl_size vl_affinepatchnormalizer_get_patch_size (VlAffinePatchNormalizer const *self);
VL_INLINE float *vl_affinepatchnormalizer_get_patch (VlAffinePatchNormalizer const *self);


/* -------------------------------------------------------------------
 *                                     Inline functions implementation
 * ---------------------------------------------------------------- */

/** @brief Get size of the result patch.
 ** @param self Affine shape normalisation object.
 ** @return Size of the square patch with norm. kpt. neighbourhood.
 **/

VL_INLINE vl_size
vl_affinepatchnormalizer_get_patch_size (VlAffinePatchNormalizer const *self)
{
  return self->patchSize ;
}

/** @brief Get patch with normalized kpt neighbourhood.
 ** @param self Affine shape normalisation object.
 ** @return Pointer to the patch of size ::vl_affinepatchnormalizer_get_patch_size
 **/

VL_INLINE float *
vl_affinepatchnormalizer_get_patch (VlAffinePatchNormalizer const *self)
{
  return self->patch ;
}


/* ---------------------------------------------------------------- */
/*                                        Covariant Frames Detector */
/* ---------------------------------------------------------------- */

/** ------------------------------------------------------------------
 ** @brief Covariant frames detector
 **
 **/

typedef struct _VlCovDet
{
  VlScaleSpace *gss ;       /**< gaussian scale space */
  VlScaleSpace *resp ;      /**< responses on the original scale space */
  VlScaleSpace *grad ;      /**< gradient of the original image scale space */

  float * patch ;           /**< feature patch. */
  float * patchGrad ;       /**< gradient of the feature patch. */
  vl_size patchSize ;       /**< size of the feature patch. */
  double patchSigma ;       /**< smoothing of the feature patch. */

  VlAffineShapeEstimator *affineShapeEstimator;   /**< affine shape detector */
  VlAffinePatchNormalizer *affineNorm; /**< affine shape normalisator */

  VlFrameType framesType;   /**< type of the computed frames */
  VlCovDetMethod method;    /**< type of the corner measure for frame loc. */
  vl_bool calcDescriptors;  /**< calc descriptors */
  vl_bool calcOrientation;  /**< calc frame orientations */

  double peakThreshold ;    /**< peak threshold. */
  double edgeThreshold ;    /**< edge threshold. */

  void *frames ;            /**< detected frames. */
  vl_size numFrames ;       /**< number of detected frames. */
  vl_size framesRes ;       /**< Size of the alloc. memory for frames */
  vl_size frameSize ;       /**< Size of the frame in bytes */

  float *descriptors;       /**< pointer to array of descriptors */
  vl_size descriptorsRes;

  double normThreshold ;    /**< norm threshold. */
  double magnification ;    /**< magnification factor. */
  double windowSize ;       /**< size of Gaussian window (in spatial bins) */

  vl_bool calcInverseSm;    /**< Output shape matrix as its inverse (ell. m.)*/

} VlCovDet;


/** @name Create and destroy
 ** @{
 **/
VL_EXPORT VlCovDet *
vl_covdet_new_disc_detector(vl_size imageWidth, vl_size imageHeight,
                            VlCovDetMethod respFunction,
                            vl_size numOctaves, vl_index firstOctave,
                            vl_size numLevels,
                            vl_bool calcOrient, vl_bool calcDescr);

VL_EXPORT VlCovDet *
vl_covdet_new_ellipse_detector(vl_size imageWidth, vl_size imageHeight,
                               VlCovDetMethod imageMeasure,
                               vl_size numOctaves, vl_index firstOctave,
                               vl_size numLevels,
                               vl_size affineWinSize, vl_size descrWinSize,
                               vl_bool calcOrient, vl_bool calcDesc);

VL_EXPORT void
vl_covdet_delete(VlCovDet *self);

/** @} */

/** @name Process data
 ** @{
 **/
VL_EXPORT void
vl_covdet_detect (VlCovDet *self, float const *image) ;

VL_EXPORT void
vl_covdet_convert_frames (VlCovDet *self,
                          float const *image,
                          const void *srcFrames,
                          vl_size srcFramesNum,
                          VlFrameType srcFramesType);

VL_EXPORT void
vl_covdet_calc_ssframe_descriptor (VlCovDet *self, float *descriptor,
                                   VlScaleSpaceFrame const* frame,
                                   double angle0);

VL_EXPORT int
vl_covdet_calc_ssframe_orientations (VlScaleSpace *grad,
                                       double angles [VL_COVDET_MAX_ANGLES],
                                       VlScaleSpaceFrame const *k);

/** @} */

/** @name Retrieve data and parameters
 ** @{
 **/
VL_INLINE vl_size vl_covdet_get_frames_num(VlCovDet const *self);
VL_INLINE vl_size vl_covdet_get_descriptors_num(VlCovDet const *self);

VL_INLINE void const *vl_covdet_get_frames_storage(VlCovDet const *self);
VL_INLINE VlFrameType vl_covdet_get_frames_type(VlCovDet const *self);

VL_INLINE VlFrameDisc const *vl_covdet_get_discs(VlCovDet const *self);
VL_INLINE VlFrameOrientedDisc const *vl_covdet_get_oriented_discs(VlCovDet const *self);
VL_INLINE VlFrameEllipse const *vl_covdet_get_ellipses(VlCovDet const *self);
VL_INLINE VlFrameOrientedEllipse const *vl_covdet_get_oriented_ellipses(VlCovDet const *self);

VL_INLINE vl_size vl_covdet_get_descriptor_size(VlCovDet const *self);
VL_INLINE float const *vl_covdet_get_descriptors(VlCovDet const *self);

VL_INLINE VlScaleSpace *vl_covdet_get_gss(VlCovDet const *self) ;
VL_INLINE VlScaleSpace const *vl_covdet_get_responses (VlCovDet const *self) ;
VL_INLINE VlScaleSpace *vl_covdet_get_gradient(VlCovDet const *self) ;

VL_INLINE VlAffineShapeEstimator* vl_covdet_get_affdet(VlCovDet *self);
VL_INLINE VlAffinePatchNormalizer* vl_covdet_get_affnorm(VlCovDet *self);

VL_INLINE double vl_covdet_get_peak_thresh (VlCovDet const *self) ;
VL_INLINE double vl_covdet_get_edge_thresh (VlCovDet const *self) ;

VL_INLINE double vl_covdet_get_norm_thresh (VlCovDet const *self) ;
VL_INLINE double vl_covdet_get_magnif (VlCovDet const *self) ;
VL_INLINE double vl_covdet_get_window_size (VlCovDet const *self) ;

VL_INLINE vl_bool vl_covdet_get_calc_inverse_sm(VlCovDet const *self) ;

/** @} */

/** @name Set parameters
 ** @{
 **/
VL_INLINE void vl_covdet_set_peak_thresh (VlCovDet *self, double t) ;
VL_INLINE void vl_covdet_set_edge_thresh (VlCovDet *self, double t) ;
VL_INLINE void vl_covdet_set_norm_thresh (VlCovDet *self, double t) ;
VL_INLINE void vl_covdet_set_magnif      (VlCovDet *self, double m) ;
VL_INLINE void vl_covdet_set_window_size (VlCovDet *self, double m) ;
VL_INLINE void vl_covdet_set_calc_inverse_sm (VlCovDet *self, vl_bool v) ;
/** @} */


/* -------------------------------------------------------------------
 *                                     Inline functions implementation
 * ---------------------------------------------------------------- */

/** ------------------------------------------------------------------
 ** @brief Get number of frames.
 ** @param self Covariant frames detector.
 ** @return number of frames.
 **/

VL_INLINE vl_size
vl_covdet_get_frames_num(VlCovDet const *self)
{
  return self-> numFrames;
}

/** ------------------------------------------------------------------
 ** @brief Get number of frame descriptors.
 ** @param self Covariant frames detector.
 ** @return number of descriptors.
 **/

VL_INLINE vl_size
vl_covdet_get_descriptors_num (VlCovDet const *self)
{
  return self->calcDescriptors ? self-> numFrames : 0;
}

/** ------------------------------------------------------------------
 ** @brief Get pointer to the frames storage.
 **
 ** @param self Covariant frames detector.
 ** @return pointer to the array of frames.
 **/

VL_INLINE void const *
vl_covdet_get_frames_storage (VlCovDet const *self)
{
  return (void const *)self->frames;
}

/** ------------------------------------------------------------------
 ** @brief Get output frames type.
 ** @param self Covariant frames detector.
 ** @return type of the detected frames.
 **
 ** Return type of frames for which the detector was configured. The
 ** type of the stored frames can be accessed by
 ** ::vl_covdet_get_frames_type() and the number of detected frames
 ** with ::vl_covdet_get_frames_num().
**
 ** @sa ::vl_covdet_get_frames_type, ::vl_covdet_get_frames_num
 **/

VL_INLINE VlFrameType
vl_covdet_get_frames_type (VlCovDet const *self)
{
  return self->framesType;
}


/** ------------------------------------------------------------------
 ** @brief Get disc frames.
 **
 ** Returns NULL pointer if this type of frames was not detected.
 **
 ** @param self Covariant frames detector.
 ** @return pointer to the list of disc frames.
 ** @sa ::vl_covdet_get_frames_num()
 **/

VL_INLINE VlFrameDisc const *
vl_covdet_get_discs (VlCovDet const *self)
{
  if (self->framesType == VL_FRAMETYPE_DISC) {
    return (VlFrameDisc const *)self->frames;
  } else {
    return NULL;
  }
}


/** ------------------------------------------------------------------
 ** @brief Get oriented disc frames.
 **
 ** Returns NULL pointer if this type of frames was not detected.
 **
 ** @param self Covariant frames detector.
 ** @return pointer to the list of oriented disc frames.
 ** @sa ::vl_covdet_get_frames_num()
 **/

VL_INLINE VlFrameOrientedDisc const *
vl_covdet_get_oriented_discs (VlCovDet const *self)
{
  if (self->framesType == VL_FRAMETYPE_ORIENTED_DISC) {
    return (VlFrameOrientedDisc const *)self->frames;
  } else {
    return NULL;
  }
}

/** ------------------------------------------------------------------
 ** @brief Get ellipse frames.
 **
 ** Return NULL pointer if this type of frames was not detected.
 **
 ** @param self Covariant frames detector.
 ** @return pointer to the list of ellipse frames.
 ** @sa ::vl_covdet_get_frames_num()
 **/

VL_INLINE VlFrameEllipse const *
vl_covdet_get_ellipses (VlCovDet const *self)
{
  if (self->framesType == VL_FRAMETYPE_ELLIPSE) {
    return (VlFrameEllipse const *)self->frames;
  } else {
    return NULL;
  }
}

/** ------------------------------------------------------------------
 ** @brief Get oriented ellipses frames.
 **
 ** Return NULL pointer if this type of frames was not detected.
 **
 ** @param self Covariant frames detector.
 ** @return pointer to the list of oriented ellipse frames.
 ** @sa ::vl_covdet_get_frames_num()
 **/

VL_INLINE VlFrameOrientedEllipse const *
vl_covdet_get_oriented_ellipses (VlCovDet const *self)
{
  if (self->framesType == VL_FRAMETYPE_ORIENTED_ELLIPSE) {
    return (VlFrameOrientedEllipse const *)self->frames;
  } else {
    return NULL;
  }
}

/** ------------------------------------------------------------------
 ** @brief Get size of the SIFT descriptor.
 ** @param self Covariant frames detector.
 ** @return size of the SIFT descriptor.
 **/

VL_INLINE vl_size vl_covdet_get_descriptor_size(VlCovDet const *self)
{
  return self->calcDescriptors ? VL_SIFT_NBO*VL_SIFT_NBP*VL_SIFT_NBP : 0 ;
}

/** ------------------------------------------------------------------
 ** @brief Get descriptors of the frames.
 **
 ** Returns pointer to array of frame descriptors where each
 ** descriptor is of size @a ::vl_covdet_get_descriptor_size.
 **
 ** @param self Covariant frames detector.
 ** @return pointer to array of descriptors or NULL if descriptors not calc.
 **/

VL_INLINE float const *
vl_covdet_get_descriptors (VlCovDet const *self)
{
  return self-> descriptors ;
}

/** ------------------------------------------------------------------
 ** @brief Get the Affine shape detector.
 **
 ** Returns the pointer to affine shape detector used for detecting
 ** anisotropic frames (ellipses).
 **
 ** @param self Covariant frames detector.
 ** @return Affine shape detector or NULL if @a self configured for discs.
 **/

VL_INLINE VlAffineShapeEstimator*
vl_covdet_get_affdet (VlCovDet *self)
{
  return self-> affineShapeEstimator;
}

/** ------------------------------------------------------------------
 ** @brief Get affine shape normalisatior.
 **
 ** Returns pointer to the affine shape normalisation object used
 ** for normalisation of frame anisotropic neighbourhood into an
 ** isotropic neighbourhood.
 **
 ** @param self Covariant frames detector.
 ** @return Affine shape normalisator or NULL if @a self configured for discs.
 **/

VL_INLINE VlAffinePatchNormalizer*
vl_covdet_get_affnorm (VlCovDet *self)
{
  return self-> affineNorm;
}

/** ------------------------------------------------------------------
 ** @brief Get peaks treashold.
 ** @param self Covariant frames detector.
 ** @return Peak threshold.
 **/

VL_INLINE double
vl_covdet_get_peak_thresh (VlCovDet const *self)
{
  return self-> peakThreshold ;
}

/** ------------------------------------------------------------------
 ** @brief Get edges threshold.
 ** @param self Covariant frames detector.
 ** @return Edge threshold.
 **/

VL_INLINE double
vl_covdet_get_edge_thresh (VlCovDet const *self)
{
  return self-> edgeThreshold ;
}

/** ------------------------------------------------------------------
 ** @brief Get gaussian scale space.
 ** @param self Covariant frames detector.
 ** @return Gaussian scale space initialised with the last image.
 **/

VL_INLINE VlScaleSpace *
vl_covdet_get_gss (VlCovDet const *self)
{
  return self-> gss ;
}

/** ------------------------------------------------------------------
 ** @brief Get polar gradient scale space.
 **
 ** Returs pointer to the scale space object of the GSS gradients where
 ** each level is array of size 2*gssLevelWidth x gssLevelHeight.
 ** The first gssLevelWidth x gssLevelHeight layer of the gradient level
 ** contains the gradient magnitude and the second the gradient angle (in
 ** radians, between 0 and @f$ 2\pi @f$) of the scale space level. The
 ** layers are aligned as [gradientMagnitude gradientAngle].
 **
 ** Returns NULL if the detector is configured for anisotropic frames
 ** detection.
 **
 ** @param self Covariant frames detector.
 ** @return Scale space formed from the polar gradient of GSS or NULL.
 **/
VL_INLINE VlScaleSpace *
vl_covdet_get_gradient (VlCovDet const *self)
{
  return self-> grad ;
}

/** ------------------------------------------------------------------
 ** @brief Get responses scale space.
 **
 ** The particular responses values are computed using function
 ** @a VlImRespFunction defined in the covdet object constructor.
 **
 ** In the case of DoG the returned scale space has one level less
 ** in each octave than the GSS. In case of other response functions
 ** the size is the same.
 **
 ** @param self Covariant frames detector.
 ** @return Scale space of responses from the last image.
 **/
VL_INLINE VlScaleSpace const *
vl_covdet_get_responses (VlCovDet const *self)
{
  return self-> resp ;
}

/** ------------------------------------------------------------------
 ** @brief Get norm threshold
 ** @param self Covariant frames detector.
 ** @return Norm threshold.
 **/

VL_INLINE double
vl_covdet_get_norm_thresh (const VlCovDet *self)
{
  return self -> normThreshold ;
}

/** ------------------------------------------------------------------
 ** @brief Get the magnification factor
 ** @param self Covariant frames detector.
 ** @return magnification factor.
 **/
VL_INLINE double
vl_covdet_get_magnif (const VlCovDet *self)
{
  return self -> magnification ;
}

/** ------------------------------------------------------------------
 ** @brief Get the Gaussian window size.
 ** @param self Covariant frames detector.
 ** @return standard deviation of the Gaussian window (in spatial bin units).
 **
 ** This value affects only descriptors calculation of isotropic frames (discs)
 **/
VL_INLINE double
vl_covdet_get_window_size (const VlCovDet *self)
{
  return self -> windowSize ;
}

/** ------------------------------------------------------------------
 ** @brief Find out if the shape matrix is exported as its inverse
 ** @param self Covariant frames detector.
 ** @return VL_TRUE if shape matrix is exported as its inverse.
 **/
VL_INLINE vl_bool
vl_covdet_get_calc_inverse_sm (const VlCovDet *self)
{
  return self -> calcInverseSm ;
}


/** ------------------------------------------------------------------
 ** @brief Set peaks threshold
 ** @param self Covariant frames detector.
 ** @param t Peaks threshold.
 **/
VL_INLINE void
vl_covdet_set_peak_thresh (VlCovDet *self, double t)
{
  self-> peakThreshold = t ;
}

/** ------------------------------------------------------------------
 ** @brief Set edges threshold
 ** @param self Covariant frames detector.
 ** @param t Edges threshold.
 **/
VL_INLINE void
vl_covdet_set_edge_thresh (VlCovDet *self, double t)
{
  self-> edgeThreshold = t ;
}

/** ------------------------------------------------------------------
 ** @brief Set norm threshold
 ** @param self Covariant frames detector.
 ** @param t Norm threshold.
 **/
VL_INLINE void
vl_covdet_set_norm_thresh (VlCovDet *self, double t)
{
  self -> normThreshold = t ;
}

/** ------------------------------------------------------------------
 ** @brief Set the magnification factor
 ** @param self Covariant frames detector..
 ** @param m Magnification factor.
 **/
VL_INLINE void
vl_covdet_set_magnif (VlCovDet *self, double m)
{
  self -> magnification = m ;
}

/** ------------------------------------------------------------------
 ** @brief Set the Gaussian window size
 ** @param self Covariant frames detector.
 ** @param x Gaussian window size (in units of spatial bin).
 **
 ** This is the parameter @f$ \hat \sigma_{\text{win}} @f$ of
 ** the standard SIFT descriptor @ref sift-tech-descriptor-std.
 **
 ** Affects only calculation of descriptors of isotropic frames (discs)
 **/
VL_INLINE void
vl_covdet_set_window_size (VlCovDet *self, double x)
{
  self -> windowSize = x ;
}

/** ------------------------------------------------------------------
 ** @brief Set whether the shape matrix should be exported as its inverse
 ** @param self Covariant frames detector.
 ** @param v VL_TRUE if shape matrix is exported as it inverse
 **
 ** This parameter affects the detector only when ellipse frames are detected.
 **/
VL_INLINE void
vl_covdet_set_calc_inverse_sm (VlCovDet *self, vl_bool v)
{
  /* Set the parameter only for ellipses */
  if (!self->calcOrientation)
    self -> calcInverseSm = v;
}


/* VL_COVDET_H */
#endif

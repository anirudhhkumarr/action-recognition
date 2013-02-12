/** @file scalespace.h
 ** @brief Scale Space (@ref scalespace)
 ** @author Andrea Vedaldi
 ** @author Karel Lenc
 ** @author Michal Perdoch
 **/

/*
Copyright (C) 2007-12 Andrea Vedaldi and Brian Fulkerson.
All rights reserved.

This file is part of the VLFeat library and is made available under
the terms of the BSD license (see the COPYING file).
*/

#ifndef VL_SCALESPACE_H
#define VL_SCALESPACE_H

#include "pgm.h"

#include <stdio.h>
#include <math.h>
#include "generic.h"
#include "imopv.h"

/** ------------------------------------------------------------------
 ** @brief Scale space feature frame
 **
 ** This structure represent a feature frame as extracted by the
 ** ::VlScaleSpace object. It represents a point in space and scale,
 ** including both continuous and approimated integer coordinates
 ** pointing to a discrete sample in the scale space structure.
 **/

typedef struct _VlScaleSpaceFrame
{
  int o ;           /**< o coordinate (octave). */

  int ix ;          /**< Integer unnormalized x coordinate. */
  int iy ;          /**< Integer unnormalized y coordinate. */
  int is ;          /**< Integer s coordinate. */

  float x ;         /**< x coordinate. */
  float y ;         /**< y coordinate. */
  float s ;         /**< s coordinate. */
  float sigma ;     /**< scale. */
} VlScaleSpaceFrame ;

/** ------------------------------------------------------------------
 ** @brief Scale space processing callback
 **
 ** Scale space callback appliable to separate planes.
 **/

typedef void (VlScaleSpaceCallback) (float const *srcImage,
                                     int src_width, int src_height,
                                     float *dstImage,
                                     int dst_width, int dst_height,
                                     int octave, int level,
                                     void *params);

/** ------------------------------------------------------------------
 ** @brief Scale space
 **
 ** Container for Scale space data.
 **/

typedef struct _VlScaleSpace
{
  double sigman ;             /**< nominal image smoothing. */
  double sigma0 ;             /**< smoothing of pyramid base. */
  double sigmak ;             /**< k-smoothing. */
  double dsigma0 ;            /**< delta-smoothing. */

  vl_size width ;             /**< image width. */
  vl_size height ;            /**< image height. */
  vl_size numOctaves ;        /**< number of octaves. */
  vl_size numLevels ;         /**< number of levels per octave. */
  vl_index firstOctave ;      /**< minimum octave index. */
  vl_index lastOctave ;       /**< maximum octave index. */
  vl_index firstLevel ;       /**< minimum level index. */
  vl_index lastLevel ;        /**< maximum level index. */

  float **octaves ;           /**< GSS data. */

  float *patch ;              /** buffer for operations on patches. */
  vl_size patchSize ;         /** size of the @c patch buffer. */

  VlScaleSpaceFrame *frames ; /**< detected feature frames. */
  vl_size framesSize ;        /**< size of the VlScaleSpace::frames buffer. */
  vl_size numFrames ;         /**< number of stored feature frames. */

  VlImSmoothFilt *smoother ;  /**< image smoother. */

} VlScaleSpace ;


/** @name Create and destroy
 ** @{
 **/
VL_EXPORT VlScaleSpace
*vl_scalespace_new (vl_size width, vl_size height,
                    vl_index numOctaves, vl_index firstOctave,
                    vl_size numLevels, vl_index firstLevel, vl_index lastLevel) ;
VL_EXPORT void  vl_scalespace_delete (VlScaleSpace *self) ;
VL_EXPORT VlScaleSpace *vl_scalespace_clone_structure (VlScaleSpace* src);
VL_EXPORT VlScaleSpace *vl_scalespace_clone (VlScaleSpace* src);
VL_EXPORT void vl_scalespace_init(VlScaleSpace *self, float const* image);
/** @} */

/** @name Process data
 ** @{
 **/
VL_EXPORT void
vl_scalespace_apply (VlScaleSpace const *self, VlScaleSpace* dst,
                     VlScaleSpaceCallback* callback, void* params);

VL_EXPORT void
vl_scalespace_diff (VlScaleSpace const* self, VlScaleSpace* dst);

VL_EXPORT void
vl_scalespace_find_local_extrema (VlScaleSpace* self,
                                  double peakThreshold,
                                  vl_size borderSize) ;

VL_EXPORT void
vl_scalespace_refine_local_extrema (VlScaleSpace * self,
                                    double peakThrehsold,
                                    double edgeThreshold,
                                    vl_size borderSize) ;


VL_EXPORT int
vl_scalespace_extract_affine_patch (VlScaleSpace *self,
                                        float *patch,
                                        vl_size patchWidth,
                                        vl_size patchHeight,
                                        double patchSigma,
                                        double t1, double t2,
                                        double a11, double a21,
                                        double a12, double a22) ;

VL_EXPORT void
vl_scalespace_frame_init (VlScaleSpace const * self,
                          VlScaleSpaceFrame * frame,
                          double x, double y, double sigma) ;
/** @} */

/** @name Retrieve data and parameters
 ** @{
 **/
VL_INLINE vl_size vl_scalespace_get_octaves_num (VlScaleSpace const *self) ;
VL_INLINE vl_index vl_scalespace_get_octave_min (VlScaleSpace const *self) ;
VL_INLINE vl_index vl_scalespace_get_octave_max (VlScaleSpace const *self) ;
VL_INLINE vl_size vl_scalespace_get_octave_width (VlScaleSpace const *self, vl_index o) ;
VL_INLINE vl_size vl_scalespace_get_octave_height (VlScaleSpace const *self, vl_index o) ;
VL_INLINE vl_size vl_scalespace_get_levels_num (VlScaleSpace const *self) ;
VL_INLINE vl_index vl_scalespace_get_level_min (VlScaleSpace const *self) ;
VL_INLINE vl_index vl_scalespace_get_level_max (VlScaleSpace const *self) ;
VL_INLINE vl_size vl_scalespace_get_frames_num (VlScaleSpace const *self) ;
VL_INLINE double vl_scalespace_get_sigma0 (VlScaleSpace const *self) ;
VL_INLINE double vl_scalespace_get_sigmak (VlScaleSpace const *self) ;
VL_INLINE float *vl_scalespace_get_octave  (VlScaleSpace const *self, vl_index o, vl_index s) ;
VL_INLINE double vl_scalespace_get_sigma_for_scale (VlScaleSpace const *self, vl_index o, vl_index s) ;
VL_INLINE VlScaleSpaceFrame const *vl_scalespace_get_frames (VlScaleSpace const *self) ;
/** @} */

/* -------------------------------------------------------------------
 *                         ScaleSpace Inline functions implementation
 * ---------------------------------------------------------------- */

/** ------------------------------------------------------------------
 ** @brief Get number of octaves.
 ** @param self ::VlScaleSpace object instance..
 ** @return number of octaves.
 **/

VL_INLINE vl_size
vl_scalespace_get_octaves_num (VlScaleSpace const *self)
{
  assert(self) ;
  return self->numOctaves ;
}

/** -------------------------------------------------------------------
 ** @brief Get index of first octave.
 ** @param self ::VlScaleSpace object instance.
 ** @return index of the first octave.
 **/

VL_INLINE vl_index
vl_scalespace_get_octave_min (VlScaleSpace const *self)
{
  assert(self) ;
  return self->firstOctave ;
}

/** -------------------------------------------------------------------
 ** @brief Get index of last octave.
 ** @param self ::VlScaleSpace object instance.
 ** @return index of the last octave.
 **/

VL_INLINE vl_index
vl_scalespace_get_octave_max (VlScaleSpace const *self)
{
  assert(self) ;
  return self->lastOctave ;
}

/** ------------------------------------------------------------------
 ** @brief Get octave width
 ** @param self ::VlScaleSpace object instance.
 ** @param o octave index.
 ** @return current octave width.
 **/

VL_INLINE vl_size
vl_scalespace_get_octave_width (VlScaleSpace const *self, vl_index o)
{
  assert(self) ;
  return VL_SHIFT_LEFT(self->width, -o) ;
}

/** ------------------------------------------------------------------
 ** @brief Get octave height
 ** @param self ::VlScaleSpace object instance.
 ** @param o octave index.
 ** @return current octave height.
 **/

VL_INLINE vl_size
vl_scalespace_get_octave_height (VlScaleSpace const *self, vl_index o)
{
  assert(self) ;
  return VL_SHIFT_LEFT(self->height, -o) ;
}

/** ------------------------------------------------------------------
 ** @brief Get the data of a scale level
 ** @param self ::VlScaleSpace object instance.
 ** @param o octave index.
 ** @param s level index.
 ** @return pointer to the data for octave @a o, level @a s.
 **
 ** The octave index @a o must be in the range @c firstOctave
 ** to @c lastOctave and the scale index @a s must be in the
 ** range @c firstLevel to @c lastLevel.
 **/

VL_INLINE float *
vl_scalespace_get_octave (VlScaleSpace const *self, vl_index o, vl_index s)
{
  vl_size width = vl_scalespace_get_octave_width(self, o) ;
  vl_size height = vl_scalespace_get_octave_height(self, o) ;
  float *octave = self->octaves[o - self->firstOctave] ;
  assert(self) ;
  assert(o >= self->firstOctave) ;
  assert(o <= self->lastOctave) ;
  assert(s >= self->firstLevel) ;
  assert(s <= self->lastLevel) ;
  return octave + width * height * (s - self->firstLevel) ;
}

/** ------------------------------------------------------------------
 ** @brief Get number of levels per octave
 ** @param self ::VlScaleSpace object instance.
 ** @return number of leves per octave.
 **/

VL_INLINE vl_size
vl_scalespace_get_levels_num (VlScaleSpace const *self)
{
  assert(self) ;
  return self->numLevels ;
}

/** ------------------------------------------------------------------
 ** @brief Get the index of the lowest level
 ** @param self ::VlScaleSpace object instance.
 ** @return index of the lowest level.
 **/

VL_INLINE vl_index
vl_scalespace_get_level_min (VlScaleSpace const *self)
{
  assert(self) ;
  return self->firstLevel ;
}

/** ------------------------------------------------------------------
 ** @brief Get the index of the top level
 ** @param self ::VlScaleSpace object instance..
 ** @return index of the top level.
 **/

VL_INLINE vl_index
vl_scalespace_get_level_max (VlScaleSpace const *self)
{
  assert(self) ;
  return self->lastLevel ;
}

/** ------------------------------------------------------------------
 ** @brief Get the number of stored feature frames.
 ** @param self ::VlScaleSpace object instance.
 ** @return number of frames.
 **/

VL_INLINE vl_size
vl_scalespace_get_frames_num (VlScaleSpace const *self)
{
  assert(self) ;
  return self->numFrames ;
}

/** ------------------------------------------------------------------
 ** @brief Get the stored feature frames.
 ** @param self ::VlScaleSpace object instance.
 ** @return pointer to the frames list.
 **/

VL_INLINE VlScaleSpaceFrame const *
vl_scalespace_get_frames (VlScaleSpace const *self)
{
  assert(self) ;
  return self->frames ;
}

/** ------------------------------------------------------------------
 ** @brief Get stddev of smoothing of pyramid base.
 ** @param self ::VlScaleSpace object instance.
 ** @return Standard deviation of smoothing of pyramid base.
 **/

VL_INLINE double
vl_scalespace_get_sigma0 (VlScaleSpace const *self)
{
  assert(self) ;
  return self->sigma0 ;
}

/** ------------------------------------------------------------------
 ** @brief Get stddev of smoothing of pyramid base.
 ** @param self ::VlScaleSpace object instance.
 ** @return Standard deviation of smoothing of pyramid base.
 **/

VL_INLINE double
vl_scalespace_get_sigmak (VlScaleSpace const *self)
{
  assert(self) ;
  return self->sigmak ;
}


/** ------------------------------------------------------------------
 ** @brief Get continous scale coordinate from scale index
 ** @param self ::VlScaleSpace object instance.
 ** @param o octave.
 ** @param s octave sub-level.
 ** @return contiunous scale coordinate.
 **/

VL_INLINE double
vl_scalespace_get_sigma_for_scale (VlScaleSpace const *self,
                                   vl_index o,
                                   vl_index s)
{
  assert(self) ;
  return self->sigma0 * pow(2.0, (double)s / self->numLevels + o) ;
}

/* ---------------------------------------------------------------- */
/*                                            Affine shape detector */
/* ---------------------------------------------------------------- */

/** ------------------------------------------------------------------
 ** @brief Affine frame
 **
 ** This structure represent a frame as extracted by the Affine
 ** shape detector.
 ** @sa ::VlAffineShapeEstimator.
 **/

typedef struct _VlAffineShapeEstimatorFrame
{
  float x ;     /**< x coordinate in original image. */
  float y ;     /**< y coordinate in original image. */
  double sigma;  /**< scale of the frame. */

  double a11 ;   /**< Affine matrix (ellipse->circle) without scale */
  double a12 ;
  double a21 ;
  double a22 ;

} VlAffineShapeEstimatorFrame ;

/** @brief Affine shape detector
 **
 ** The Affine computes the Affine (elliptic) shape of a frame.
 **
 ** This structure is @ref main-glossary "opaque".
 **
 ** @sa @ref affine
 **/
typedef struct _VlAffineDet
{
  /** @name Inner data planes @internal */
  /** @{ */
  float *img ;        /**< Work buffer for warped neighbourhood. */
  float *mask ;       /**< Gaussian mask for the smm calculation. */
  float *fx ;         /**< x-gradient of the patch. */
  float *fy ;         /**< y-gradient of the patch. */
  /** @} */

  int win_size ;   /**< Size of the window for smm calculation. */
  int win_nel;     /**< Number of pixels in the window */

  /** @name Configuration @internal */
  /** @{ */
  int  maxNumIterations ;   /**< Max. number of iterations. */
  float convergenceThreshold ;   /**< Convergence threshold. */
  /** @} */
} VlAffineShapeEstimator ;


/** @name Create and destroy
 ** @{ */
VL_EXPORT VlAffineShapeEstimator * vl_affineshapeestimator_new (vl_size win_size) ;
VL_EXPORT void vl_affineshapeestimator_delete (VlAffineShapeEstimator *self) ;
/** @} */

/** @name Process data
 ** @{ */
VL_EXPORT int
vl_affineshapeestimator_estimate (VlAffineShapeEstimator *self,
                                  VlScaleSpace* scsp,
                                  const VlScaleSpaceFrame *frame,
                                  VlAffineShapeEstimatorFrame* aff_frame) ;

VL_EXPORT void
vl_affineshapeestimator_frame_init_from_aff (VlAffineShapeEstimator const *self,
                                             VlAffineShapeEstimatorFrame *frm,
                                             double x, double y,
                                             double a11, double a12,
                                             double a21, double a22) ;

VL_EXPORT void
vl_affineshapeestimator_frame_init_from_ell (VlAffineShapeEstimator const *self,
                                             VlAffineShapeEstimatorFrame *frm,
                                             double x, double y,
                                             double e11, double e12, double e22) ;

VL_EXPORT int
vl_affineshapeestimator_interpolate_bilinear (const float * image, vl_size width, vl_size height,
                                              double tx, double ty,
                                              double a11, double a12, double a21, double a22,
                                              float* dst, vl_size dstWidth, vl_size dstHeight) ;
/** @} */

/** @name Retrieve data and parameters
 ** @{ */
VL_INLINE vl_size vl_affineshapeestimator_get_window_size (VlAffineShapeEstimator const *self) ;
VL_INLINE vl_size vl_affineshapeestimator_get_max_iter (VlAffineShapeEstimator const *self) ;
VL_INLINE double vl_affineshapeestimator_get_conv_thresh (VlAffineShapeEstimator const *self) ;

/** @} */

/** @name Set parameters
 ** @{ */
VL_INLINE void vl_affineshapeestimator_set_max_iter (VlAffineShapeEstimator *self, int m) ;
VL_INLINE void vl_affineshapeestimator_set_conv_thresh (VlAffineShapeEstimator *self, double t) ;
/** @} */

/* -------------------------------------------------------------------
 *            Affine shape estimation Inline functions implementation
 * ---------------------------------------------------------------- */

/** @brief Get size of the window for smm calculation.
 ** @param f Affine filter.
 ** @return Size of the window for smm calculation.
 **/

VL_INLINE vl_size
vl_affineshapeestimator_get_window_size (VlAffineShapeEstimator const *f)
{
  return f-> win_size ;
}

/** @brief Get maximal number of iterations.
 ** @param f Affine filter.
 ** @return Maximal number of iterations.
 **/

VL_INLINE vl_size
vl_affineshapeestimator_get_max_iter (VlAffineShapeEstimator const *self)
{
  return self-> maxNumIterations ;
}

/** @brief Get convergance threshold.
 ** @param f Affine filter.
 ** @return Convergence threshold.
 **/

VL_INLINE double
vl_affineshapeestimator_get_conv_thresh (VlAffineShapeEstimator const *self)
{
  return self-> convergenceThreshold ;
}

/** @brief Set maximal number of iterations.
 ** @param f Affine filter.
 ** @param m maximal number of iterations.
 **/

VL_INLINE void
vl_affineshapeestimator_set_max_iter (VlAffineShapeEstimator *self, int m)
{
  self -> maxNumIterations = m ;
}

/** ------------------------------------------------------------------
 ** @brief Set convergence threshold.
 ** @param f Affine filter.
 ** @param t Convergence threshold.
 **/

VL_INLINE void
vl_affineshapeestimator_set_conv_thresh (VlAffineShapeEstimator *self, double t)
{
  self->convergenceThreshold = t ;
}

/* VL_SCALESPACE_H */
#endif


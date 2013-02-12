/** @file imopv.h
 ** @brief Vectorized image operations
 ** @author Andrea Vedaldi
 **/

/*
Copyright (C) 2007-12 Andrea Vedaldi and Brian Fulkerson.
All rights reserved.

This file is part of the VLFeat library and is made available under
the terms of the BSD license (see the COPYING file).
*/

#ifndef VL_IMOPV_H
#define VL_IMOPV_H

#include "generic.h"

/** @name Image convolution flags
 ** @{ */
#define VL_PAD_BY_ZERO       (0x0 << 0) /**< @brief Pad with zeroes. */
#define VL_PAD_BY_CONTINUITY (0x1 << 0) /**< @brief Pad by continuity. */
#define VL_PAD_MASK          (0x3)      /**< @brief Padding field selector. */
#define VL_TRANSPOSE         (0x1 << 2) /**< @brief Transpose result. */
/** @} */

/** @name Image convolution
 ** @{ */
VL_EXPORT
void vl_imconvcol_vf (float* dst, int dst_stride,
                      float const* src,
                      int src_width, int src_height, int src_stride,
                      float const* filt, int filt_begin, int filt_end,
                      int step, unsigned int flags) ;

VL_EXPORT
void vl_imconvcol_vd (double* dst, int dst_stride,
                      double const* src,
                      int src_width, int src_height, int src_stride,
                      double const* filt, int filt_begin, int filt_end,
                      int step, unsigned int flags) ;

VL_EXPORT
void vl_imconvcoltri_f (float * dest, vl_size destStride,
                        float const * image,
                        vl_size imageWidth, vl_size imageHeight, vl_size imageStride,
                        vl_size filterSize,
                        vl_size step, int unsigned flags) ;

VL_EXPORT
void vl_imconvcoltri_d (double * dest, vl_size destStride,
                        double const * image,
                        vl_size imageWidth, vl_size imageHeight, vl_size imageStride,
                        vl_size filterSize,
                        vl_size step, int unsigned flags) ;
/** @} */

/** @name Integral image
 ** @{ */
VL_EXPORT
void vl_imintegral_f (float * integral,  vl_size integralStride,
                      float const * image,
                      vl_size imageWidth, vl_size imageHeight, vl_size imageStride) ;

VL_EXPORT
void vl_imintegral_d (double * integral,  vl_size integralStride,
                      double const * image,
                      vl_size imageWidth, vl_size imageHeight, vl_size imageStride) ;

VL_EXPORT
void vl_imintegral_i32 (vl_int32 * integral,  vl_size integralStride,
                        vl_int32 const * image,
                        vl_size imageWidth, vl_size imageHeight, vl_size imageStride) ;

VL_EXPORT
void vl_imintegral_ui32 (vl_uint32 * integral,  vl_size integralStride,
                         vl_uint32 const * image,
                         vl_size imageWidth, vl_size imageHeight, vl_size imageStride) ;
/** @} */

/** @name Distance transform */
/** @{ */

VL_EXPORT void
vl_image_distance_transform_d (double const * image,
                               vl_size numColumns,
                               vl_size numRows,
                               vl_size columnStride,
                               vl_size rowStride,
                               double * distanceTransform,
                               vl_uindex * indexes,
                               double coeff,
                               double offset) ;

VL_EXPORT void
vl_image_distance_transform_f (float const * image,
                               vl_size numColumns,
                               vl_size numRows,
                               vl_size columnStride,
                               vl_size rowStride,
                               float * distanceTransform,
                               vl_uindex * indexes,
                               float coeff,
                               float offset) ;

/** @} */


/** @name Image smoothing */
/** @{ */

/** @internal @brief Image smoothing filter
 **/
typedef struct _VlImSmoothFilt
{
  int       width ;             /**< image width. */
  int       height ;            /**< image height. */

  float    *temp ;              /**< temporary pixel buffer. */
  size_t    temp_sz;

  float    *gaussFilter ;       /**< current Gaussian filter */
  double    gaussFilterSigma ;  /**< current Gaussian filter std */
  vl_size   gaussFilterWidth ;  /**< current Gaussian filter width */

} VlImSmoothFilt ;


VL_EXPORT
VlImSmoothFilt*  vl_imsmooth_new    (int width, int height) ;

VL_EXPORT
void         vl_imsmooth_delete (VlImSmoothFilt *f) ;

VL_EXPORT
void
vl_imsmooth_smooth_sub_image (VlImSmoothFilt *self, float *outputImage, vl_size outStride,
                     float const *inputImage, vl_size inWidth, vl_size inHeight,
                     vl_size inStride, double sigma);


VL_INLINE
void
vl_imsmooth_smooth_image (VlImSmoothFilt *self, float *outputImage,
                 float const *inputImage, vl_size inWidth, vl_size inHeight,
                 double sigma);

/** @} */


/** @name Image gradients */
/** @{ */

VL_EXPORT
void
vl_imgradient_polar_f (float* amplitudeGradient, float* angleGradient,
                       vl_size gradWidthStride, vl_size gradHeightStride,
                       float const* image,
                       vl_size imageWidth, vl_size imageHeight,
                       vl_size imageStride);

VL_EXPORT
void
vl_imgradient_polar_d (double* amplitudeGradient, double* angleGradient,
                       vl_size gradWidthStride, vl_size gradHeightStride,
                       double const* image,
                       vl_size imageWidth, vl_size imageHeight,
                       vl_size imageStride);

VL_EXPORT
void
vl_imgradient_f (float* xGradient, float* yGradient,
                 vl_size gradWidthStride, vl_size gradHeightStride,
                 float* image,
                 vl_size imageWidth, vl_size imageHeight, vl_size imageStride);

VL_EXPORT
void
vl_imgradient_d(double* xGradient, double* yGradient,
                vl_size gradWidthStride, vl_size gradHeightStride,
                double* image,
                vl_size imageWidth, vl_size imageHeight, vl_size imageStride);

VL_EXPORT
void
vl_imgradient_polar_f_callback(float const *sourceImage,
                                 int sourceImageWidth, int sourceImageHeight,
                               float *dstImage,
                               int dstWidth, int dstHeight,
                               int octave, int level,
                               void *params);

/** @} */


/* -------------------------------------------------------------------
 *                                     Inline functions implementation
 * ---------------------------------------------------------------- */


/** ------------------------------------------------------------------
 ** @brief Smooth an image with a gaussian filter.
 **
 ** This functions smoothes inputImage with gaussian filter with standard
 ** deviation sigma into outputImage.
 **
 ** Please note that the filter is faster if images are smaller or equal
 ** size defined in vl_smooth_new function.
 **
 ** @param self        Smooth filter.
 ** @param outputImage output image buffer.
 ** @param inputImage  input image buffer.
 ** @param inWidth    input image width.
 ** @param inHeight   input image height.
 ** @param sigma       smoothing.
 **/

VL_INLINE void
vl_imsmooth_smooth_image(VlImSmoothFilt *self, float *outputImage,
                const float *inputImage, vl_size inWidth, vl_size inHeight,
                double sigma)
{
  vl_imsmooth_smooth_sub_image(self, outputImage, inWidth,
                      inputImage, inWidth, inHeight, inWidth,
                      sigma);
}


/* VL_IMOPV_H */
#endif

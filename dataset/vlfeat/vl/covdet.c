/** @file framedet.c
 ** @brief Covariant feature detectors - Definition
 ** @author Andrea Vedaldi
 ** @author Karel Lenc
 **/

/*
Copyright (C) 2007-12 Andrea Vedaldi and Brian Fulkerson.
All rights reserved.

This file is part of the VLFeat library and is made available under
the terms of the BSD license (see the COPYING file).
*/

/**
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@page framedet Covariant image frames detector
@author Andrea Vedaldi
@author Karel Lenc
@par "Credits:" May people have contributed with suggestions and bug
reports. Although the following list is certainly incomplete, we would
like to thank: Wei Dong, Loic, Giuseppe, Liu, Erwin, P. Ivanov, and
Q. S. Luo.
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

@ref framedet.h implements a @ref framedet-usage "Covariant frames detector", a
reusable object to extract covariant image features (in the following text
referred as frames) @cite{lowe99object} @cite{mikolajczyk04scaff} from one or
multiple images.

- @ref framedet-intro
  - @ref framedet-intro-frame-classes
  - @ref framedet-intro-frame-detection
  - @ref framedet-intro-disc-detection
  - @ref framedet-intro-ellipse-detection
  - @ref framedet-intro-descriptor
- @ref framedet-usage
  - @ref framedet-usage-sift
  - @ref framedet-usage-hessaff
  - @ref framedet-usage-conversion
- @ref framedet-tech
  - @ref framedet-tech-disc-detector
    - @ref framedet-tech-disc-peak
    - @ref framedet-tech-disc-edge
    - @ref framedet-tech-disc-orientation
    - @ref framedet-tech-disc-hessian
  - @ref framedet-tech-descriptor
    - @ref framedet-tech-descriptor-can
    - @ref framedet-tech-descriptor-image
    - @ref framedet-tech-descriptor-std
  - @ref framedet-tech-ellipse-detector
    - @ref framedet-tech-ellipse-frm
    - @ref framedet-tech-ellipse-affshape
    - @ref framedet-tech-ellipse-smm
    - @ref framedet-tech-ellipse-trsmm
    - @ref framedet-tech-ellipse-covm
    - @ref framedet-tech-ellipse-iter
    - @ref framedet-tech-ellipse-norm

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@section framedet-intro Overview
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

Image frame is a set of attributes for a selected image region (also called
keypoint) usually with an associated descriptor. These frames are extraced
using <b>@ref framedet-intro-detector "Image frames detector"</b> and
extracted frames can vary in its covariance to image transformation depending
on the class of the detected frames. The frame detection can also differ in
the image response function which was used for detection of distinct image
regions. Special case of these frames are <b>@ref framedet-usage-sift
SIFT frames</b>. All the detected frames can be accompanied with SIFT
descriptors.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsection framedet-intro-frame-classes Image frame classes
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@image html frame-types.png "Classes of image frames."

Image frames can differ in the attributes which are computed for the
selected image region (blob). Based on these attributes these frames gain
covariance to certain image transformations which are summarised in the
following table:

<table>
 <caption>Attributes of image frame classes.</caption>
 <tr style="font-weight:bold;">
 <td>Frame class</td>
 <td>Frame attributes</td>
 <td>Covariant to</td>
 <td>See also</td>
 </tr>
 <tr>
 <td> Disc</td>
 <td> Frame coordinates \f$ (x,y) \f$, scale (radius) \f$ \sigma \f$ </td>
 <td> Translations and scalings</td>
 <td> @ref framedet-intro-disc-detection</td>
 </tr>
 <tr>
 <td> Oriented Disc</td>
 <td> Frame coordinates \f$ (x,y) \f$, scale \f$ \sigma \f$ and orientation \f$ \theta \f$ </td>
 <td> Translations, scalings and rotation (similarities) </td>
 <td> @ref framedet-intro-disc-detection</td>
 </tr>
 <tr>
 <td> Ellipse</td>
 <td> Frame coordinates \f$ (x,y) \f$ and shape matrix \f$ E \f$ with three degrees of freedom </td>
 <td> Translation and affinities up to residual rotations</td>
 <td> @ref framedet-intro-ellipse-detection</td>
 </tr>
 <tr>
 <td> Oriented Ellipse</td>
 <td> Frame coordinates \f$ (x,y) \f$ and affine transfomration \f$ A \f$ between ellipse and circle </td>
 <td> Translation and affinities up to residual rotations</td>
 <td> @ref framedet-intro-ellipse-detection</td>
 </tr>
</table>

The name of the frame classes are selected in order to get graphical intuition
of the blob structure and they share the way how they behave under linear
transformations.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsection framedet-intro-frame-detection Image frames detection
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

Covariant frames detector can be configured to detect any of the mentioned classes.
Because only some of the attributes are shared among classes the
detection of each frame type differs.

@image html frame-detection.png "Pipeline of frames detection."

In all cases the detection start by detecting disc frames using the
::VlScaleSpace object. Detecting the major orientations extends the
disc frame into several oriented disc frames.

If the detector is configured to detect ellipses, image region of each
disc is examined for its affine shape based on its second moment matrix.
In order to assign orientations, local anisotropic structures has to
be transformed into isotrpic circular regions. This step, in fact, transforms
the elliptic blob into a circular blob and thereafter the orientation can be
assigned in the same way as for disc frames.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsection framedet-intro-disc-detection Disc detection
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

@sa
@ref scalespace-tech-ss "Scale space technical details",
@ref framedet-tech-detector "Detector technical details"

A disc frame is a circular image region and in case of oriented disx
with assinged orientation. It is described by a geometric <em>frame</em>
of three/four parameters: the keypoint center coordinates @e x and @e y,
its @e scale (the radius of the region), and its @e orientation (an angle
expressed in radians). The disc detector uses as keypoints image
structures which resemble &ldquo;blobs&rdquo;. By searching for blobs
at multiple scales and positions, the disc detector is invariant (or,
more accurately, covariant) to translation, rotations, and rescaling
of the image.

The keypoint orientation is also determined from the local image
appearance and is covariant to image rotations. Depending on the
symmetry of the keypoint appearance, determining the orientation can
be ambiguous. In this case, the Covariant frames detector returns a list of up to
four possible orientations, constructing up to four frames (differing
only by their orientation) for each detected image blob.

@image html sift-frame.png "Disc frames are circular image regions with an orientation."

There are several parameters that influence the detection of disc
keypoints. First, searching keypoints at multiple scales is obtained
by constructing a so-called &ldquo;Gaussian scale space&rdquo;. The
scale space is just a collection of images obtained by progressively
smoothing the input image, which is analogous to gradually reducing
the image resolution. Conventionally, the smoothing level is called
<em>scale</em> of the image. The construction of the scale space is
influenced by the following parameters, set when creating the SIFT
filter object by ::vl_sift_new():

- <b>Number of octaves</b>. Increasing the scale by an octave means
  doubling the size of the smoothing kernel, whose effect is roughly
  equivalent to halving the image resolution. By default, the scale
  space spans as many octaves as possible (i.e. roughly <code>
  log2(min(width,height)</code>), which has the effect of searching
  keypoints of all possible sizes.
- <b>First octave index</b>. By convention, the octave of index 0
  starts with the image full resolution. Specifying an index greater
  than 0 starts the scale space at a lower resolution (e.g. 1 halves
  the resolution). Similarly, specifying a negative index starts the
  scale space at an higher resolution image, and can be useful to
  extract very small features (since this is obtained by interpolating
  the input image, it does not make much sense to go past -1).
- <b>Number of levels per octave</b>. Each octave is sampled at this
  given number of intermediate scales (by default 3). Increasing this
  number might in principle return more refined keypoints, but in
  practice can make their selection unstable due to noise (see [1]).

Each image of the scale-space is then transformed using different image
response functions in order to detect corner-like blobs in the image
which are found to be stable under several image transformations. The
Covariant frames detector supports the following response functions:

- <b>Difference of Gaussian</b>. Computionaly effective approximation
  of Laplacian of Gaussian (LoG) - scalar value bases on the second order
  derivative of the image brightness function.
- <b>Hessian response</b>. Image response based on determinant of
  hessian matrix which is constructed from second order partial
  derivatives of the image brigtness. Contrary to LoG it also counts
  with mixed derivatives.

Based on the response function new gaussian scale-space is constructed
where the distinct image regions (keypoints) are detected as local
extrema of the response function across spatial coordinates an scale.
Keypoints are further refined by eliminating those that are likely to
be unstable, either because they are selected nearby an image edge,
rather than an image blob, or are found on image structures with low
contrast. Filtering is controlled by the follow:

- <b>Peak threshold.</b> This is the minimum amount of contrast to
  accept a keypoint. It is set by configuring the SIFT filter object
  by ::vl_covdet_set_peak_thresh().
- <b>Edge threshold.</b> This is the edge rejection threshold. It is
  set by configuring the SIFT filter object by
  ::vl_covdet_set_edge_thresh().

<table>
 <caption>Summary of the parameters influencing the disc Covariant frames detector.</caption>
 <tr style="font-weight:bold;">
 <td>Parameter</td>
 <td>See also</td>
 <td>Controlled by</td>
 <td>Comment</td>
 </tr>
 <tr>
 <td>image repsonse function</td>
 <td> @ref framedet-intro-disc-detection </td>
 <td>::vl_covdet_new_disc_detector</td>
 <td>Values: @c VL_IMRESP_DOG or @c VL_IMRESP_HESSIAN</td>
 </tr>
 <tr>
 <td>number of octaves</td>
 <td> @ref framedet-intro-disc-detection </td>
 <td>::vl_covdet_new_disc_detector</td>
 <td></td>
 </tr>
 <tr>
 <td>first octave index</td>
 <td> @ref framedet-intro-disc-detection </td>
 <td>::vl_covdet_new_disc_detector</td>
 <td>set to -1 to extract very small features</td>
 </tr>
 <tr>
 <td>number of scale levels per octave</td>
 <td> @ref framedet-intro-disc-detection </td>
 <td>::vl_covdet_new_disc_detector</td>
 <td>can affect the number of extracted keypoints</td>
 </tr>
 <tr>
 <td>edge threshold</td>
 <td> @ref framedet-intro-disc-detection </td>
 <td>::vl_covdet_set_edge_thresh</td>
 <td>decrease to eliminate more keypoints</td>
 </tr>
 <tr>
 <td>peak threshold</td>
 <td> @ref framedet-intro-disc-detection </td>
 <td>::vl_covdet_set_peak_thresh</td>
 <td>increase to eliminate more keypoints</td>
 </tr>
</table>


<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsection framedet-intro-ellipse-detection Ellipse detection
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

As it was said before, oriented discs are covariant only to similarities.
Although in the case of the viewpoint change the image blob not only changes
its scale and position but also its affine shape. This is why attributes
covariant to affine transformation are introduced.

Estimation of the affine shape of the detected disc frame is based on
iterative procedure based on the Second Moment Matrix (SMM). This matrix
is used as a local image texture descriptor which is covariant to affine
tranformation of the Image domain.

The used iterative procedure sequentialy tries to estimate the affine shape
of the frame region by computing an affine transformation which transforms
the local anisotrpic structure into an isotropic one. Using the graphic
representations it tries to transform the elliptic region into a circular one.
The size of the regions which is used for calculating the SMM can be controlled
by parameter <b> Affine Window Size </b>.

The procedure stops when the isotropy measure is closer to ideal measure of
circular region than the parameter <b> convergence threshold </b>. The
procedure is also limited by the maximal number of iterations, controlled
by parameter <b> Maximum Iterations </b>.

In order to obtain affine invariant descriptors and residual rotation these
attributes are calculated from the transformed isotropic structure. In order
to obtain higher precission, the descriptors are usually calculated from larger
windows then the windows used for the affine shape estimation. The size of
these windows can be controlled by parameter <b> Descriptor Window Size </b>.

Because the parameters <b> Affine Window Size </b> and
<b> Descriptor Window Size </b> influence the size of the allocated memory
they have to be defined in the constructor of the Covariant frames detector.

<table>
 <caption>Summary of the parameters influencing the disc Covariant frames detector.</caption>
 <tr style="font-weight:bold;">
 <td>Parameter</td>
 <td>See also</td>
 <td>Controlled by</td>
 <td>Comment</td>
 </tr>
 <tr>
 <td colspan="4" align="center">
 <i>Disc detector parameters, see @ref framedet-intro-disc-detection</i>
 </td>
 </tr>
 <tr>
 <td>Size of the window for affine shape estimation</td>
 <td> @ref framedet-intro-ellipse-detection </td>
 <td>::vl_covdet_new_ellipse_detector</td>
 <td>affine shape estimation precision</td>
 </tr>
 <tr>
 <td>Size of the window for descriptor and orientation calculation</td>
 <td> @ref framedet-intro-ellipse-detection </td>
 <td>::vl_covdet_new_ellipse_detector</td>
 <td>descriptor calculation precision</td>
 </tr>
 <tr>
 <td>Maximum number of iteration for aff. shape est.</td>
 <td> @ref framedet-intro-ellipse-detection </td>
 <td>::vl_affineshapeestimator_set_max_iter</td>
 <td>affine shape estimation precision</td>
 </tr>
 <tr>
 <td>Convergence criterion threshold of aff. shape est.</td>
 <td> @ref framedet-intro-ellipse-detection </td>
 <td>::vl_affineshapeestimator_set_conv_thresh</td>
 <td>affine shape estimation precision</td>
 </tr>
</table>

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsection framedet-intro-descriptor SIFT Descriptor
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

@sa @ref framedet-tech-descriptor "Descriptor technical details"

A SIFT descriptor is a 3-D spatial histogram of the image gradients in
characterizing the appearance of a keypoint. The gradient at each
pixel is regarded as a sample of a three-dimensional elementary
feature vector, formed by the pixel location and the gradient
orientation. Samples are weighed by the gradient norm and accumulated
in a 3-D histogram @em h, which (up to normalization and clamping)
forms the SIFT descriptor of the region. An additional Gaussian
weighting function is applied to give less importance to gradients
farther away from the keypoint center. Orientations are quantized into
eight bins and the spatial coordinates into four each, as follows:

@image html sift-descr-easy.png "The SIFT descriptor is a spatial histogram of the image gradient."

SIFT descriptors are computed by either calling
::vl_sift_calc_keypoint_descriptor or
::vl_sift_calc_raw_descriptor. They accept as input a disc
frame, which specifies the descriptor center, its size, and its
orientation on the image plane. In the case of elliptic frames
the descriptor is calculated from the affine-normalised image
region of the image based on the affine shape of the frame.

The following parameters influence the descriptor calculation:

- <b>magnification factor</b>. The descriptor size is determined by
multiplying the keypoint scale by this factor. It is set by
::vl_sift_set_magnif.
- <b>Gaussian window size</b>. The descriptor support is determined by
a Gaussian window, which discounts gradient contributions farther away
from the descriptor center. The standard deviation of this window is
set by ::vl_sift_set_window_size and expressed in unit of bins.

VLFeat SIFT descriptor uses the following convention. The @em y axis
points downwards and angles are measured clockwise (to be consistent
with the standard image convention). The 3-D histogram (consisting of
@f$ 8 \times 4 \times 4 = 128 @f$ bins) is stacked as a single
128-dimensional vector, where the fastest varying dimension is the
orientation and the slowest the @em y spatial coordinate. This is
illustrated by the following figure.

@image html sift-conv-vlfeat.png "VLFeat conventions"

@note Keypoints (frames) D. Lowe's SIFT implementation convention is
slightly different: The @em y axis points upwards and the angles are
measured counter-clockwise.

@image html sift-conv.png "D. Lowes' SIFT implementation conventions"

<table>
 <caption>Summary of the parameters influencing the SIFT descriptor.</caption>
 <tr style="font-weight:bold;">
 <td>Parameter</td>
 <td>See also</td>
 <td>Controlled by</td>
 <td>Comment</td>
 </tr>
 <tr>
 <td>magnification factor</td>
 <td> @ref sift-intro-descriptor </td>
 <td>::vl_siftdesc_set_magnif</td>
 <td>increase this value to enlarge the image region described</td>
 </tr>
 <tr>
 <td>Gaussian window size</td>
 <td> @ref sift-intro-descriptor </td>
 <td>::vl_siftdesc_set_window_size</td>
 <td>smaller values let the center of the descriptor count more</td>
 </tr>
</table>


<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@section framedet-intro-extensions Extensions
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

<b>Eliminating low-contrast descriptors.</b> Near-uniform patches do
not yield stable keypoints or descriptors. ::vl_sift_set_norm_thresh()
can be used to set a threshold on the average norm of the local
gradient to zero-out descriptors that correspond to very low contrast
regions. By default, the threshold is equal to zero, which means that
no descriptor is zeroed. Normally this option is useful only with
custom keypoints, as detected keypoints are implicitly selected at
high contrast image regions.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@section framedet-usage Using the VlFrameDet object
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

The code provided in this module can be used in different ways.  You
can instantiate and use a <b>VlFrameDet</b> object to extract
arbitrary type of frames and descriptors from one or multiple images.

To use a <b>VlFrameDet</b> object:

- Initialize a VlFrameDet object with ::vl_covdet_new_ellipse_detector
  or ::vl_covdet_new_disc_detector. The filter can
  be reused for multiple images of the same size (e.g. for an entire
  video sequence).
- Call ::vl_covdet_detect() to detect frames in an image. This
  function also calculates descriptors if it was defined in the constructor.
- Number of detected keypoints can be accesed by ::vl_covdet_get_num_frames
- With methods ::vl_covdet_get_discs, ::vl_covdet_get_oriented_discs,
  ::vl_covdet_get_ellipses or ::vl_covdet_get_oriented_ellipses. Or
  the array storing frames can be accessed directly with method
  ::vl_covdet_get_frames_storage and the type of the frames by
  ::vl_covdet_get_frames_type
- Descriptors can be accessed with method ::vl_covdet_get_descriptors(),
  the size of a frame descriptor is ::vl_covdet_get_descriptor_size()
- Delete the SIFT filter by ::vl_covdet_delete().

Please note that for each image the frame storage can be reallocated therefore
the frames can be stored in different mamory location. That is why the new
pointer to frames and descriptors should be get after each detection.

To compute SIFT descriptors of custom keypoints, use
::vl_sift_calc_descriptor().

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsection framedet-usage-sift SIFT keypoint detector
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

The ::VlFrameDet can be used for example as a SIFT frames detector by construction
of the ::VlFrameDet object and settings the parameters as follows:

@code{.c}
VlImRespFunction respFunction = VL_IMRESP_DOG;
vl_bool calcOrientation = VL_TRUE;
vl_bool calcDescriptor = VL_TRUE;
VlFrameDet siftDet = vl_covdet_new_disc_detector(w, h, respFunction,
                                                   O, firstOctave, S,
                                                   calcOrientation, calcDescriptor);
@endcode

Then the detection is fairly simple:

@code{.c}
vl_covdet_detect(siftDet, image);

vl_size framesNum = vl_covdet_get_frames_num(siftDet);
VlFrameOrientedDisc const *frames = vl_covdet_get_oriented_discs(siftDet);

vl_size descriptorSize = vl_covdet_get_descriptor_size(siftDet);
float const *descriptors = vl_covdet_get_descriptors(siftDet);
@endcode

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsection framedet-usage-hessaff Hessian-Affine keypoint detector
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

Similalrly the Hessian-Affine detector can be contructed as well
simply by constructing the frames detector with the following parameters:

@code{.c}
VlImRespFunction respFunction = VL_IMRESP_HESSIAN;
vl_bool calcOrientation = VL_TRUE;
vl_bool calcDescriptor = VL_TRUE;
VlFrameDet hessaffDet = vl_covdet_new_ellipse_detector(w, h, respFunction,
                                                         O, firstOctave, S,
                                                         calcOrientation, calcDescriptor);
@endcode

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsection framedet-usage-conversion Keypoint conversion
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

The constructed ::VlFrameDet object can be also used for conversion
between different classes of frames using the method
::vl_covdet_convert_frames() where the output class of frames
is specified by the class of frames detected by the detector.
The missing attributes if the frames are computed in the same
way as they are computed by the detector. The following code
shows how to convert disc frames into oriented ellipses and to
calculate their decriptors:

@code{.c}
VlFrameOrientedEllipse* dstFrames;
float const *dstDescriptors;

VlImRespFunction respFunction = VL_IMRESP_DOG; // Irrelevant value
vl_bool dstCalcOrientation = VL_TRUE;
vl_bool dstCalcDescriptor = VL_TRUE;
VlFrameDet detector = vl_covdet_new_ellipse_detector(w, h, respFunction,
                         O, firstOctave, S, dstCalcOrientation, dstCalcDescriptor);

VlFrameDisc *srcFrames = ...
vl_size srcFramesNum = ...
VlFrameType srcFrameType = VL_FRAME_DISC;
vl_covdet_convert_frames(detector, image, (const void*)srcFrames, srcFramesNum,
                           srcFramesType);

dstFrames = vl_covdet_get_oriented_ellipses(detector);
dstDescriptors = vl_covdet_get_descriptors(detector);

@endcode


<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@section framedet-tech Technical details
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsection framedet-tech-disc-detector Disc Detector
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

The SIFT frames (keypoints) are extracted based on local extrema
(peaks) of the DoG scale space. Numerically, local extrema are
elements whose @f$ 3 \times 3 \times 3 @f$ neighbors (in space and
scale) have all smaller (or larger) value.  Once extracted, local
extrema are quadratically interpolated (this is very important
especially at the lower resolution scales in order to have accurate
keypoint localization at the full resolution).  Finally, they are
filtered to eliminate low-contrast responses or responses close to
edges and the orientation(s) are assigned, as explained next.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsubsection framedet-tech-disc-peak Eliminating low contrast responses
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

Peaks which are too short may have been generated by noise and are
discarded.  This is done by comparing the absolute value of the DoG
scale space at the peak with the <b>peak threshold</b> @f$t_p@f$ and
discarding the peak its value is below the threshold.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsubsection framedet-tech-disc-edge Eliminating edge responses
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

Peaks which are too flat are often generated by edges and do not yield
stable features. These peaks are detected and removed as follows.
Given a peak @f$x,y,\sigma@f$, the algorithm evaluates the @em x,@em y
Hessian of of the DoG scale space at the scale @f$\sigma@f$.  Then the
following score (similar to the Harris function) is computed:

@f[
\frac{(\mathrm{tr}\,D(x,y,\sigma))^2}{\det D(x,y,\sigma)},
\quad
D =
\left[
\begin{array}{cc}
\frac{\partial^2 \mathrm{DoG}}{\partial x^2} &
\frac{\partial^2 \mathrm{DoG}}{\partial x\partial y} \\
\frac{\partial^2 \mathrm{DoG}}{\partial x\partial y} &
\frac{\partial^2 \mathrm{DoG}}{\partial y^2}
\end{array}
\right].
 @f]

This score has a minimum (equal to 4) when both eigenvalues of the
Jacobian are equal (curved peak) and increases as one of the
eigenvalues grows and the other stays small. Peaks are retained if the
score is below the quantity @f$(t_e+1)(t_e+1)/t_e@f$, where @f$t_e@f$
is the <b>edge threshold</b>. Notice that this quantity has a minimum
equal to 4 when @f$t_e=1@f$ and grows thereafter. Therefore the range
of the edge threshold is @f$[1,\infty)@f$.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsection framedet-tech-disc-orientation Orientation assignment
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

A peak in the DoG scale space fixes 2 parameters of the keypoint: the
position and scale. It remains to choose an orientation. In order to
do this, SIFT computes an histogram of the gradient orientations in a
Gaussian window with a standard deviation which is 1.5 times bigger
than the scale @f$\sigma@f$ of the keypoint.

@image html sift-orient.png

This histogram is then smoothed and the maximum is selected.  In
addition to the biggest mode, up to other three modes whose amplitude
is within the 80% of the biggest mode are retained and returned as
additional orientations.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsubsection framedet-tech-disc-hessian Hessian keypoint detector
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

Scale adapted hessian matrix is defined as:
\f[
  H(\mathbf{x}, \sigma_D) =
  \sigma_D^2
  \left[\begin{array}{cc}
  L_{xx}(\mathbf{x}, \sigma_D) & L_{xy}(\mathbf{x}, \sigma_D) \\
  L_{xy}(\mathbf{x}, \sigma_D) & L_{yy}(\mathbf{x}, \sigma_D)
  \end{array}\right]
\f]
Where \f$ \sigma_D \f$ is scale of the current level.

Where \f$ L(\mathbf{x}, \sigma) \f$ is given as:
\f[
  L(\mathbf{x}, \sigma) = g(\mathbf{x}, \sigma) * I(\mathbf{x})
\f]
Where \f$ I \f$ stands for the image intensities.

And \f$ L_{ab} \f$ stands for second derivations as:
\f[
  L_{ab} = \frac{\partial L}{\partial a} \frac{\partial L}{\partial b}
\f]

Note also that normalisation factor \f$ \sigma_{D} \f$ is needed due to
maintain scale invariance of the response. Scale normalised derivative
of order \f$ i \f$ is defined as:
\f[
  \sigma^{m} L_{i_1 \dotsc i_m}(\mathbf{x},\sigma) =
  \sigma^{m} g_{i_1 \dotsc i_m}(\mathbf{x},\sigma) * I(\mathbf{x})
\f]
Due to the properties of partial derivative of Gaussian function
\cite{mikolajczyk01index}.

Having the scale adapted Hessian matrix, the Hessian response is given as
determinant of the Hessian matrix:
\f[
  \mathrm{Resp}(\sigma_D) = |H(\sigma_D))| =
  \sigma_D^4(L_{xx}(\sigma_D) L_{yy}(\sigma_D) - L_{xy}^2(\sigma_D))
\f]

In the used implementation the Hessian matrix is calculated only in window
of \f$ 3 \times 3 \f$ pixels and the second order partial derivatives
are approximated as follows:

\f[
  L_{xx} = \left[\begin{array}{ccc}1 & -2 & 1\end{array}\right] * L \quad
  L_{xx} = \left[\begin{array}{ccc}
              -\frac{1}{4} & 0 &  \frac{1}{4} \\
              0            & 0 & 0            \\
               \frac{1}{4} & 0 & -\frac{1}{4} \\
           \end{array}\right] * L \quad
  L_{yy} = \left[\begin{array}{c}1 \\ -2 \\ 1\end{array}\right] * L
\f]

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsection framedet-tech-descriptor Descriptor
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

A SIFT descriptor of a local region (keypoint) is a 3-D spatial
histogram of the image gradients.  The gradient at each pixel is
regarded as a sample of a three-dimensional elementary feature vector,
formed by the pixel location and the gradient orientation. Samples are
weighed by the gradient norm and accumulated in a 3-D histogram @em h,
which (up to normalization and clamping) forms the SIFT descriptor of
the region. An additional Gaussian weighting function is applied to
give less importance to gradients farther away from the keypoint
center.

<!-- @image html sift-bins.png -->

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsubsection framedet-tech-descriptor-can Construction in the canonical frame
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

Denote the gradient vector field computed at the scale @f$ \sigma @f$ by
@f[
  J(x,y) = \nalba I_\sigma(x,y)
  =
  \left[\begin{array}{cc}
  \frac{\partial I_\sigma}{\partial x} &
  \frac{\partial I_\sigma}{\partial y} &
  \end{array}\right]
@f]

The descriptor is a 3-D spatial histogram capturing the distribution
of @f$ J(x,y) @f$. It is convenient to describe its construction in
the <em>canonical frame</em>. In this frame, the image and descriptor
axes coincide and each spatial bin has side 1. The histogram has @f$
N_\theta \times N_x \times N_y @f$ bins (usually @f$ 8 \times 4 \times
4 @f$), as in the following figure:

@image html sift-can.png Canonical SIFT descriptor and spatial binning functions

Bins are indexed by a triplet of indexes <em>t, i, j</em> and their
centers are given by

@f{eqnarray*}
 \theta_t &=& \frac{2\pi}{N_\theta} t, \quad t = 0,\dots,N_{\theta}-1, \\
 x_i &=& i - \frac{N_x -1}{2}, \quad i = 0,\dots,N_x-1, \\
 y_j &=& j - \frac{N_x -1}{2}, \quad j = 0,\dots,N_y-1. \\
@f}

The histogram is computed by using trilinear interpolation, i.e.  by
weighing contributions by the <em>binning functions</em>

@f{eqnarray*}
  \displaystyle
  w(z) &=& \mathrm{max}(0, 1 - |z|),
  \\
  \displaystyle
  w_\mathrm{ang}(z) &=& \sum_{k=-\infty}^{+\infty}
  w\left(
  \frac{N_\theta}{2\pi} z + N_\theta k
  \right).
@f}

The gradient vector field is transformed in a three-dimensional
density map of weighed contributions

@f[
   f(\theta, x, y) = |J(x,y)| \delta(\theta - \angle J(x,y))
@f]

The historam is localized in the keypoint support by a Gaussian window
of standard deviation @f$ \sigma_{\mathrm{win}} @f$. The histogram is
then given by

@f{eqnarray*}
 h(t,i,j) &=& \int
 g_{\sigma_\mathrm{win}}(x,y)
 w_\mathrm{ang}(\theta - \theta_t) w(x-x_i) w(y-y_j)
 f(\theta,x,y)
 d\theta\,dx\,dy
\\
&=& \int
 g_{\sigma_\mathrm{win}}(x,y)
 w_\mathrm{ang}(\angle J(x,y) - \theta_t) w(x-x_i) w(y-y_j)
 |J(x,y)|\,dx\,dy
@f}

In post processing, the histogram is @f$ l^2 @f$ normalized, then
clamped at 0.2, and @f$ l^2 @f$ normalized again.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsubsection framedet-tech-descriptor-image Calculation in the image frame
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

Invariance to similarity transformation is attained by attaching
descriptors to SIFT keypoints (or other similarity-covariant frames).
Then projecting the image in the canonical descriptor frames has the
effect of undoing the image deformation.

In practice, however, it is convenient to compute the descriptor
directly in the image frame. To do this, denote with a hat quantities
relative to the canonical frame and without a hat quantities relative
to the image frame (so for instance @f$ \hat x @f$ is the @e
x-coordinate in the canonical frame and @f$ x @f$ the x-coordinate in
the image frame). Assume that canonical and image frame are
related by an affinity:

@f[
  \mathbf{x} = A \hat\mathbf{x} + T,
  \qquad
  \mathbf{x} =
  \left[\begin{array}{cc}
    x \\
    y
  \end{arraty}\right],
  \quad
  \hat\mathbf{x} =
  \left[\begin{array}{cc}
    \hat x \\
    \hat y
  \end{arraty}\right].
@f]

@image html sift-image-frame.png

Then all quantities can be computed in the image frame directly. For instance,
the image at infinite resolution in the two frames are related by

@f[
 \hat I_0(\hat{\mathbf{x}})  = I_0(\mathbf{x}),
 \qquad
 \mathbf{x} = A \hat{\mathbf{x}} + T.
@f]

The canonized image at scale @f$ \hat \sigma @f$ is in relation with the scaled image

@f[
 \hat I_{\hat{\sigma}}(\hat{\mathbf{x}})  = I_{A\hat{\sigma}}(\mathbf{x}),
 \qquad \mathbf{x} = A \hat{\mathbf{x}} + T
@f]

where, by generalizing the previous definitions, we have

@f[
 I_{A\hat \sigma}(\mathbf{x}) = (g_{A\hat\sigma} * I_0)(\mathbf{x}),
\quad
 g_{A\hat\sigma}(\mathbf{x})
 =
 \frac{1}{2\pi|A|\hat \sigma^2}
 \exp
 \left(
 -\frac{1}{2}
 \frac{\mathbf{x}^\top A^{-\top}A^{-1}\mathbf{x}}{\hat \sigma^2}
 \right)
@f]

Deriving shows that the gradient fields are in relation

@f[
  \hat J(\hat{\mathbf{x}}) = J(\mathbf{x}) A,
 \quad J(\mathbf{x}) = (\nabla I_{A\hat\sigma})(\mathbf{x}),
 \qquad \mathbf{x} = A \hat{\mathbf{x}} + T.
@f]

Therefore we can compute the descriptor either in the image or canonical frame as:

@f{eqnarray*}
 h(t,i,j)
 &=&
 \int
 g_{\hat \sigma_\mathrm{win}}(\hat{\mathbf{x}})\,
 w_\mathrm{ang}(\angle \hat J(\hat{\mathbf{x}}) - \theta_t)\,
 w_{ij}(\hat{\mathbf{x}})\,
 |\hat J(\hat{\mathbf{x}})|\,
 d\hat{\mathbf{x}}
 \\
 &=& \int
 g_{A \hat \sigma_\mathrm{win}}(\mathbf{x} - T)\,
 w_\mathrm{ang}(\angle J(\mathbf{x})A - \theta_t)\,
 w_{ij}(A^{-1}(\mathbf{x} - T))\,
 |J(\mathbf{x})A|\,
 d\mathbf{x}.
@f}

where we defined the product of the two spatial binning functions

@f[
 w_{ij}(\hat{\mathbf{x}}) = w(\hat x - \hat x_i) w(\hat y - \hat y_j)
@f]


In the actual implementation, this integral is computed by visiting a
rectangular area of the image that fully contains the keypoint grid
(along with half a bin border to fully include the bin windowing
function). Since the descriptor can be rotated, this area is a
rectangle of sides @f$m/2\sqrt{2} (N_x+1,N_y+1)@f$ (see also the
illustration).

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsubsection framedet-tech-descriptor-std Standard SIFT descriptor
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

For a SIFT-detected keypoint of center @f$ T @f$, scale @f$ \sigma @f$
and orientation @f$ \theta @f$, the affine transformation @f$ (A,T)
@f$ reduces to the similarity transformation

@f[
     \mathbf{x} = m \sigma R(\theta) \hat{\mathbf{x}} + T
@f]

where @f$ R(\theta) @f$ is a counter-clockwise rotation of @f$ \theta
@f$ radians, @f$ m \mathcal{\sigma} @f$ is the size of a descriptor
bin in pixels, and @e m is the <b>descriptor magnification factor</b>
which expresses how much larger a descriptor bin is compared to
the scale of the keypoint @f$ \sigma @f$
(the default value is @e m = 3). Moreover, the
standard SIFT descriptor computes the image gradient at the scale of
the keypoints, which in the canonical frame is equivalent to a
smoothing of @f$ \hat \sigma = 1/m @f$. Finally, the default
Gaussian window size is set to have standard deviation
 @f$ \hat \sigma_\mathrm{win} = 2 @f$. This yields the formula

@f{eqnarray*}
 h(t,i,j)
 &=&
 m \sigma \int
 g_{\sigma_\mathrm{win}}(\mathbf{x} - T)\,
 w_\mathrm{ang}(\angle J(\mathbf{x}) - \theta  - \theta_t)\,
 w_{ij}\left(\frac{R(\theta)^\top \mathbf{x} - T}{m\sigma}\right)\,
 |J(\mathbf{x})|\,
 d\mathbf{x},
\\
\sigma_{\mathrm{win}} &=& m\sigma\hat \sigma_{\mathrm{win}},
\\
 J(\mathbf{x})
 &=& \nabla (g_{m \sigma \hat \sigma} * I)(\mathbf{x})
 = \nabla (g_{\sigma} * I)(\mathbf{x})
 = \nabla I_{\sigma} (\mathbf{x}).
@f}

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsection framedet-tech-ellipse-detector Ellipse Detector
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsection framedet-tech-ellipse-frm Scale and affine invariant keypoint
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@image html affine-frame.png "Affine invariant keypoint."

Affine invariant keypoint (frame) consists from from the following values:

\f[
  f_K = \left(
  \begin{array}{ccccc}
    x_K & y_K & e_{11} & e_{12} & e_{22}
  \end{array} \right)^T
\f]

Affine shape is described by symetric shape matrix \f$ E^{-1} \f$
where \f$ \mathbf{x}^T E \mathbf{x} = 1 \f$ is the equation of the elllipse
matrix representation.
The shape matrix is given as:
\f[
  E^{-1} = \left(
  \begin{array}{cc}
    e_{11} & e_{12} \\
    e_{12} & e_{22}
  \end{array} \right)
\f]

Elipse is used in a similar manner as circles in scale invariant features
because it can represent the affine shape of the feature - it is affine
transformation of a circle. In fact this affine transformation \f$ A^{-1} \f$
(it is inverted matrix \f$ A \f$ because in the following text, \f$ A \f$ is
used as a transformation which transforms ellipse to circle) can be converted
into the ellipse matrix as:
\f[
  E = A^T A
\f]

Ellipse matrix, because it is real, symmetric and positive deffinite,
can be also decomposed using eigen decomposition:

\f[
  E = Q^{-1} \left(
  \begin{array}{cc}
    \lambda_1 & 0 \\
    0 & \lambda_2
  \end{array} \right) Q
\f]

The eigen-values can be relted to the size of ellipse axes where the
smaller eigen value \f$ \lambda_{max} \f$ with its eigen vector represents
the direction of the fastest change and the bigger eigen value \f$ \lambda_{max} \f$
with its eig. vector the direction of the slowest change.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsection framedet-tech-ellipse-affshape Affine shape estimation
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

In order to obtain frames which are covariant to affine transformation
we have to obtain its attributes which would be covariant to affine
transformation. These attributes can be also related as 'Affine shape'
of the image region because then the blob can be normalised into
isotropic blob using affine transformation encoded in these attributes
(which is an implication from the affine covariance).

The basic idea of affine shape estimation is to find an affine transformation
which would normalise the image blob. For a simplicity at first it would be
shown how to find a transformation which would normalise ellipse into a
circle.

@image html affine-transf.png "Affine shape estimation."

Anisotrpic image blob can be represented by an elipse with equation:
\f[
  \bar{\mathbf{x}}^T E \bar{\mathbf{x}} = 1
\f]

Where \f$ \bar{\mathbf{x}} = \mathbf{x} - \mathbf{x}_0 \f$ is a coordinate
system based on the kypoint centre \f$ \mathbf{x}_0 \f$ and \f$ \mathbf{x} \f$w
is coordinate in the image \f$ I(\mathbf{x}) \f$.

The affine transformation can be defined as a positive-definite matrix
\f$ A \f$ of size \f$ 2 \times 2 \f$ such that:
\f[
  \hat{\mathbf{x}} = A \bar{\mathbf{x}}
\f]

Where \f$ \hat{\mathbf{x}} \f$ are coordinates in the transformed image
\f$ \hat{I}(\hat{\mathbf{x}}) \f$.

We want to find affine transformation which would transform points on the
ellipse \f$ \mathbf{x} \f$ to points on a circle \f$ \hat{\mathbf{x}} \f$
s.t.
\f[
  \hat{\mathbf{x}}^T \hat{\mathbf{x}} = 1
\f]

This can be done by decomposing matrix \f$ E \f$ to:
\f{eqnarray*}
  E = E^{\frac{T}{2}} R^T R E^{\frac{1}{2}}
\f}

Where \f$ R \f$ is an arbitrary rotation matrix. Then when substituted to the
ellipse equation we can see that the affine transformation can be expressed as:
\f[
  A = R E^{\frac{1}{2}}
\f]
Therefore it transforms points on the ellipse to points on a circle.

NOTE - this is wrong somewhere because we have \f$ A^T A = E^{-1} \f$ not E.
It should also be true from the ellipse equation comparing the scales.

In order to describe the anisotropic structure we need a descriptor
\f$ \mu: \mathbb{R}^2 \rightarrow \mathbb{R}^2 \f$which
covariant to affine transformations such that:
\f[
  \mu (\hat{I}(\hat{\mathbf{x}})) = \mu_0 = \mu (I( A^{-1} \hat{\mathbf{x}}))
\f]

It was observed that a windowed second moment matrix has got some interesting
properties which makes is suited for estimating local linear distortion of the
image blob.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsubsection framedet-tech-ellipse-smm Windowed Second moment matrix
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

To show that second moment matrix can be used as a descriptor of the
affine shape of local anisotrpic structure we have to extend the
classic definition to our scale space and circular Gaussian distribution
must be replaced by multivariate gaussian distribution in order to
be able to detect anisotropic structures in our basically isotropic
scale space.

Multivariate 2D gaussian distribution is given as:
\f[
  g(\mathbf{x}, \Sigma) =
  \frac{1}{2 \pi | \Sigma |^{1/2}}
   \exp \left( -\frac{1}{2} \mathbf{x}^T \Sigma^{-1} \mathbf{x} \right)
\f]
Where \f$ \Sigma \f$ is covariance matrix of size \f$ 2 \times 2 \f$.
In the 1D case \f$ \Sigma \f$ would be a real number, variance of the
distribution \f$ \sigma^2 \f$.

Extending our framework we can now define the second moment matrix (SMM):
\f[
  \mu(\mathbf{x}, \Sigma_I, \Sigma_D) =
  |\Sigma_D| g(\Sigma_I) *
  (\nabla L(\mathbf{x}, \Sigma_D) \nabla L(\mathbf{x}, \Sigma_D)^T)
\f]

Where we define two covariance matrices, \f$ \Sigma_D \f$ which determines
'derivation' (or local) gaussian kernel. This covariance matrix describes
what initial smoothing was used for the data on which the derivations were
calculated. In our scale space framework we can see that it would be in
relation with the scale of the current level.
The covariance matrix \f$ \Sigma_I \f$ define 'integration' gaussian kernel
which define the window function over which the SMM is calculated,
in our case the window is selected as Gaussian function.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsubsection framedet-tech-ellipse-trsmm Transformation of SMM
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

Now having extended the definition of the second moment matrix, we
would like to see how SMM behaves when the image was transformed
using some affine transformation \f$ A \f$. Let have coordinates
\f$ \bar{\mathbf{x}} \f$ in 'left' image \f$ I(\bar{\mathbf{x}}) \f$ and
a cordinate \f$ \hat{\mathbf{x}} \f$ in 'right' image
\f$ \hat{I}(\hat{\mathbf{x}}) \f$. The image coordinate are in relation
\f$ \hat{\mathbf{x}}  = A \bar{\mathbf{x}}\f$ and image intensity values
are in relation \f$ \hat{I}(\hat{\mathbf{x}}) = I(A^{-1} \hat{\mathbf{x}}) \f$

The SMM \f$ M \f$ computed in point \f$ \bar{\mathbf{x}} \f$ in the left
image is in relation with the SMM \f$ \hat{M} \f$ computed in point
\f$ \hat{\mathbf{x}} \f$ in the right image:
\f{eqnarray*}
  M = \mu(\mathbf{x},\Sigma_{I}, \Sigma_{D}) &=&
    A^T \mu(\hat{\mathbf{x}}, \hat{\Sigma}_{I}, \hat{\Sigma}_{D}) A \\
    &=&
    A^T \mu(\hat{\mathbf{x}}, A \Sigma_{I} A^T, A \Sigma_{D} A^T) A
    = A^T \hat{M} A
\f}
Where the gaussian kernels are in relation: \f$ \hat{\Sigma} = A \Sigma A^T \f$

This is an important property of the SMM which is demanded by previously
defined affine shape estimation framework. Therefore we can use the SMM
in a similar way as the Ellipse matrix and use it to find the affine
transformation which would normalise the anisotropic image blob.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsubsection framedet-tech-ellipse-covm Covariance matrices and the SMM
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

In the previous section we have not considered the shape of the covariance
matrices in the SMM computation. General variation gives us six degrees
of freedom. However in this algorithm these matrices are chosen to be
in direct proportion to the second moment matrix because it yields
quite agreeable properties:
\f[
  \hat{\Sigma} = A ~\Sigma~ A^T = R M^{\frac{1}{2}} ~\Sigma~ M^{\frac{T}{2}} R^T
\f]
And when we choose \f$ \Sigma = \sigma M^{-1} \f$ we get:
\f[
  \hat{\Sigma} = \sigma ~ R M^{\frac{1}{2}}
                 ~M^{-1}~ M^{\frac{T}{2}} R^T = \sigma I
\f]
Which complies with the fact that the isotropic structure in the right image
should be sampled with isotropic, circular, gaussian kernels. Therefore we
define integration scale \f$ \sigma_I \f$ and derivation scale \f$ \sigma_D \f$
s.t. \f$ \Sigma_I = \sigma_I M^{-1} \f$ and \f$ \Sigma_D = \sigma_D M^{-1} \f$.

This decision has got also intuitive explanation that the image blob is sampled
with a gaussian window which has the same shape as its affine region.


<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsubsection framedet-tech-ellipse-iter Iterative procedure
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

We have shown that SMM can be used as an approximation descriptor of the
affine shape of an image blob. However for selecting the
covariance kernels of the SMM computation we need to know the SMM.
That is why the first covariance matrices are estimated as circular
and then we use iterative procedure which in each step selects more precise
affine region which is defined by shape adaptation matrix \f$ U_k \f$.

The convergence criterion of the iterative procedure is based on simple
isotropy measure based on second moment matrix \f$ M \f$ eigen values
\f$ \lambda_{min}(M) \f$ and \f$ \lambda_{max}(M) \f$:
\f[
  \mathcal{Q} = \frac{\lambda_{min}(M)}{\lambda_{max}(M)}
\f]
Where the \f$ \mathcal{Q} \in \langle 0, 1 \rangle \f$ with \f$ 1 \f$ for
perfect isotropic structure when the SMM is close to pure rotation. Then
the convergence criterion is defined as:
\f[
  1 - \mathcal{Q} < \epsilon_C
\f]
Where usually \f$ \epsilon_C = 0.05 \f$.

The iterative procedure of estimation of affine shape of point \f$ \mathbf{x} \f$
goes as follows:

Input: Scale invariant frame \f$ K \f$ with spatial coordinate
\f$ \mathbf{x}_K = (x_k, y_k) \f$ and scale \f$ \sigma_K \f$

-# initialize \f$ U_0 \f$ to identity matrix and set \f$ k = 0 \f$
-# set \f$ k = k+1 \f$ and if \f$ k \geq k_{max} \f$ reject the frame
   \f$ K \f$ as unstable.
-# normalize the actual window \f$ W(\hat{\mathbf{x}}) = I(\mathbf{x}) \f$
   centered in \f$ \sigma_K U_{k-1} \hat{\mathbf{x}}_K = \mathbf{x}_K \f$.
   Window \f$ W \f$ now contains affine region normalised into a circular
   one based on the actual estimation of the affine shape \f$ U_{k-1} \f$
   and the detected derivation scale \f$ \sigma_K \f$ of the frame
   (because \f$ |U_{k-1}| = 1 \f$).
   The window size is defined by argument window_size.
-# calculate the second moment matrix of the normalised window \f$ W \f$
   \f[
    M_k = \mu_W(\hat{\mathbf{c}}, \sigma_I I, \sigma_K I)
   \f]
   This is easily achieved by weighting the window values by circular
   gaussian window which standard deviation is
   \f$ \sigma_W = \frac{0.5 \mathrm{win\_size}}{3} \f$
   based on the empirical 3-sigma rule. We can see that because the frame
   affine regions is always fitted to the window of same size the integration
   scale is a constant multiply of the frame derivation scale.

-# concatenate transformation
   \f[
   U_K = \frac{M_k^{1/2}}{\left|M_k^{1/2}\right|} ~ U_{k-1}
   \f]
   Where the matrix \f$ M_k^{1/2} \f$ is normalised by its determinant
   in order not to scale the affine region size and only change its
   shape. In this way the actual affine transformation is iteratively
   improved as the measurement region gets more adapted to the real
   affine shape of the image structure.
-# In the case of divergence:
   \f[
    \frac{\lambda_{min}(M_k^{1/2})}{\lambda_{min}(M_k^{1/2})} > 6
   \f]
   e.g. when the point \f$ K \f$ lies on an edge, reject the point as
   unstable.
-# go to Step 2 if
  \f[
  1 - \frac{\lambda_{min}(M_k^{1/2})}{\lambda_{min}(M_k^{1/2})} \geq \epsilon_C
  \f] where the default value is \f$ \epsilon_C = 0.05 \f$. Please note that
  it is calculated from the SMM of the estimated normalised isotropic shape
  whereas the divergence criterion is calculated from the estimated affine
  transformation.
  If not, the estimated affine transformation is equal \f$ A = U_k \f$.

In the current implementation the derivation kernel used in the windowed SMM
is the kernel used in the GSS level where the point was detected, contrary to
the definition, because of the performance issues. This also fixes the
scale for which the frame was detected.

The affine transformation can be transformed to the shape matrix which also
contains the scale of the frame:
\f[
  E^{-1} = \sigma_{K}^{2} (A^T A)
\f]

In the case when oriented ellipse frames are detected, the affine image
blob is firstly normalised into a larger window. This normalisation is also
performed in more precise manner from the original image rather than from
the scale space level which is usually blurred more than it is neccessary.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsection framedet-tech-ellipse-norm Affine shape normalisation
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

In order to calculate descriptor of the affine invariant frame we have to
normalise the associated image blob according to its shape in order to be
able to calculate affine invariant frame descriptor. Generally we normalise
frame neighbourhood multiplied by measurement scale (parameter mr_scale)
into a square patch \f$P\f$ (normalised region is isotropic, therefore circular)
of size patch_size.

Because the scale space generally simulate downsampling of the input image and
is smoothed with circular kernels it does not take into account the affine
shape of the frame. Therefore if we would use scale space values for
descriptor calculation, the data would be smoothed by wrong kernel which would
not reflect the affine shape of the point neighbourhood (where the distances
between pixels are not constant).

That is why the data from the original image are used. Because the
downsampling is done using bilinear interpolation which takes into account
only nearest pixels, the image has to be smoothed accordingly as described
in the followng procedure. This can be also viewed from the sampling theorem
point of view, where the smoothing is just a low-pass 2D filter suppressing
higher frequencies in order to prevent aliasing.

Input: Affine invariant frame \f$ K \f$ with spatial coordinate
\f$ \mathbf{x}_K \f$, scale of the original scale invariant frame
\f$ \sigma_K \f$ and its affine shape described by the affine transformation
\f$ A_K, |A_K| = 1 \f$ and shape matrix \f$ E_K \f$.

-# Test if the frame measurement region touches the image boundary,
   if so, it is not possible to calculate the frame descriptor and
   the kyepoint is rejected.
-# Calaculate the radius \f$ r_c \f$ of circumscribed circle of the
   anistropic image structure as:
   \f[
    r_{c_K} = \sigma_K \mathrm{mr\_scale}
   \f]
   and the scale between the patch and the circumscribed circle:
   \f[
    s_{p \rightarrow c} = \frac{2 ~ r_{c_K}}{\mathrm{patch\_sz}}
   \f]
-# If the \f$ s_{p \rightarrow c} > 0.4 \f$, i.e. that the distance between
   the patch pixels in image is bigger than \f$ 0.4 \f$ so the image must be
   smoothed in order to perform correct bilinear interpolation
   -# Warp the measurement region into a temporary patch of size
      \f$ 2 * r_{c_K} \f$ with the affine transformation \f$ A_K \f$ using
      bilinear interpolation. Image is not anyhow scaled (\f$ |A_K| = 1 \f$)
      so the pixel distance remains the same.
   -# Smooth the \f$ P_t \f$ with circular gaussian kernel with
      \f$ \sigma = 1.5 s_{p \rightarrow c} \f$. The multiplication factor
      \f$ 1.5 \f$ has been chosen due to properties of the bilinear
      transformation so that pixels value influence spread accordingly
      to its neighbourhood.
    -# Downsample the \f$ P_t \f$ to \f$ P \f$ using bilinear interpolation.
-# If the \f$ s_{p \rightarrow c} \leq 0.4 \f$ then the smoothing is not
   needed and generally the measurement region is oversampled.
   Warp the affine measurement region directly to the patch \f$ P \f$
   using bilinear interpolation.

**/

#include "covdet.h"
#include "mathop.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <float.h>

/* ---------------------------------------------------------------- */
/*                                             Local definitions    */
/* ---------------------------------------------------------------- */

/** @internal @brief Use bilinear interpolation to compute orientations */
#define VL_SIFT_BILINEAR_ORIENTATIONS 1

#define EXPN_SZ  256          /**< ::fast_expn table size @internal */
#define EXPN_MAX 25.0         /**< ::fast_expn table max  @internal */
double expn_tab [EXPN_SZ+1] ; /**< ::fast_expn table      @internal */

#define log2(x) (log(x)/VL_LOG_OF_2)

/** ------------------------------------------------------------------
 ** @internal
 ** Parameters needed by hessian response function.
 **
 ** These parameters are passed to function ::vl_scalespace_apply.
 **/
typedef struct _HessianParams{
  double sigma0; /**< smoothing of pyramid base. */
  double sigmak; /**< k-smoothing */
} HessianParams;

/* ---------------------------------------------------------------- */
/*                                            Exported variables    */
/* ---------------------------------------------------------------- */

const char* vlFrameNames[VL_FRAMETYPE_NUM] = {
  "INVALID",
  "disc",
  "oriented disc",
  "ellipse",
  "oriented ellipse"
} ;

VlEnumerator vlFrameTypes [VL_FRAMETYPE_NUM] = {
  {"disc" ,             (vl_index)VL_FRAMETYPE_DISC             },
  {"oriented-disc",     (vl_index)VL_FRAMETYPE_ORIENTED_DISC    },
  {"ellipse",           (vl_index)VL_FRAMETYPE_ELLIPSE          },
  {"oriented-ellipse",  (vl_index)VL_FRAMETYPE_ORIENTED_ELLIPSE },
  {0,                   0                              }
} ;

VlEnumerator vlCovdetMethods [VL_COVDET_METHOD_NUM] = {
  {"dog" ,    (vl_index)VL_COVDET_METHOD_DOG     },
  {"hessian", (vl_index)VL_COVDET_METHOD_HESSIAN },
  {0,         0                      }
} ;

/* ---------------------------------------------------------------- */
/*                                               Local functions    */
/* ---------------------------------------------------------------- */

/** @internal @brief Fast @f$exp(-x)@f$ approximation
 ** @param x argument.
 ** @return approximation of @f$exp(-x)@f$.
 **
 ** The argument @a x must be in the range <code>[0, ::EXPN_MAX]</code>.
 **/

VL_INLINE double
fast_expn (double x)
{
  double a,b,r ;
  int i ;
  /*assert(0 <= x && x <= EXPN_MAX) ;*/

  if (x > EXPN_MAX) return 0.0 ;

  x *= EXPN_SZ / EXPN_MAX ;
  i = vl_floor_d (x) ;
  r = x - i ;
  a = expn_tab [i    ] ;
  b = expn_tab [i + 1] ;
  return a + r * (b - a) ;
}

/** @internal @brief Initialize tables for ::fast_expn
 **/

VL_INLINE void
fast_expn_init ()
{
  int k  ;
  for(k = 0 ; k < EXPN_SZ + 1 ; ++ k) {
    expn_tab [k] = exp (- (double) k * (EXPN_MAX / EXPN_SZ)) ;
  }
}

/** ------------------------------------------------------------------
 ** @internal @brief Compute det Hessian response of an image.
 ** @param outputImage output imgae buffer.
 ** @param inputImage input image buffer.
 ** @param width input image width.
 ** @param height input image height.
 ** @param sigma current amount of smoothing.
 **/

VL_INLINE void
vl_hessian_response (float * outputImage,
                     float const * inputImage,
                     vl_size width,
                     vl_size height,
                     double sigma)
{
  float norm2   = (float)(sigma * sigma);
  int const xo  = 1;      /* x-stride */
  int const yo  = width;      /* y-stride */
  vl_size r, c;

  float p11, p12, p13, p21, p22, p23, p31, p32, p33;

  /* setup input pointer to be centered at 0,1 */
  float const *in = inputImage + yo;

  /* setup output pointer to be centered at 1,1 */
  float *out = outputImage + xo + yo;

  /* move 3x3 window and convolve */
  for (r = 1; r < height - 1; ++r)
  {
    /* fill in shift registers at the beginning of the row */
    p11 = in[-yo]; p12 = in[xo - yo];
    p21 = in[  0]; p22 = in[xo     ];
    p31 = in[+yo]; p32 = in[xo + yo];
  /* setup input pointer to (2,1) of the 3x3 square */
  in += 2;
    for (c = 1; c < width - 1; ++c)
    {
      float Lxx, Lyy, Lxy;
    /* fetch remaining values (last column) */
      p13 = in[-yo]; p23 = *in; p33 = in[+yo];

    /* Compute 3x3 Hessian values from pixel differences. */
    Lxx = (-p21 + 2*p22 - p23);
    Lyy = (-p12 + 2*p22 - p32);
    Lxy = ((p11 - p31 - p13 + p33)/4.0f);

    /* normalize and write out */
    *out = (Lxx * Lyy - Lxy * Lxy)*norm2;

      /* move window */
      p11=p12; p12=p13;
      p21=p22; p22=p23;
      p31=p32; p32=p33;

    /* move input/output pointers */
      in++; out++;
    }
    out += 2;
  }

  /* Copy the computed values to borders */
  in  = outputImage + yo + xo ;
  out = outputImage + xo ;

  /* Top row without corners */
  memcpy(out, in, (width - 2)*sizeof(float));
  out--;
  in -= yo;

  /* Left border columns without last row */
  for (r = 0; r < height - 1; r++){
    *out = *in;
    *(out + yo - 1) = *(in + yo - 3);

    in += yo;
    out += yo;
  }

  /* Bottom corners */
  in -= yo;
  *out = *in;
  *(out + yo - 1) = *(in + yo - 3);

  /* Bottom row without corners */
  out++;
  memcpy(out, in, (width - 2)*sizeof(float));

}

/** ------------------------------------------------------------------
 ** @internal @brief Hessian operator callback
 **/

VL_INLINE void
vl_hessian_response_callback (float const *src_image,
                              int src_width, int src_height,
                              float *dst_image,
                              int VL_UNUSED(dst_w), int VL_UNUSED(dst_h),
                              int VL_UNUSED(octave), int level,
                              void* params)
{
  HessianParams* hp = (HessianParams*) params;
  double sigma = hp->sigma0 * pow(hp->sigmak, (double)level);

  vl_hessian_response(dst_image, src_image, src_width, src_height, sigma*sigma);
}

/* ---------------------------------------------------------------- */
/*                                            Image frames detector */
/* ---------------------------------------------------------------- */

/** ------------------------------------------------------------------
 ** @brief Create an isotropic frames (discs) detector
 ** @param imageWidth     image width.
 ** @param imageHeight    image height.
 ** @param respFunction   image response function used for frame detection
 ** @param numOctaves     number of octaves.
 ** @param numLevels      number of levels per octave.
 ** @param firstOctave    first octave index.
 ** @param calcOrient     do calc frame orientation
 ** @param calcDescr      do calc frame SIFT descriptor
 **
 ** The function allocates and returns a new Covariant frames detector for the
 ** specified image and scale space geometry.
 **
 ** When @a calcOrient is true, detector detects @a VlFrameOrientedDisc
 ** otherwise it detects @a VlFrameDisc. Each type has got its own get
 ** methods.
 **
 ** Setting @a numOctaves to a negative value sets the number of octaves
 ** to the maximum possible value depending on the size of the image.
 **
 ** @return the new disc frame Covariant frames detector.
 ** @sa ::vl_covdet_delete().
 **/

VL_EXPORT VlCovDet *
vl_covdet_new_disc_detector (vl_size imageWidth,
                             vl_size imageHeight,
                             VlCovDetMethod method,
                             vl_size numOctaves,
                             vl_index firstOctave,
                             vl_size numLevels,
                             vl_bool calcOrient,
                             vl_bool calcDescr)
{
  VlCovDet* self  = vl_calloc (1, sizeof(VlCovDet));
  if (self == NULL) goto err_alloc_self ;

  self->calcDescriptors = calcDescr ;
  self->calcOrientation = calcOrient ;
  self->method = method ;
  self->resp = vl_scalespace_new (imageWidth, imageHeight,
                                  numOctaves, firstOctave,
                                  numLevels, -1, numLevels) ;
  if (self == NULL) goto err_alloc_resp ;

  if (calcOrient) {
    self->framesType = VL_FRAMETYPE_ORIENTED_DISC ;
  } else {
    self->framesType = VL_FRAMETYPE_DISC ;
  }
  self->frameSize = vl_frame_size(self->framesType) ;

  if (method == VL_COVDET_METHOD_DOG) {
    /*
     * Initialise the GSS with one more scale level per octave
     * in order to have numLevels levels in the DoG.
     */
    self->gss = vl_scalespace_new (imageWidth, imageHeight,
                                   numOctaves, firstOctave,
                                   numLevels, -1, numLevels + 1) ;
    /* TODO: handle error */

    /*
     * The gaussian sc. space should have the same properties as DoG,
     * the one plane more there is only in order to get requested
     * number of planes in DoG.
     */

    self->gss->sigmak = self->resp->sigmak ;
    self->gss->sigma0 = self->resp->sigma0 ;
    self->gss->dsigma0 = self->resp->dsigma0 ;
    self->gss->numLevels = self->resp->numLevels ;
    self->peakThreshold = 0.0 ;

  } else if (method == VL_COVDET_METHOD_HESSIAN) {
    self->gss = vl_scalespace_new (imageWidth, imageHeight,
                                   numOctaves, firstOctave,
                                   numLevels, -1, numLevels) ;
    self->peakThreshold = 16.0f/3.0f*16.0f/3.0f ;
  } else {
    VL_ASSERT(0,"Invalid response function.");
  }

  if (calcOrient || calcDescr) {
    /* Get GSS number of octaves - gradient has got different size so it can
     * yield different number of octaves.
     */
    numOctaves = vl_scalespace_get_octaves_num(self->gss);
    self->grad = vl_scalespace_new(imageWidth*2, imageHeight,
                                   numOctaves, firstOctave,
                                   numLevels, 0, numLevels - 1);
  } else {
    self-> grad = 0;
  }

  self-> affineShapeEstimator = 0;
  self-> affineNorm = 0;
  self-> patchGrad = 0;

  self-> frames = 0;
  self-> numFrames = 0;
  self-> framesRes = 0;
  self-> descriptors = 0;
  self-> descriptorsRes = 0;

  self-> normThreshold = 0.0 ;
  self-> magnification = 3.0 ;
  self-> windowSize = VL_SIFT_NBP / 2 ;
  self-> calcInverseSm = VL_FALSE;
  self-> edgeThreshold = 10.0 ;

  /* initialize fast_expn stuff */
  fast_expn_init () ;

  return self;
err_alloc_resp:
  vl_free(self) ;
err_alloc_self:
  return NULL ;
}

/** ------------------------------------------------------------------
 ** @brief Create a new ellipse frames detector
 ** @param imageWidth     image width.
 ** @param imageHeight    image height.
 ** @param respFunction   image response function used for frame detection
 ** @param numOctaves     number of octaves.
 ** @param firstOctave    first octave index.
 ** @param numLevels      number of levels per octave.
 ** @param affineWinSize  size of the window for affine shape detection
 ** @param descrWinSize   size of the patch for sift calculation
 ** @param calcOrient     do calc frames orientations.
 ** @param calcDesc       do calc frames descriptors
 **
 ** The function allocates and returns a new Covariant frames detector for the
 ** specified image and scale space geometry.
 **
 ** When @a calcOrient is true, detector detects @a VlFrameOrientedEllipse
 ** otherwise it detects @a VlFrameEllipse. Each type has got its own get
 ** methods.
 **
 ** Setting @a numOctaves to a negative value sets the number of octaves
 ** to the maximum possible value depending on the size of the image.
 **
 ** @return the new ellipse frame Covariant frames detector.
 ** @sa ::vl_covdet_delete().
 **/

VL_EXPORT VlCovDet *
vl_covdet_new_ellipse_detector (vl_size imageWidth, vl_size imageHeight,
                                VlCovDetMethod respFunction,
                                vl_size numOctaves, vl_index firstOctave, vl_size numLevels,
                                vl_size affineWindowSize, vl_size descrWindowSize,
                                vl_bool calcOrient, vl_bool calcDesc)
{
  VlCovDet *self = vl_covdet_new_disc_detector (imageWidth, imageHeight, respFunction,
                                                numOctaves, firstOctave, numLevels,
                                                VL_FALSE, VL_FALSE);

  self->calcDescriptors = calcDesc;
  self->calcOrientation = calcOrient;
  self->affineShapeEstimator = vl_affineshapeestimator_new (affineWindowSize);

  if (calcOrient) {
    self->framesType = VL_FRAMETYPE_ORIENTED_ELLIPSE ;
  } else {
    self->framesType = VL_FRAMETYPE_ELLIPSE ;
  }
  self->frameSize = vl_frame_size(self->framesType);

  if (calcDesc || calcOrient) {
    self-> affineNorm = vl_affinepatchnormalizer_new (imageWidth, imageHeight, descrWindowSize);
    self->patchSize = descrWindowSize ;
    self->patch = vl_malloc(descrWindowSize * descrWindowSize * sizeof(float)) ;
    self->patchGrad = vl_malloc(descrWindowSize * descrWindowSize * 2 * sizeof(float)) ;
    self->patchSigma = 1.5 ;
    self->magnification = (double)descrWindowSize / VL_SIFT_NBP ;
  }

  return self ;
}

/** -------------------------------------------------------------------
 ** @brief Delete a ::VlFrameDet object
 ** @param self object to delete.
 ** @sa ::vl_covdet_new_disc_detector()
 ** @sa ::vl_covdet_new_ellipse_detector()
 **/

VL_EXPORT void
vl_covdet_delete (VlCovDet* self)
{
  if (self) {
    if (self->gss) vl_scalespace_delete (self-> gss);
    if (self->resp) vl_scalespace_delete (self->resp);
    if (self->grad) vl_scalespace_delete (self-> grad);
    if (self->descriptors) vl_free (self->descriptors);
    if (self->frames) vl_free (self->frames);
    if (self->patch) vl_free(self->patch) ;
    if (self->patchGrad) vl_free (self-> patchGrad);
    if (self->affineShapeEstimator) vl_affineshapeestimator_delete (self-> affineShapeEstimator);
    if (self->affineNorm) vl_affinepatchnormalizer_delete (self-> affineNorm);
    vl_free(self);
  }
}

/** ------------------------------------------------------------------
 ** @internal @brief Store ::VlScaleSpaceFrame into the frames storage
 ** @param self Covariant frames detector object.
 ** @param ssFrame Pointer to scale space frame.
 ** @param angle Angle of the frame, used only for oriented discs.
 ** @param frame Pointer to frame storage.
 ** @return Pointer to the next frame in the frames storage.
 **
 ** Stores scale space frame to the frames storage based on the detector
 ** settings.
 **/

VL_INLINE void*
_vl_covdet_store_scalespace_frame(VlCovDet const *self,
                                    VlScaleSpaceFrame const *ssFrame,
                                    double angle, void *frame)
{
  VlFrameType frm_type = self->framesType;

  if (frm_type == VL_FRAMETYPE_DISC) {
    VlFrameDisc* disc_frm = (VlFrameDisc*) frame;
    disc_frm->x = ssFrame->x;
    disc_frm->y = ssFrame->y;
    disc_frm->sigma = ssFrame->sigma;
    frame = (void*)(disc_frm + 1);
  } else {
    VlFrameOrientedDisc* odisc_frm = (VlFrameOrientedDisc*) frame;
    odisc_frm->x = ssFrame->x;
    odisc_frm->y = ssFrame->y;
    odisc_frm->sigma = ssFrame->sigma;
    odisc_frm->angle = angle;
    frame = (void*)(odisc_frm + 1);
  }

  return frame;
}

/** ------------------------------------------------------------------
 ** @internal @brief Stores ::VlAffineShapeEstimatorFrame to frames storage.
 ** @param self Covariant frames detector object.
 ** @param affineFrame Pointer to affine detector frame.
 ** @param angle Angle of the frame, used only oriented ellipses.
 ** @param frame Pointer to frame storage.
 ** @return Pointer to the next frame in the frames storage.
 **
 ** Converts the frame into apropriate output frame based on the detector
 ** settings.
 **/

VL_INLINE void*
_vl_covdet_store_affinedet_frame (VlCovDet *self,
                                    VlAffineShapeEstimatorFrame const *affineFrame,
                                    double angle, void *frame)
{
  VlFrameType frm_type    = self->framesType;

  if (frm_type == VL_FRAMETYPE_ELLIPSE) {
    VlFrameEllipse* el_frm = (VlFrameEllipse*)frame;
    double a = affineFrame->a11 * affineFrame->sigma ;
    double b = affineFrame->a12 * affineFrame->sigma ;
    double c = affineFrame->a21 * affineFrame->sigma ;
    double d = affineFrame->a22 * affineFrame->sigma ;
    double e11 = a*a + b*b;
    double e12 = a*c + b*d;
    double e22 = c*c + d*d;

    el_frm->x = affineFrame->x;
    el_frm->y = affineFrame->y;
    if (!self->calcInverseSm) {
      el_frm->e11 = e11;
      el_frm->e12 = e12;
      el_frm->e22 = e22;
    } else {
      /* Calc the inverse of the shape matrix */
      double t, r, l1, l2;
      if (e12 != 0)
      {
        r = (double)(e22-e11)/(2.0 * e12);
        if (r>=0) t =  1.0/( r+sqrt(1+r*r));
         else     t = -1.0/(-r+sqrt(1+r*r));
        r = 1.0/sqrt(1+t*t); /* c */
        t = t*r;             /* s */
      } else {
        r = 1;
        t = 0;
      }

      l1 = 1./(r*r*e11-2*r*t*e12+t*t*e22);
      l2 = 1./(t*t*e11+2*r*t*e12+r*r*e22);
      el_frm->e11   =   (float) r*r*l1+t*t*l2;
      el_frm->e12   =   (float)-r*t*l1+t*r*l2;
      el_frm->e22   =   (float) t*t*l1+r*r*l2;
    }

    frame = (void*)(el_frm + 1);

  } else if (frm_type == VL_FRAMETYPE_ORIENTED_ELLIPSE) {
    VlFrameOrientedEllipse* oel_frm = (VlFrameOrientedEllipse*)frame;
    double sc = affineFrame->sigma;
    double a11 = affineFrame->a11 * sc;
    double a12 = affineFrame->a12 * sc;
    double a21 = affineFrame->a21 * sc;
    double a22 = affineFrame->a22 * sc;
    double sin_a, cos_a;
    /* Rotate by -angle in order to make the major angle=0 */
    sin_a = sin(angle);
    cos_a = cos(angle);

    oel_frm->x = affineFrame->x;
    oel_frm->y = affineFrame->y;
    /* A is affine transformation from ellipse to circle,
     * then we have to first transform and then rotate,
     * therefore A_n = R' * A
     */
    oel_frm->a11 =  a11*cos_a + a21*sin_a;
    oel_frm->a12 =  a12*cos_a + a22*sin_a;
    oel_frm->a21 = -a11*sin_a + a21*cos_a;
    oel_frm->a22 = -a12*sin_a + a22*cos_a;

    frame = (void*)(oel_frm + 1);
  }
  return frame;
}

/** ------------------------------------------------------------------
 ** @internal @brief Compute disc frames based on the detected scale space frames.
 **
 ** Computes orientations and descriptors if needed and exports
 ** the results into the frames and descriptors storage.
 **
 ** @param self Covariant frames detector object.
 **/

VL_INLINE void
_vl_covdet_discs_detection(VlCovDet* self)
{
  VlScaleSpace *gss       = self->gss;
  VlScaleSpace *resp      = self->resp;
  VlScaleSpace *grad      = self->grad;
  vl_bool calc_desc       = self->calcDescriptors;
  vl_bool calc_orient     = self->calcOrientation;
  int desc_size           = vl_covdet_get_descriptor_size(self);
  float *desc             = self->descriptors;
  int scsp_kpts_num       =  vl_scalespace_get_frames_num(resp);
  VlScaleSpaceFrame const *ssfrm = vl_scalespace_get_frames(resp);
  void* frames            = self->frames;
  int kpts_num            = 0;
  int i;

  /* Compute gradient when needed */
  if (calc_desc || calc_orient) {
    vl_scalespace_apply(gss, grad, vl_imgradient_polar_f_callback, 0);
  }

  for (i = 0; i < scsp_kpts_num; ++i)
  {
    double angles[VL_COVDET_MAX_ANGLES];
    int a, angles_num;
    angles_num = 1;
    angles[0] = 0.;

    if (calc_orient) {
      angles_num = vl_covdet_calc_ssframe_orientations(grad, angles, ssfrm);
    }

    for (a = 0; a < angles_num; ++a) {
      frames =
          _vl_covdet_store_scalespace_frame(self, ssfrm, angles[a], frames);

      if (calc_desc) {
        vl_covdet_calc_ssframe_descriptor(self, desc, ssfrm, angles[a]);
        desc += desc_size;
      }
      kpts_num++;
    }
    ssfrm++;
  }

  self->numFrames = kpts_num;
}


/** ------------------------------------------------------------------
 ** @internal @brief Compute ellipse frames from the scale space frames
 ** @param self Covariant frames detector object
 ** @param image Pointer to the original image
 **
 ** Compute the affine shape of the scale space frame and compute its
 ** orientations and/or descriptors eventually based on the detector
 ** settings. Export the frames into frames storage.
 **/

VL_INLINE void
_vl_covdet_ellipses_detection (VlCovDet *self)
{
  float *descriptor = self->descriptors;
  vl_size descriptorSize = vl_covdet_get_descriptor_size(self) ;
  VlScaleSpaceFrame const *scaleSpaceFrame = vl_scalespace_get_frames(self->resp) ;
  VlScaleSpaceFrame const *scaleSpaceFrameEnd = scaleSpaceFrame + vl_scalespace_get_frames_num(self->resp) ;
  void *frame = self->frames ;
  vl_size numFrames = 0 ;

  for ( ; scaleSpaceFrame < scaleSpaceFrameEnd ; scaleSpaceFrame++) {
    VlAffineShapeEstimatorFrame affineFrame ;
    vl_bool hasAffineShape ;
    double angles [4] = {0.0, 0.0, 0.0, 0.0} ;
    int numAngles = 1 ;
    int a ;
    double factor = 2.0 * self->magnification / self->patchSize ;

    /*
     patchSize pixels are mapped in the rectangle [-2,2] representing
     containing the unit circle (standard reference frame).
     */

    hasAffineShape = vl_affineshapeestimator_estimate (self->affineShapeEstimator,
                                                       self->gss,
                                                       scaleSpaceFrame,
                                                       &affineFrame) ;
    if (! hasAffineShape) continue ;
    if (self->calcOrientation || self->calcDescriptors) {
      int error ;
      error = vl_scalespace_extract_affine_patch
      (self->gss,
       self->patch, self->patchSize, self->patchSize, self->patchSigma,
       affineFrame.x,
       affineFrame.y,
       factor * affineFrame.sigma * affineFrame.a11,
       factor * affineFrame.sigma * affineFrame.a21,
       factor * affineFrame.sigma * affineFrame.a12,
       factor * affineFrame.sigma * affineFrame.a22) ;

      vl_imgradient_polar_f (self->patchGrad, self->patchGrad + 1,
                           2, self->patchSize * 2,
                           self->patch,
                           self->patchSize, /* width */
                           self->patchSize, /* hieght */
                           self->patchSize) ; /* stride */

      if (self->calcOrientation) {
        numAngles = vl_calc_frame_orientations (self->patchGrad,
                                                self->patchSize,
                                                self->patchSize,
                                                angles,
                                                self->patchSize / 2.0,
                                                self->patchSize / 2.0,
                                                self->patchSize / 8.0,
                                                0);
      }
    }
#if 0
    {
      char name [256] ;
      snprintf(name, 256, "/Users/vedaldi/Desktop/bla/%10.0f%10.0f.pgm",
               affineFrame.x, affineFrame.y) ;
      vl_pgm_write_f(name, self->patch, self->patchSize, self->patchSize) ;
    }
#endif

    for (a = 0 ; a < numAngles ; ++ a) {
      frame = _vl_covdet_store_affinedet_frame (self, &affineFrame,
                                                  angles[a], frame) ;
      ++ numFrames ;

      if (self->calcDescriptors) {
        /* set window size so it the whole patch is sampled */
        const double windowSize = (double)VL_SIFT_NBP/2.;
        vl_sift_calc_descriptor (self->patchGrad, self->patchSize, self->patchSize,
                                 descriptor,
                                 self->patchSize / 2.0, self->patchSize / 2.0,
                                 1.0, angles[a],
                                 self->magnification, windowSize,
                                 self->normThreshold, 0) ;
        descriptor += descriptorSize;
      }
    }
  }
  self->numFrames = numFrames ;
}

/** ------------------------------------------------------------------
 ** @internal @brief Allocate enough space for frames and descriptors.
 **
 ** Allocates frames and descriptors storage based on the demanded type
 ** of frames.
 **
 ** @param self Covariant frames detector object.
 ** @param numFrames Number of frames
 **/

VL_INLINE void
_vl_covdet_alloc_frame_storage(VlCovDet* self, vl_size numFrames)
{
  vl_size kpt_size  = self->frameSize;
  vl_bool calc_desc = self->calcDescriptors;
  vl_size need_kpts_res, need_descs_res;


  need_kpts_res = kpt_size * numFrames;
  if(self-> framesRes < need_kpts_res) {
    if (self-> frames)
      vl_free(self->frames);
    self->frames = vl_malloc(need_kpts_res);
    self->framesRes = need_kpts_res;
  }

  if (calc_desc) {
    vl_size desc_size = vl_covdet_get_descriptor_size(self);
    need_descs_res = sizeof(float) * desc_size * numFrames;
    if(self-> descriptorsRes < need_descs_res) {
      if (self->descriptors)
        vl_free(self->descriptors);
      self->descriptors = vl_malloc(need_descs_res);
      self->descriptorsRes = need_descs_res;
    }
  }
}

/** -------------------------------------------------------------------
 ** @brief Detect frames in an image using defined Covariant frames detector
 **
 ** @param self Covariant frames detector
 ** @param image Pointer to image of size defined in @a self constructor
 **
 ** This functions detects frames in input image and if defined in the
 ** constructor also their SIFT descriptors.
 ** Detects frames in the @a image using the configured @a VlFrameDet
 ** object. If defined in the constructor of the object, also computes
 ** descriptors of the detected frames.
 **
 ** The results are accessible with functions depending on the detector
 ** definition in the constructor:
 ** - ::vl_covdet_get_discs() for discs ::VlFrameDisc when
 **     constructed with
 **     @code
 **       vl_covdet_new_disc_detector(..., calc_orient=VL_FALSE, ...)
 **     \endcode
 ** - ::vl_covdet_get_oriented_discs() for oriented discs ::VlFrameDisc when
 **     constructed with
 **     @code
 **       vl_covdet_new_disc_detector(..., calc_orient=VL_TRUE, ...)
 **     \endcode
 ** - ::vl_covdet_get_ellipses() for ellipses ::VlFrameEllipse when
 **     constructed with
 **     @code
 **       vl_covdet_new_ellipse_detector(..., calc_orient=VL_FALSE, ...)
 **     \endcode
 ** - ::vl_covdet_get_oriented_ellipses() for oriented ellipses
 **   ::VlFrameOrientedEllipse when constructed with
 **     @code
 **       vl_covdet_new_ellipse_detector(..., calc_orient=VL_TRUE, ...)
 **     \endcode
 ** and ::vl_covdet_get_frames_num() to get number of detected frames.
 **
 ** For accessing the SIFT descriptors of the frames use
 ** ::vl_covdet_get_descriptors(), ::vl_covdet_get_descriptors_num() and
 ** ::vl_covdet_get_descriptor_size().
 **/

VL_EXPORT void
vl_covdet_detect (VlCovDet *self, float const *image)
{
  vl_size maxNumFrameOrientations = self->calcOrientation ? VL_COVDET_MAX_ANGLES : 1;
  self->numFrames = 0;

  vl_scalespace_init(self->gss, image);

  switch (self->method) {
    case VL_COVDET_METHOD_DOG:
      vl_scalespace_diff (self->gss, self->resp) ;
      break ;

    case VL_COVDET_METHOD_HESSIAN: {
      HessianParams params;
      params.sigma0 = vl_scalespace_get_sigma0(self->gss) ;
      params.sigmak = vl_scalespace_get_sigmak(self->gss) ;
      vl_scalespace_apply (self->gss, self->resp,
                           &vl_hessian_response_callback,
                           &params);
      } break ;

    default:
      assert(0) ;
  }

  vl_scalespace_find_local_extrema (self->resp, self->peakThreshold, 1) ;
  vl_scalespace_refine_local_extrema (self->resp, self->peakThreshold, self->edgeThreshold, 1) ;
  self->numFrames = vl_scalespace_get_frames_num (self->resp) ;

  /*
   * Size of the memory needed is rather overestimated because not all
   * frames has got all 4 major angles when frame orientation
   * is demanded.
   */
  _vl_covdet_alloc_frame_storage(self, self->numFrames * maxNumFrameOrientations);

  if (self->framesType == VL_FRAMETYPE_DISC || self->framesType == VL_FRAMETYPE_ORIENTED_DISC) {
    _vl_covdet_discs_detection (self) ;
  } else {
    _vl_covdet_ellipses_detection (self) ;
  }
}

/* ---------------------------------------------------------------- */
/*                                                 Frame conversion */
/* ---------------------------------------------------------------- */

/** ------------------------------------------------------------------
 ** @internal
 ** @brief Convert a frame into a VlAffinDetFrame based on its type.
 **
 ** Covert any output frame into VlAffineShapeEstimator frame which is needed for
 ** orientation assingment or descriptor calculation of elliptic frames.
 **
 ** If needed the affine shape of the frame neighbourhood is calculated.
 **
 ** @param self Covariant frames detector object.
 ** @param dstFrame Pointer to the output affine frame
 ** @param srcFramePtr Pointer to the input frame
 ** @param srcFrameType Type of the input frame
 **
 ** @returns VL_TRUE if the frame was succesfully converted.
 **/

VL_INLINE vl_bool
_vl_covdet_convert_to_affinedet_frame(VlCovDet *self,
                                        VlAffineShapeEstimatorFrame *dstFrame,
                                        void const *srcFramePtr,
                                        VlFrameType srcFrameType)
{
  VlAffineShapeEstimator *aff_det = self->affineShapeEstimator;
  VlScaleSpace *gss = self->gss;

  /* Initialise VlAffineShapeEstimator frame from points without affine shape */
  if (srcFrameType == VL_FRAMETYPE_DISC || srcFrameType == VL_FRAMETYPE_ORIENTED_DISC){
    VlScaleSpaceFrame ssfrm;
    double x, y, sigma;
	vl_bool has_aff_shape ;

    if (srcFrameType == VL_FRAMETYPE_DISC) {
      VlFrameDisc const *frm = (VlFrameDisc const *)srcFramePtr;
      x = frm->x;
      y = frm->y;
      sigma = frm->sigma;
      srcFramePtr = (void*)(frm + 1);
    } else {
      /* In case of oriented disc, the orientation is irrelevant
       * because is calculated over isotropic region.
       */
      VlFrameOrientedDisc const *frm =
          (VlFrameOrientedDisc const *)srcFramePtr;
      x = frm->x;
      y = frm->y;
      sigma = frm->sigma;
    }

    vl_scalespace_frame_init(gss, &ssfrm, x, y, sigma);
    /* Detect the affine shape */
    has_aff_shape =
        vl_affineshapeestimator_estimate(aff_det, gss, &ssfrm, dstFrame);
    if(!has_aff_shape) {
      return VL_FALSE;
    }
  } else {
    if (srcFrameType == VL_FRAMETYPE_ELLIPSE) {
      /* If the input frame is ellipse we need to find the affine
       * transofrmation.
       */
      VlFrameEllipse const *frm = (VlFrameEllipse const *)srcFramePtr;
      vl_affineshapeestimator_frame_init_from_ell(aff_det, dstFrame,frm->x, frm->y,
                                       frm->e11, frm->e12, frm->e22);
    } else {
      /* Otherwise initialise from the affine transformation.
       * Please note that the angle is part of the affine transform.
       * so the patch will be already normalised by rotation.
       */
      VlFrameOrientedEllipse const *frm =
          (VlFrameOrientedEllipse const *)srcFramePtr;
      vl_affineshapeestimator_frame_init_from_aff(aff_det, dstFrame,frm->x, frm->y,
                                       frm->a11, frm->a12,
                                       frm->a21, frm->a22);
    }
  }
  return VL_TRUE;
}

/** ------------------------------------------------------------------
 ** @internal @brief Convert a frame into a VlScaleSpaceFrame.
 ** @param self Covariant frames detector object.
 ** @param dstFrame Pointer to the output ScaleSpace frame
 ** @param srcFrame Pointer to the input frame
 ** @param srcFrameType Type of the input frame
 ** @return Angle if the input frame was oriented disc, 0 otherwise.
 **
 ** If the frame type has got affine shape, only the scale is extracted.
 ** In case of the oriented elliptic frame, its orientation is irrelevant
 ** for an isotropic frame and therefore it is ignored.
 **/

VL_INLINE float
_vl_covdet_convert_to_scalespace_frame(VlCovDet* self,
                                         VlScaleSpaceFrame* dstFrame,
                                         void const * srcFrame,
                                         VlFrameType srcFrameType)
{
  double x = 0, y = 0, sigma = 0, angle = 0 ;

  switch (srcFrameType) {
  case VL_FRAMETYPE_DISC: {
    VlFrameDisc const *frm = (VlFrameDisc const *)srcFrame;
    x = frm->x;
    y = frm->y;
    sigma = frm->sigma;
  } break;
  case VL_FRAMETYPE_ORIENTED_DISC: {
    VlFrameOrientedDisc const *frm = (VlFrameOrientedDisc const *)srcFrame;
    x = frm->x;
    y = frm->y;
    sigma = frm->sigma;
    angle = frm->angle;
  } break;
  case VL_FRAMETYPE_ELLIPSE: {
    VlFrameEllipse const *frm = (VlFrameEllipse const *)srcFrame;
    double det;

    x = frm->x;
    y = frm->y;

    det = frm->e11 * frm->e22 - frm->e12*frm->e12;
    /* Because det(a A) = a^n det(A) and E = A'A */
    sigma = sqrt(sqrt(det));
  } break;
  case VL_FRAMETYPE_ORIENTED_ELLIPSE: {
    VlFrameOrientedEllipse const *frm =
        (VlFrameOrientedEllipse const *)srcFrame;
    double det;

    x = frm->x;
    y = frm->y;

    det = frm->a11 * frm->a22 - frm->a12 * frm->a21;
    /* Because det(a A) = a^n det(A) */
    sigma = sqrt(det);
  } break;
  default:
    VL_ASSERT(0, "Invalid frame type.");
    break;
  }

  vl_scalespace_frame_init(self->gss, dstFrame, x, y, sigma);
  return angle;
}

/** ------------------------------------------------------------------
 ** @brief Convert covariant image frames
 ** @param self Covariant frames detector object
 ** @param image Input image
 ** @param srcFrames Pointer to the input frames storage
 ** @param srcFramesNum Number of input frames
 ** @param srcFramesType Type of the input frames
 **
 ** Convert array of frames into frames of type defined in the
 ** ::VlFrameDet object constructor.
 **
 ** Converted frames are stored in the same storage as detected frames
 ** and are accessible using the same methods as for accessing detected
 ** frames
 **
 ** @sa ::vl_covdet_new_disc_detector(), ::vl_covdet_new_ellipse_detector()
 ** @sa ::vl_covdet_get_frames_num(), ::vl_covdet_get_frames_type()
 **/

VL_EXPORT void
vl_covdet_convert_frames(VlCovDet *self, float const *image,
                         void const *srcFrames, vl_size srcFramesNum,
                         VlFrameType srcFramesType)
{

  /* TODO the conversion of oriented-ellipse does not work properly in case
   * of descriptor calculation because calculating the descriptor with an
   * angle and calculating the descriptor on rotated patch yields different
   * results...
   */

  VlScaleSpace *gss = self->gss;
  VlScaleSpace *grad = self->grad;
  VlAffineShapeEstimator *aff_det = self->affineShapeEstimator;
  VlAffinePatchNormalizer *aff_norm = self->affineNorm;
  vl_size desc_sz = vl_covdet_get_descriptor_size(self);

  VlFrameType out_frms_type = self->framesType;
  vl_bool in_has_orientation = (srcFramesType == VL_FRAMETYPE_ORIENTED_DISC ||
                                srcFramesType == VL_FRAMETYPE_ORIENTED_ELLIPSE);
  vl_bool out_has_affine = (out_frms_type == VL_FRAMETYPE_ELLIPSE ||
                         out_frms_type == VL_FRAMETYPE_ORIENTED_ELLIPSE);
  vl_bool out_has_orient = self->calcOrientation;
  vl_bool out_calc_descs = self->calcDescriptors;

  vl_size in_frm_sz = vl_frame_size(srcFramesType);
  /* When oriented keypoint type is changed, calc orientation again. */
  vl_bool in_calc_orient = out_has_orient && (srcFramesType != out_frms_type);
  /* Discard orientation of input frames always when it is recalculated
   * or when the output ignores it. */
  vl_bool in_discard_orient = in_has_orientation &&
                              (in_calc_orient || !out_has_orient);
  int out_frms_sz_mult = in_calc_orient ? VL_COVDET_MAX_ANGLES : 1;

  void *out_frms_ptr = NULL;
  void const *in_frms_ptr = NULL;
  float *out_descrs_ptr = NULL;
  int nframes = 0;
  unsigned int i;

  _vl_covdet_alloc_frame_storage(self, srcFramesNum * out_frms_sz_mult);
  out_frms_ptr = self->frames;
  out_descrs_ptr = self->descriptors;

  in_frms_ptr = srcFrames;

  if((srcFramesType == out_frms_type) && !out_calc_descs) {
    /* There is actually nothing to do, the frames are the same... */
    self->numFrames = srcFramesNum;
    memcpy(out_frms_ptr, in_frms_ptr, srcFramesNum * in_frm_sz);
    return;
  }

  vl_scalespace_init(gss, image);

  /* Ellipses .......................................... */
  if (out_has_affine) {
    float* patch = 0;
    float* patch_grad = 0;
    int patch_sz = 0;
    double pt_pos = 0.;
    VlAffineShapeEstimatorFrame last_aff_frm;
    double factor = 2.0 * self->magnification / self->patchSize ;

    if (in_discard_orient){
      vl_affineshapeestimator_frame_init_from_aff(aff_det, &last_aff_frm,
                                       0., 0., 0., 0., 0., 0.);
    }

    if (in_calc_orient || out_calc_descs) {
      patch = vl_affinepatchnormalizer_get_patch(aff_norm);
      patch_sz = vl_affinepatchnormalizer_get_patch_size(aff_norm);
      pt_pos = (double)patch_sz/2.;
    }

    for (i = 0; i < srcFramesNum; ++i)
    {
      double angles[4];
      int nangles, a;
      VlAffineShapeEstimatorFrame aff_frm;
      angles[0] = 0.;
      nangles = 1;

/* Macro for incrementing void pointer */
#define charp(x) (*(char**)&x)

      if (!_vl_covdet_convert_to_affinedet_frame(self, &aff_frm, in_frms_ptr,
                                                   srcFramesType)) {
        charp(in_frms_ptr) += in_frm_sz;
        continue;
      }

      if(in_discard_orient) {
        /* Discard ellipses with the same shape as the previous ellipse. */
        if (aff_frm.x == last_aff_frm.x &&
            aff_frm.y == last_aff_frm.y &&
            aff_frm.a11 == last_aff_frm.a11 &&
            aff_frm.a12 == last_aff_frm.a12 &&
            aff_frm.a21 == last_aff_frm.a21 &&
            aff_frm.a22 == last_aff_frm.a22 &&
            aff_frm.sigma == last_aff_frm.sigma) {
          charp(in_frms_ptr) += in_frm_sz;
          continue;
        } else {
          last_aff_frm = aff_frm;
        }
      }

      /* Normalise the patch and calc gradient when needed */
      if (in_calc_orient || out_calc_descs) {
        int err = vl_scalespace_extract_affine_patch(
                    gss, self->patch, self->patchSize, self->patchSize,
                    self->patchSigma,
                    aff_frm.x, aff_frm.y,
                    factor * aff_frm.sigma * aff_frm.a11,
                    factor * aff_frm.sigma * aff_frm.a21,
                    factor * aff_frm.sigma * aff_frm.a12,
                    factor * aff_frm.sigma * aff_frm.a22);
        if (err == VL_ERR_OK) {
          patch_grad = self->patchGrad;

          vl_imgradient_polar_f(patch_grad, patch_grad + 1, 2, patch_sz*2,
                              patch, patch_sz, patch_sz, patch_sz);

          /*
           * If not, angle is part of the affine transformation and is not
           * needed to be recalculated.
           */
          if (in_calc_orient){
            nangles = vl_calc_frame_orientations(patch_grad,
                                                 patch_sz, patch_sz,
                                                 angles, pt_pos, pt_pos,
                                                 patch_sz/8., 0);
          }

        } else {
          charp(in_frms_ptr) += in_frm_sz;
          continue;
        }
      }

      /* Store supplemented frames to storage for export to matlab ... */
      for(a = 0; a < nangles; ++a) {
        nframes++;

        if (out_frms_type == srcFramesType) {
          /* When the keypoint type is the same, simply copy the data */
          memcpy(out_frms_ptr, in_frms_ptr, in_frm_sz);
          charp(out_frms_ptr) += in_frm_sz;
        } else {
          out_frms_ptr =
              _vl_covdet_store_affinedet_frame(self, &aff_frm,
                                                   angles[a], out_frms_ptr);
        }

        if (out_calc_descs) {
          /* Set window size that the whole patch is accounted */
          const double windowSize = (double)VL_SIFT_NBP/2.;
          vl_sift_calc_descriptor(patch_grad, patch_sz, patch_sz,
                                  out_descrs_ptr,
                                  pt_pos, pt_pos, 1., angles[a],
                                  self->magnification, windowSize,
                                  self->normThreshold, 0);
          out_descrs_ptr += desc_sz;
        }
      }

      charp(in_frms_ptr) += in_frm_sz;
    }
  } else {
    /* Discs .......................................... */
    VlScaleSpaceFrame last_ssfrm;

    if (in_discard_orient) {
      vl_scalespace_frame_init(gss, &last_ssfrm, 0., 0., 0.);
    }

    if (in_calc_orient || out_calc_descs) {
      vl_scalespace_apply(gss, grad, vl_imgradient_polar_f_callback, 0);
    }

    for (i = 0; i < srcFramesNum; ++i)  {
      VlScaleSpaceFrame ssfrm;
      double angles[VL_COVDET_MAX_ANGLES];
      int nangles;
      int a;

      nangles = 1;

      angles[0] =
        _vl_covdet_convert_to_scalespace_frame(self, &ssfrm, in_frms_ptr,
                                                srcFramesType);

      /* When converting from oriented to unoriented, discard duplicates */
      if (in_discard_orient) {
        if (ssfrm.x == last_ssfrm.x &&
            ssfrm.y == last_ssfrm.y &&
            ssfrm.sigma == last_ssfrm.sigma){
          charp(in_frms_ptr) += in_frm_sz;
          continue;
        } else {
          last_ssfrm = ssfrm;
          angles[0] = 0.;
        }
      }

      if (in_calc_orient) {
        nangles = vl_covdet_calc_ssframe_orientations(grad, angles, &ssfrm);
      }

      /* Store supplemented frames to storage for export to matlab ... */
      for(a = 0; a < nangles; ++a) {
        nframes++;

        if (out_frms_type == srcFramesType) {
          /* When the keypoint type is the same, simply copy the data */
          memcpy(out_frms_ptr, in_frms_ptr, in_frm_sz);
          charp(out_frms_ptr) += in_frm_sz;
        } else {
          out_frms_ptr = _vl_covdet_store_scalespace_frame(self, &ssfrm,
                                                angles[a], out_frms_ptr);
        }

        if (out_calc_descs) {
          vl_covdet_calc_ssframe_descriptor (self, out_descrs_ptr,
                                               &ssfrm, angles [a]) ;
          out_descrs_ptr += desc_sz;
        }
      }
      charp(in_frms_ptr) += in_frm_sz;
    }
  }
  self->numFrames = nframes;
}


/* ---------------------------------------------------------------- */
/*                                               Frame orientations */
/* ---------------------------------------------------------------- */


/** ------------------------------------------------------------------
 ** @brief Calculate the isotropic frame orientation(s)
 **
 ** @param gradient     Pointer to polar gradient image.
 ** @param imageWidth   Original image width (grad. width is @c 2*imageWidth).
 ** @param imageHeight  Original image height.
 ** @param angles       Orientations (output) with size @a VL_KPTS_MAX_ANGLES.
 ** @param frameX       Coordinate x of the frame in @a gradient.
 ** @param frameY       Coordinate y of the frame in @a gradient.
 ** @param frameSigma   Scale of the frame in @a gradient.
 ** @param frameOctave  Octave of the frame.
 **
 ** The function computes the orientation(s) of the frame @a k.
 ** The function returns the number of orientations found (up to
 ** four). The orientations themselves are written to the vector @a
 ** angles.
 **
 ** @return number of orientations found.
 **/

VL_EXPORT int
vl_calc_frame_orientations(float* gradient, int imageWidth, int imageHeight,
                           double angles[VL_COVDET_MAX_ANGLES],
                           float frameX, float frameY,
                           float frameSigma, int frameOctave)
{
  double const winf   = 1.5 ;
  double       xper   = pow (2.0, frameOctave) ;
  int const    xo     = 2 ;         /* x-stride */
  int const    yo     = 2 * imageWidth ;     /* y-stride */
  double       x      = frameX     / xper ;
  double       y      = frameY     / xper ;
  double       sigma  = frameSigma / xper ;

  int          xi     = (int) (x + 0.5) ;
  int          yi     = (int) (y + 0.5) ;

  double const sigmaw = winf * sigma ;

  int          W      = VL_MAX(floor (3.0 * sigmaw), 1) ;

  int          nangles= 0 ;

  enum {nbins = 36} ;

  double hist [nbins], maxh ;
  float const * pt ;
  int xs, ys, iter, i ;

  if(xi < 0            ||
     xi > imageWidth - 1        ||
     yi < 0            ||
     yi > imageHeight - 1 ) {
    return 0;
  }

  /* clear histogram */
  memset (hist, 0, sizeof(double) * nbins) ;

  /* compute orientation histogram */
  pt = gradient;
  pt  += xi*xo + yi*yo ;

#undef  at
#define at(dx,dy) (*(pt + xo * (dx) + yo * (dy)))

  for(ys  =  VL_MAX (- W,       - yi) ;
      ys <=  VL_MIN (+ W, imageHeight - 1 - yi) ; ++ys) {

    for(xs  = VL_MAX (- W,       - xi) ;
        xs <= VL_MIN (+ W, imageWidth - 1 - xi) ; ++xs) {


      double dx = (double)(xi + xs) - x;
      double dy = (double)(yi + ys) - y;
      double r2 = dx*dx + dy*dy ;
      double wgt, mod, ang, fbin ;

      /* limit to a circular window */
      if (r2 >= W*W + 0.6) continue ;

      wgt  = fast_expn (r2 / (2*sigmaw*sigmaw)) ;
      mod  = *(pt + xs*xo + ys*yo    ) ;
      ang  = *(pt + xs*xo + ys*yo + 1) ;
      fbin = nbins * ang / (2 * VL_PI) ;

#if defined(VL_SIFT_BILINEAR_ORIENTATIONS)
      {
        int    bin  = vl_floor_d (fbin - 0.5) ;
        double rbin = fbin - bin - 0.5 ;
        hist [(bin + nbins) % nbins] += (1 - rbin) * mod * wgt ;
        hist [(bin + 1    ) % nbins] += (    rbin) * mod * wgt ;
      }
#else
      {
        int    bin  = vl_floor_d (fbin) ;
        bin = vl_floor_d (nbins * ang / (2*VL_PI)) ;
        hist [(bin) % nbins] += mod * wgt ;
      }
#endif

    } /* for xs */
  } /* for ys */

  /* smooth histogram */
  for (iter = 0; iter < 6; iter ++) {
    double prev  = hist [nbins - 1] ;
    double first = hist [0] ;
    int i ;
    for (i = 0; i < nbins - 1; i++) {
      double newh = (prev + hist[i] + hist[(i+1) % nbins]) / 3.0;
      prev = hist[i] ;
      hist[i] = newh ;
    }
    hist[i] = (prev + hist[i] + first) / 3.0 ;
  }

  /* find the histogram maximum */
  maxh = 0 ;
  for (i = 0 ; i < nbins ; ++i)
    maxh = VL_MAX (maxh, hist [i]) ;

  /* find peaks within 80% from max */
  nangles = 0 ;
  for(i = 0 ; i < nbins ; ++i) {
    double h0 = hist [i] ;
    double hm = hist [(i - 1 + nbins) % nbins] ;
    double hp = hist [(i + 1 + nbins) % nbins] ;

    /* is this a peak? */
    if (h0 > 0.8*maxh && h0 > hm && h0 > hp) {

      /* quadratic interpolation */
      double di = - 0.5 * (hp - hm) / (hp + hm - 2 * h0) ;
      double th = 2 * VL_PI * (i + di + 0.5) / nbins ;
      angles [ nangles++ ] = th ;
      if( nangles == VL_COVDET_MAX_ANGLES )
        break;
    }
  }

  return nangles ;
}


/** ------------------------------------------------------------------
 ** @brief Calculate the frame orientation(s)
 **
 ** @param grad     Gradient scale space of double width.
 ** @param angles   orientations (output)  with min size @a VL_KPTS_MAX_ANGLES.
 ** @param k        frame.
 **
 ** The function computes the orientation(s) of the frame @a k.
 ** The function returns the number of orientations found (up to
 ** four). The orientations themselves are written to the vector @a
 ** angles.
 **
 ** @remark The function requires the frame scale level @c k->s to
 ** be in the range @c s_min+1 and @c s_max-1 (where usually @c
 ** s_min=-1 and @c s_max=S+1) and frame octave @c k->o to be in the
 ** range @c o_min and @c o_max. If this is not the case, the function
 ** returns zero orientations.
 **
 ** This method is here just from compatibility reasons, can be
 ** removed when framedet driver would replace sift driver.
 **
 ** @return number of orientations found.
 **/

int
vl_covdet_calc_ssframe_orientations (VlScaleSpace *grad,
                                       double angles [VL_COVDET_MAX_ANGLES],
                                       VlScaleSpaceFrame const *k)
{
  int          o      = k->o;

  int          w      = vl_scalespace_get_octave_width(grad, o) >> 1 ;
  int          h      = vl_scalespace_get_octave_height(grad, o) ;
  int          smin   = vl_scalespace_get_level_min (grad) ;
  int          smax   = vl_scalespace_get_level_max (grad) ;

  int          si     = k-> is ;

  /* skip the frame if it is out of bounds */
  if(si < smin    ||
     si > smax  ) {
    return 0 ;
  }

  return vl_calc_frame_orientations(vl_scalespace_get_octave(grad, o, si),
                                       w, h, angles,
                                       k->x, k->y, k->sigma, k->o);

}

/* ---------------------------------------------------------------- */
/*                                                  SIFT descriptor */
/* ---------------------------------------------------------------- */


/** ------------------------------------------------------------------
 ** @internal
 ** @brief Normalizes in norm L_2 a descriptor
 ** @param begin begin of histogram.
 ** @param end   end of histogram.
 **/

VL_INLINE float
normalize_histogram
(float *begin, float *end)
{
  float* iter ;
  float  norm = 0.0 ;

  for (iter = begin ; iter != end ; ++ iter)
    norm += (*iter) * (*iter) ;

  norm = vl_fast_sqrt_f (norm) + VL_EPSILON_F ;

  for (iter = begin; iter != end ; ++ iter)
    *iter /= norm ;

  return norm;
}

/** ------------------------------------------------------------------
 ** @brief Compute descriptor of a scale space frame
 **
 ** @param self         Covariant frames detector object
 ** @param descriptor   SIFT descriptor (output)
 ** @param frame        ScaleSpace frame.
 ** @param angle0       Frame orientation.
 **
 ** The function computes the SIFT descriptor of the frame @a frame of
 ** orientation @a angle0. The function fills the buffer @a descriptor
 ** which must be large enough to hold the descriptor.
 **
 ** @note This method is here just from compatibility reasons, can be
 ** removed when framedet driver would replace sift driver.
 **/

VL_EXPORT void
vl_covdet_calc_ssframe_descriptor (VlCovDet *self,
                                   float *descriptor,
                                   VlScaleSpaceFrame const* frame,
                                   double angle0)
{
  VlScaleSpace const *grad = self->grad;
  int          o           = frame->o;
  double       xper        = pow (2.0, o) ;

  int          w           = vl_scalespace_get_octave_width(grad, o) >> 1 ;
  int          h           = vl_scalespace_get_octave_height(grad, o) ;
  double       x           = frame-> x     / xper ;
  double       y           = frame-> y     / xper ;
  double       sigma       = frame-> sigma / xper ;
  int          si          = frame-> is ;
  float const* grad_octave;

  /* check bounds */
  if(si    <  vl_scalespace_get_level_min(grad)    ||
     si    >  vl_scalespace_get_level_max(grad)     )
    return ;

  grad_octave = vl_scalespace_get_octave(grad, o, si);

  vl_sift_calc_descriptor(grad_octave, w, h, descriptor, x, y, sigma, angle0,
                          self->magnification, self->windowSize,
                          self->normThreshold, 1);

}

/** ------------------------------------------------------------------
 ** @brief Calculate SIFT descriptor from raw data
 **
 ** @param gradient         image gradients.
 ** @param descriptor       SIFT descriptor (output).
 ** @param imageWidth       image width.
 ** @param imageHeight      image height.
 ** @param x                frame x coordinate.
 ** @param y                frame y coordinate.
 ** @param sigma            frame scale.
 ** @param angle0           frame orientation.
 ** @param magnification    magnification factor
 ** @param windowSize       size of Gaussian window (in spatial bins)
 ** @param normThreshold    norm threshold
 ** @param borderSize       Size of the border of the grad. planes to ignore
 **
 ** The function runs the SIFT descriptor on raw data. Here @a gradient
 ** is a 2 x @a imageWidth x @a imageHeight array (by convention, the memory
 ** layout is as such the first index is the fastest varying
 ** one). The first @a imageWidth x @a imageHeight layer of the array contains
 ** the gradient magnitude and the second the gradient angle (in
 ** radians, between 0 and @f$ 2\pi @f$). @a x, @a y and @a sigma give
 ** the frame center and scale respectively.
 **
 ** In order to be equivalent to a standard SIFT descriptor the image
 ** gradient must be computed at a smoothing level equal to the scale
 ** of the frame. In practice, the actual SIFT algorithm makes the
 ** following additional approximation, which influence the result:
 **
 ** - Scale is discretized in @c S levels.
 ** - The image is downsampled once for each octave (if you do this,
 **   the parameters @a x, @a y and @a sigma must be
 **   scaled too).
 **/

void
vl_sift_calc_descriptor (float const* gradient, int imageWidth, int imageHeight,
                         float *descriptor,
                         double x, double y, double sigma, double angle0,
                         double magnification, double windowSize,
                         double normThreshold,
                         int borderSize)
{
  /*
     The SIFT descriptor is a three dimensional histogram of the
     position and orientation of the gradient.  There are VL_SIFT_NBP bins for
     each spatial dimension and VL_SIFT_NBO bins for the orientation dimension,
     for a total of VL_SIFT_NBP x VL_SIFT_NBP x VL_SIFT_NBO bins.

     The support of each spatial bin has an extension of SBP = 3sigma
     pixels, where sigma is the scale of the frame.  Thus all the
     bins together have a support SBP x VL_SIFT_NBP pixels wide. Since
     weighting and interpolation of pixel is used, the support extends
     by another half bin. Therefore, the support is a square window of
     SBP x (VL_SIFT_NBP + 1) pixels. Finally, since the patch can be
     arbitrarily rotated, we need to consider a window 2W += sqrt(2) x
     SBP x (VL_SIFT_NBP + 1) pixels wide.
  */

  int          w      = imageWidth ;
  int          h      = imageHeight ;
  int const    xo     = 2 ;         /* x-stride */
  int const    yo     = 2 * w ;     /* y-stride */

  int          xi     = (int) (x + 0.5) ;
  int          yi     = (int) (y + 0.5) ;

  double const st0    = sin (angle0) ;
  double const ct0    = cos (angle0) ;
  double const SBP    = magnification * sigma + VL_EPSILON_D ;
  int    const W      = floor
    (sqrt(2.0) * SBP * (VL_SIFT_NBP + 1) / 2.0 + 0.5) ;

  int const binto = 1 ;          /* bin theta-stride */
  int const binyo = VL_SIFT_NBO * VL_SIFT_NBP ;  /* bin y-stride */
  int const binxo = VL_SIFT_NBO ;        /* bin x-stride */

  int bin, dxi, dyi ;
  float const *pt ;
  float       *dpt ;

  /* check bounds */
  if(xi    <  0               ||
     xi    >= w               ||
     yi    <  0               ||
     yi    >= h -    1        )
    return ;

  /* clear descriptor */
  memset (descriptor, 0, sizeof(float) * VL_SIFT_NBO*VL_SIFT_NBP*VL_SIFT_NBP) ;

  /* Center the scale space and the descriptor on the current frame.
   * Note that dpt is pointing to the bin of center (SBP/2,SBP/2,0).
   */
  pt  = gradient + xi*xo + yi*yo ;
  dpt = descriptor + (VL_SIFT_NBP/2) * binyo + (VL_SIFT_NBP/2) * binxo ;

#undef atd
#define atd(dbinx,dbiny,dbint) *(dpt + (dbint)*binto + (dbiny)*binyo + (dbinx)*binxo)

  /*
   * Process pixels in the intersection of the image rectangle
   * (1,1)-(M-1,N-1) and the frame bounding box.
   */
  for(dyi =  VL_MAX(- W,  borderSize - yi   ) ;
      dyi <= VL_MIN(+ W, h - yi - (1 + borderSize)) ; ++ dyi) {

    for(dxi =  VL_MAX(- W,  borderSize - xi   ) ;
        dxi <= VL_MIN(+ W, w - xi - (1 + borderSize)) ; ++ dxi) {

      /* retrieve */
      float mod   = *( pt + dxi*xo + dyi*yo + 0 ) ;
      float angle = *( pt + dxi*xo + dyi*yo + 1 ) ;
      float theta = vl_mod_2pi_f (angle - angle0) ;

      /* fractional displacement */
      float dx = xi + dxi - x;
      float dy = yi + dyi - y;

      /* get the displacement normalized w.r.t. the frame
         orientation and extension */
      float nx = ( ct0 * dx + st0 * dy) / SBP ;
      float ny = (-st0 * dx + ct0 * dy) / SBP ;
      float nt = VL_SIFT_NBO * theta / (2 * VL_PI) ;

      /* Get the Gaussian weight of the sample. The Gaussian window
       * has a standard deviation equal to VL_SIFT_NBP/2. Note that dx and dy
       * are in the normalized frame, so that -VL_SIFT_NBP/2 <= dx <=
       * VL_SIFT_NBP/2. */
      float const wsigma = windowSize ;
      float win = fast_expn
        ((nx*nx + ny*ny)/(2.0 * wsigma * wsigma)) ;

      /* The sample will be distributed in 8 adjacent bins.
         We start from the ``lower-left'' bin. */
      int         binx = vl_floor_f (nx - 0.5) ;
      int         biny = vl_floor_f (ny - 0.5) ;
      int         bint = vl_floor_f (nt) ;
      float rbinx = nx - (binx + 0.5) ;
      float rbiny = ny - (biny + 0.5) ;
      float rbint = nt - bint ;
      int         dbinx ;
      int         dbiny ;
      int         dbint ;

      /* Distribute the current sample into the 8 adjacent bins*/
      for(dbinx = 0 ; dbinx < 2 ; ++dbinx) {
        for(dbiny = 0 ; dbiny < 2 ; ++dbiny) {
          for(dbint = 0 ; dbint < 2 ; ++dbint) {

            if (binx + dbinx >= - (VL_SIFT_NBP/2) &&
                binx + dbinx <    (VL_SIFT_NBP/2) &&
                biny + dbiny >= - (VL_SIFT_NBP/2) &&
                biny + dbiny <    (VL_SIFT_NBP/2) ) {
              float weight = win
                * mod
                * vl_abs_f (1 - dbinx - rbinx)
                * vl_abs_f (1 - dbiny - rbiny)
                * vl_abs_f (1 - dbint - rbint) ;

              atd(binx+dbinx, biny+dbiny, (bint+dbint) % VL_SIFT_NBO) += weight ;
            }
          }
        }
      }
    }
  }

  /* Standard SIFT descriptors are normalized, truncated and normalized again */
  if(1) {

    /* normalize L2 norm */
    float norm = normalize_histogram (descriptor, descriptor + VL_SIFT_NBO*VL_SIFT_NBP*VL_SIFT_NBP) ;

    /*
       Set the descriptor to zero if it is lower than our
       norm_threshold.  We divide by the number of samples in the
       descriptor region because the Gaussian window used in the
       calculation of the descriptor is not normalized.
     */
    /*
     * NOTE in the former code numSamples
     * was in vl_sift_calc_keypoint_descriptor ignored...
     */
    int numSamples = borderSize ?
      (VL_MIN(W, w - xi -1) - VL_MAX(-W, - xi) + borderSize) *
      (VL_MIN(W, h - yi -1) - VL_MAX(-W, - yi) + borderSize) : 1 ;

    if(normThreshold && norm < normThreshold * numSamples) {
        for(bin = 0; bin < VL_SIFT_NBO*VL_SIFT_NBP*VL_SIFT_NBP ; ++ bin)
            descriptor [bin] = 0;
    }
    else {
      /* truncate at 0.2. */
      for(bin = 0; bin < VL_SIFT_NBO*VL_SIFT_NBP*VL_SIFT_NBP ; ++ bin) {
        if (descriptor [bin] > 0.2) descriptor [bin] = 0.2;
      }

      /* normalize again. */
      normalize_histogram (descriptor, descriptor + VL_SIFT_NBO*VL_SIFT_NBP*VL_SIFT_NBP) ;
    }
  }
}


/* ---------------------------------------------------------------- */
/*                                       Affine shape normalisation */
/* ---------------------------------------------------------------- */

/** ------------------------------------------------------------------
 ** @brief Create a new Affine shape normalisation object
 ** @param imageWidth width of the original image.
 ** @param imageHeight height of the original image.
 ** @param patchSize size of the square patch.
 ** @return the new Affine shape normalisation object.
 **
 ** This class handles normalisation of frame affine anisotropic
 ** neighbourhood into an isotropic one based on the detected
 ** affine shape matrix. The neighbourhood is normalised into
 ** a square patch of size @a patchSize based on the data from
 ** original image of size @a imageWidth @a imageHeight
 **
  ** @sa ::vl_affinepatchnormalizer_delete().
 **/

VL_EXPORT VlAffinePatchNormalizer*
vl_affinepatchnormalizer_new (vl_size imageWidth, vl_size imageHeight, vl_size patchSize)
{
  VlAffinePatchNormalizer* self = vl_malloc(sizeof(VlAffinePatchNormalizer));
  vl_size nel, wss ;

  self->width = imageWidth;
  self->height = imageHeight;
  nel = imageWidth * imageHeight ;

  wss = sizeof(float) * nel;
  self->workspace = vl_malloc(wss);
  self->wss = wss;

  self->patch  = vl_malloc(sizeof(float) * patchSize * patchSize);
  self->patchSize = patchSize;

  self->smoothFilter = vl_imsmooth_new(imageWidth, imageHeight);

  return self;
}

/** -------------------------------------------------------------------
 ** @brief Delete Affine shape normalisation object
 ** @param self Affine shape normalisation object to be deleted
 **
 ** The function frees the resources allocated by ::vl_affinepatchnormalizer_new().
 **/

VL_EXPORT void
vl_affinepatchnormalizer_delete (VlAffinePatchNormalizer *self)
{
  if(self) {
    if(self->workspace) vl_free(self->workspace);
    if(self->patch) vl_free(self->patch);
    if(self->smoothFilter) vl_imsmooth_delete(self->smoothFilter);
    vl_free(self);
  }
}

/** ------------------------------------------------------------------
 ** @brief Normalize the affine shape into a square patch
 ** @param self Affine shape normalisation object.
 ** @param image Image data. Image is of size defined in constructor.
 ** @param frame Frame which neighbourhood is normalized.
 ** @param magnif The magnification of keypoint scale.
 ** @return error code. The function returns 0 if it was not able to normalize
 ** due to fact, that the neighbourhood went over the image boundary.
 **
 ** The function normalises anisotropic neighbourhood of a keypoint into
 ** an isotropic patch based on its affine frame (@sa affine). This step is
 ** needed for calculation of affine invariant SIFT descriptor and
 ** orientations of the keypoint neughbourhood.
 **
 ** The shape of the neighbourhood is determined by the affine transformation
 ** of the frame (parameters frame->axx) and its size is defined by the frame
 ** scale (frm->sigma) and the magnification factor
 ** (@sa vl_affinepatchnormalizer_get_magnif() and @sa vl_affinepatchnormalizer_set_magnif()).
 **
 ** The result is stored into patch ::vl_affinepatchnormalizer_get_patch and
 ** ::vl_affinepatchnormalizer_get_patch_size.
 **
 ** @sa covdet
 **/

VL_EXPORT int
vl_affinepatchnormalizer_normalise (VlAffinePatchNormalizer *self,
                                    float const *image,
                                    VlAffineShapeEstimatorFrame const *frame,
                                    double magnif)
{
  /* TODO implement the more effective way of normalisation where the
   * preblurred levels from the scale space are used.
   */
  float* workspace  = self->workspace;
  size_t wss        = self->wss;
  int    width      = self->width;
  int    height     = self->height;
  float* patch      = self->patch;

  /* half patch size in pixels of image */
  int patchRadiusImage = ceil(frame->sigma * magnif) ;
  int patchSizeImage = 2 * patchRadiusImage + 1 ;
  double patchToImageScale = (double)patchSizeImage / (double)self->patchSize ;

  /* Test if the patch in the image touches boundary of the image,
   * if so, do not normalize */
  {
    int i;
    const int half_patch_size  = self->patchSize >> 1;
    double x_b[4];
    double y_b[4];
    x_b[0] = -half_patch_size; x_b[1] = -half_patch_size;
    x_b[2] = +half_patch_size; x_b[3] = +half_patch_size;
    y_b[0] = -half_patch_size; y_b[1] = +half_patch_size;
    y_b[2] = -half_patch_size; y_b[3] = +half_patch_size;
    for (i=0; i<4; i++)
    {
      double imx = frame->x + frame->a11 * x_b[i] + frame->a12 * y_b[i] ;
      double imy = frame->y + frame->a21 * x_b[i] + frame->a22 * y_b[i] ;
      if (floor(imx) <= 0 ||
          floor(imy) <= 0 ||
          ceil(imx) >= width - 2 ||
          ceil(imy) >= height - 2)
      {
        return 0;
      }
    }
  }

  if (patchToImageScale > 0.4)
  {
    size_t wss_needed;

    /* The pixels in the image are 0.4 apart => the affine deformation
     * leaves +1 border for the bilinear interpolation */
    patchSizeImage += 2;
    wss_needed = patchSizeImage * patchSizeImage * sizeof(float);

    /* Test if the workspace for smoothed image is big enough */
    if (wss_needed >= wss) {
      if (workspace)
        vl_free(workspace);
      workspace = vl_malloc(wss_needed);
      self->workspace = workspace;
      self->wss = wss_needed;
    }

    if (vl_affineshapeestimator_interpolate_bilinear (image, width, height,
                                                      frame->x, frame->y,
                                                      frame->a11, frame->a12,
                                                      frame->a21, frame->a22,
                                                      workspace, patchSizeImage, patchSizeImage))
    {
      vl_imsmooth_smooth_image (self->smoothFilter, workspace, workspace,
                                patchSizeImage, patchSizeImage,
                                1.5f * patchToImageScale) ;

      if (vl_affineshapeestimator_interpolate_bilinear(workspace, patchSizeImage, patchSizeImage,
                                                       (float)(patchSizeImage>>1),
                                                       (float)(patchSizeImage>>1),
                                                       patchToImageScale, 0.,
                                                       0., patchToImageScale,
                                                       patch, self->patchSize, self->patchSize))
      {
        return 1;
      }
    } else {
      return 0;
    }
  } else
  {
    /* If imageToPatchScale is small (i.e. lot of oversampling),
     * affine normalize without smoothing */
    if (vl_affineshapeestimator_interpolate_bilinear (image, width, height,
                                                      frame->x, frame->y,
                                                      patchToImageScale * frame->a11,
                                                      patchToImageScale * frame->a12,
                                                      patchToImageScale * frame->a21,
                                                      patchToImageScale * frame->a22,
                                                      patch, self->patchSize, self->patchSize))
    {
      return 1;
    }
  }
  return 0;
}



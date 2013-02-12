% VL_DEMO_COVDET   Demo: covariant detectors

imagePath = fullfile(vl_root,'data','checker.jpg') ;
im = im2single(rgb2gray(imread(imagePath))) ;

%im = im2single(checkerboard(20)) ;
%im = im(:,round(1:1.5:end)) ;

vl_demo_print = @(x) x ;

% SIFT
peakThreshold = 0.05 ;
peakThreshold = 0.005 ;
edgeThreshold = 100 ;

[frames, descrs, imGss, imResp] = ...
    vl_covdet(im, ...
              'AffineAdaptation', false, ...
              'Orientation', true, ...
              'PeakThreshold', peakThreshold, ...
              'EdgeThreshold', edgeThreshold, ...
              'FirstOctave',-1, ...
              'Method', 'Hessian', ...
              'verbose') ;

figure(1) ; clf ;
imagesc(im) ; colormap gray ; axis image off ; hold on ;
vl_plotframe(frames,'linewidth',3,'color', 'k') ;
vl_plotframe(frames,'linewidth',1,'color', 'y') ;
vl_demo_print('covdet_basic_1') ;

figure(2) ; clf ; vl_plotss(imGss) ; colormap(gray(256)) ;
vl_demo_print('covdet_basic_1_gss') ;

figure(3) ; clf ; vl_plotss(imResp,'uniform',false) ; colormap(gray(256)) ;
vl_demo_print('covdet_basic_1_resp') ;

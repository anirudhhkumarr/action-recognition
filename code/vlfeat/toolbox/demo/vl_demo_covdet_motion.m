function mov = vl_demo_covdet_motion()
% VL_DEMO_COVDET_MOTION  Demo: covariance movie

affine = false ;
peakThreshold = 0.05 ;
edgeThreshold = 1000 ;

if affine
  motion = affMotion() ;
else
  motion = simMotion() ;
end

for t = 1:numel(motion)
  im = geni(motion{t}) ;
  frames = vl_covdet(im, ...
                     'AffineAdaptation', affine, ...
                     'Orientation', true, ...
                     'PeakThreshold', peakThreshold, ...
                     'EdgeThreshold', edgeThreshold, ...
                     'FirstOctave',-1, ...
                     'verbose') ;
  figure(1) ; clf ;
  imagesc(im) ; colormap gray ; axis image off ; hold on ;
  vl_plotframe(frames,'linewidth',3,'color', 'k') ;
  vl_plotframe(frames,'linewidth',1,'color', 'y') ;
  drawnow ;
  if nargout > 0, mov(t) = getframe() ; end
end

function H = simMotion()
tr = linspace(0,1,200) ;
for t = 1:numel(tr)
  th = 4*tr(t) ;
  s = tr(t)/2+.5 ;
  H{t} = [s*cos(th) s*sin(th) 0 ;
          s*sin(th) -s*cos(th) 0 ;
          0 0 1];
end

function H = affMotion()
tr = linspace(0,1,200) ;
for t = 1:numel(tr)
  th = 4*tr(t) ;
  s = tr(t)/2+.5 ;
  H{t} = [1*cos(th) 1*sin(th) 0 ;
          s*sin(th) -s*cos(th) 0 ;
          0 0 1];
end

function im = geni(H)
xr = linspace(-10,10,200) ;
yr = linspace(-10,10,200) ;
[x,y] = meshgrid(xr,yr) ;

[xp,yp] = vl_waffine(H(1:2,1:2),H(1:2,3),x,y) ;
xp = xp ./ (H(3,1) * x + H(3,2) * y + 1) ;
yp = yp ./ (H(3,1) * x + H(3,2) * y + 1) ;
xp = mod(xp+1,2)-1;
yp = mod(yp+1,2)-1;

sigma = .4 ;
im = zeros(numel(yr),numel(xr),'single') ;
im = im + exp(-(xp.^2 + yp.^2) / sigma) ;




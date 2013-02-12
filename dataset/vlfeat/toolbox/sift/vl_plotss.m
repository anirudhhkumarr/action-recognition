function h = vl_plotss(ss, varargin)
% VL_PLOTSS  Plot scale space
%   VL_PLOTSS(SS) plots the scale space SS into the current
%   figure. The scale space structure SS has the format:
%
%   SS.o_idxs:: octave range
%   SS.s_idxs:: octave subdivisions range
%   SS.data:: data for all octaves and subdivisions
%     This is a cell array with NUMEL(SS.o_idxs) rows and
%     NUMEL(SS.s_idxs) columns.
%
%   VL_PLOTSS() accepts the following options:
%
%   Uniform:: true
%     If true, rescales the range of all the scale levles to a common
%     range for visualization. If set to false, rescales each level
%     independently.
%
%   See also: VL_COVDET().

% AUTORIGHTS

opts.uniform = true ;
opts = vl_argparse(opts, varargin) ;

numOctaves = numel(ss.o_idxs) ;
totNumLevels = numel(ss.s_idxs) ;

h = [] ;
for o = ss.o_idxs
  for s = ss.s_idxs
    oi = find(o == ss.o_idxs) ;
    si = find(s == ss.s_idxs) ;
    vl_tightsubplot(numOctaves, totNumLevels, (oi-1)*totNumLevels+si) ;
    h(end+1) = imagesc(ss.data{oi,si}) ;
    axis image off ; hold on ;
    text(0,0,sprintf('%d/%d',o,s),...
         'backgroundcolor','w', ...
         'verticalalignment','top') ;
  end
end

if opts.uniform
  minv = min(cellfun(@(x) min(min(x)), ss.data(:))) ;
  maxv = max(cellfun(@(x) max(max(x)), ss.data(:))) ;
  clim = [minv maxv] ;
  set(cell2mat(get(h, 'parent')), 'clim', [minv maxv]) ;
end

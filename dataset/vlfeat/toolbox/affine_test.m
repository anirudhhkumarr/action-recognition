vl_setup;

%pfx = fullfile(vl_root,'data','graffA.png') ;
%pfx = fullfile(vl_root,'data','a.jpg') ;
%pfx = fullfile(vl_root,'data','patch.pgm') ;
pfx = fullfile(vl_root,'data','blobaff.png') ;
orig_I = imread(pfx) ;

figure(1);
clf;
c_map = gray(255);
colormap(c_map);

plot_frames = true;

if size(orig_I,3) == 3
    I = single(rgb2gray(orig_I)) ;
    image(orig_I);
else
    I = single(orig_I);
    imagesc(orig_I);
end

%% SIFT

%[s_f,s_d,gss,resp,grad] = vl_sift(I) ;
%[fs_f,fs_d,gss,resp,grad] = vl_covdet('dog','oriented_disc',I) ;

%s_h1 = vl_plotframe(s_f) ;
%s_h2 = vl_plotframe(s_f) ;
%set(s_h1,'color','k','linewidth',3) ;
%set(s_h2,'color','y','linewidth',2) ;

%%

[fah_f fah_d] = vl_covdet(I,'Method','hessian','AffineAdaptation',true,'Orientation',true,'PeakThresh',(85.0/3.0)^2) ;

%%

if plot_frames

  fah_h1 = vl_plotframe(fah_f) ;
  fah_h2 = vl_plotframe(fah_f) ;
  set(fah_h1,'color','k','linewidth',3) ;
  set(fah_h2,'color','m','linewidth',2) ;
end

hold on
plot(fah_f(1,:), fah_f(2,:),'g.');

%%

tmpName = 'temp_file';
imgFile = [tmpName '.png'];
featFile = [tmpName '.png.hesaff.sift'];

imwrite(orig_I,imgFile);
args = sprintf(' "%s" ',imgFile);
binPath = '/home/kaja/projects/c/vlfeat_benchmarks/+affineDetectors/thirdParty/cmpHessian/hesaff'
cmd = [binPath ' ' args];

[status,msg] = system(cmd);
if status
error('%d: %s: %s', status, cmd, msg) ;
end

frames = vl_ubcread(featFile,'format','oxford');
%delete(imgFile); delete(featFile);

if plot_frames
    %ah_h1 = vl_plotframe(frames) ;
    %ah_h2 = vl_plotframe(frames) ;
    %set(ah_h1,'color','k','linewidth',3) ;
    %set(ah_h2,'color','r','linewidth',2) ;
end
plot(frames(1,:), frames(2,:),'rx');

%%

pt = (85.0/3.0)^2;
frm_types = {'disc','oriented-disc','ellipse','oriented-ellipse'};
frm_types_num = size(frm_types,1);
frames = cell(frm_types_num,1);
descs = cell(frm_types_num,1);

figure(2);
clf;
colormap(c_map);

for i=1:4
  affine_adapt = i > 2;
  calc_orient = mod(i - 1,2) == 1;
  
  [frames{i} descs{i}] = vl_covdet(I,'Method','hessian','AffineAdaptation',affine_adapt,'Orientation',calc_orient,'PeakThresh',pt);
  
  subplot(2,2,i);
  imagesc(I);
  fah_h1 = vl_plotframe(frames{i}) ;
  set(fah_h1,'color','m','linewidth',2) ;
end

%%

conv_frames = cell(frm_types_num,frm_types_num);
conv_descs = cell(frm_types_num,frm_types_num);

frm_diff = cell(frm_types_num,frm_types_num);
desc_diff = cell(frm_types_num,frm_types_num);

for tp_to=1:4
  for tp_from=1:4
    fprintf('Testing conversion from %s to %s...\n',...
      frm_types{tp_from},frm_types{tp_to});
    
    affine_adapt = false;
    orientation = false;
    
    if tp_to == 1
      affine_adapt = false;
      orientation = false;
    elseif tp_to == 2
      affine_adapt = false;
      orientation = true;
    elseif tp_to == 3
      affine_adapt = true;
      orientation = false;
    elseif tp_to == 4
      affine_adapt = true;
      orientation = true;
    end
    
    [conv_frames{tp_from,tp_to} conv_descs{tp_from,tp_to}] = ... 
      vl_covdet(I,'method','hessian','AffineAdaptation',affine_adapt, ...
                'Orientation',orientation,'Frames',frames{tp_from});
    if size(conv_frames{tp_from,tp_to}) == size(frames{tp_to})
      frm_diff{tp_from,tp_to} = ...
        sum(sum(abs(conv_frames{tp_from,tp_to} - frames{tp_to})));
      desc_diff{tp_from,tp_to} = ...
        sum(sum(abs(conv_descs{tp_from,tp_to} - descs{tp_to})));
      fprintf('\tFrame diff: %f\n\tDesc. diff: %f\n',...
        frm_diff{tp_from,tp_to},desc_diff{tp_from,tp_to});
    else
      fprintf('\tDifferent size \n');
    end
  end
end



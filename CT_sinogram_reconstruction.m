%% Disk and Sinogram data provided by class

% README
% For extra credit asked to display given disk phantom and sinogram
%   backproject the sinogram, filter, and try to reproduce disk phantom 
% 
% The data is stored in MATLAB variable "phantom"
% The sinogram is stored in corresponding MATLAB variable named "sg1" 

clear;

circ = [  0   0 110  2; 
        -65   0  20  1; 
          0   0  35  0; 
         65 -25  25  4
         50  50  7   8];

% Image parameters: number of pixels, size, etc.
nx = 128; ny = 128;
dx = 2;		                        % 2 mm / pixel
x = dx * ((1:nx)'-(nx+1)/2);
y = -dx * ((1:ny)'-(ny+1)/2);
xx = x(:,ones(1,ny));
yy = y(:,ones(1,ny))';

% Generate data for disk phantom
phantom = zeros(nx,ny);
for ii=1:size(circ,1)
    cx = circ(ii,1); cy = circ(ii,2); rad = circ(ii,3); amp = circ(ii,4);
    t = find( ((xx-cx)/rad).^2 + ((yy-cy)/rad).^2 <= 1 );
    phantom(t) = amp * ones(size(t));
end
  
% Image the phantom
figure(1)
imagesc(x, y, phantom')               % NOTE the transpose (') here and the x and y values
colormap('gray')
axis('square')
title('Disk Phantom')
xlabel('Position')
ylabel('Position')

% Geometry parameters
nr = 128;	dr = 2;		            % number of radial samples and ray spacing
na = nr*2;          	            % number of angular samples
r = dr * ((1:nr)'-(nr+1)/2);	    % radial sample positions
angle = (0:(na-1))'/na * pi;	    % angular sample positions

% Compute sinogram for the phantom
rr = r(:,ones(1,na));
sg2 = zeros(nr, na);
for ii=1:size(circ,1)
    cx = circ(ii,1); cy = circ(ii,2); rad = circ(ii,3); amp = circ(ii,4);
    tau = cx * cos(angle) + cy * sin(angle);
    tau = tau(:,ones(1,nr))';
    t = find( (rr-tau).^2 <= rad.^2 );
    if ii > 1, amp = amp - circ(1,4); end	% small disks embedded
    sg2(t) = sg2(t)+amp*2*sqrt(rad^2-(rr(t)-tau(t)).^2);
end

% Sinogram of the phantom
figure(2)
imagesc(r, angle/pi*180, sg2')   % NOTE the transpose (') here and
colormap('gray')                 % the fact that angle is displayed in degrees
title('Sinogram: Disk Phantom')
xlabel('Position')
ylabel('Angle')

sinogram = sg2;                         % disk phantom

fprintf('number of rays = %g', nr);
fprintf('\n');
fprintf('number of views = %g', na);
fprintf('\n');

%% Computing the zeroth moment
projection = sum(sinogram,1);

figure(3)
projection_max=max(projection);
plot(projection);
axis([0,na,0,2*projection_max]);

%% Backprojecting to create laminograms

lamin = zeros(nx,ny);
% backproject when theta is 0
for ia = 1:na
    fprintf('angle %g of %g', ia, na, '\n');

    % smear through ia
    backproj = repmat(sinogram(:,ia),[1,[size(sinogram,1),1]]);
    % rotate projection
    temp = imrotate(backproj, angle(ia)*(180/pi), 'bicubic','crop'); 
    lamin = temp+lamin;
end

figure(4)
imagesc(x, y, lamin'); colormap('gray'); axis('image')  
title('Simple Backprojection Image')
xlabel('X Position')
ylabel('Y Position')

%% Filtering

sinogrampad = sinogram;

% applying wrapped filter
sinogramfilt = real( ifft(fft(sinogram).*(triang(nr)*ones(1,size(sinogrampad,2)))));

figure(5)
plot(r, sinogram(:,64)./max(sinogram(:,64)), '-',...
   r, sinogramfilt(:,64)./max(sinogramfilt(:,64)),':');
title('Filtered sinogram')
xlabel('Position');
ylabel('Angle');

% Display Filtered Sinogram 
figure(6)
imagesc(r, angle/pi *180, sinogramfilt'); colormap('gray'); axis('image')  
title('Filtered sinogram')
xlabel('Postion')
ylabel('Angle')

%% Backproject filtered sinogram

bpf_recon = zeros(nx,ny);
for ia = 1:na
    fprintf('angle %g of %g', ia, na);
    fprintf('\n');

    % backproject at theta = 0
    backprojectfilt = repmat(sinogramfilt(:,ia)',nr,1);
    
    % full rotation
    temp = imrotate(backprojectfilt, angle(ia)*(180/pi), 'bicubic','crop'); 
    bpf_recon = bpf_recon+temp;
end

figure(7)
imagesc(x, y, max(bpf_recon,0)); colormap('gray'); axis('image')  
title('Filtered Backprojection Image')
xlabel('X Position')
ylabel('Y Position')

%% Reconstruction comparison

var = find(r==-7);
% Prepare lines for comparison 
line_phantom = phantom(:,var);
line_lamin = lamin(var,:);
line_bpf = bpf_recon(var,:);

figure(8)
plot(y, line_phantom/max(line_phantom),'r--', ...
  y, line_bpf/max(line_bpf),'b-', ...
  y, line_lamin/max(line_lamin), 'k-');
legend('Original', 'Filtered BP', 'Laminogram')
title('Reconstruction Comparison')
xlabel('Horizontal position (mm)');
ylabel('Normalized attenuation');


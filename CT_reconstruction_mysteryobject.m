%% Recreating object based on given mystery sinogram ('sinogram2')

clear;

load sinogram2; sinogram = sg2; 
circ = [  0   0 110  2; 
        -65   0  20  1; 
          0   0  35  0; 
         65 -25  25  4
         50  50  7   8];
     
% Paramaters
dr=1;
[nr, na] = size(sinogram);
nx=nr; ny=nr; dx=dr;  % 
angle = (0:(na-1))'/na * pi;
r = dr * ((1:nr)'-(nr+1)/2);
fprintf('number of rays = %g', nr);
fprintf('number of views = %g', na);

x = dx * ((1:nx)'-(nx+1)/2);
y = -dx * ((1:ny)'-(ny+1)/2);
xx = x(:,ones(1,ny));
yy = y(:,ones(1,ny))';

% phantom sinogram
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

figure(11)
imagesc(r, angle/pi*180, sg2)   % NOTE the transpose (') here and
colormap('gray')                 % the fact that angle is displayed in degrees
title('Sinogram: Disk Phantom')
xlabel('Position (i.e., Rays)')
ylabel('Angle (i.e., Views)')

%% Implement backprojection 

lamin = zeros(nx,ny);
for ia = 1:na
    % backproject through theta
    backproj2 = repmat(sinogram(:,ia),[1,[size(sinogram,1),1]]);
    % rotation
    temp = imrotate(backproj2, angle(ia)*(180/pi), 'bicubic','crop'); 
    lamin = temp+lamin;
end

lamin=imrotate(lamin,90);
figure(12)
imagesc(x, y, lamin); colormap('gray'); axis('image') 
title('Simple Backprojection Image')
xlabel('X Position')
ylabel('Y Position')


%% Filtering

sinogrampad2 = sinogram;
% full filter
sinogramfilt2 = real( ifft(fft(sinogram).*(triang(nr)*ones(1,size(sinogrampad2,2)))));

figure(13)
plot(r, sinogram(:,64)./max(sinogram(:,64)), '-',...
   r, sinogramfilt2(:,64)./max(sinogramfilt2(:,64)),':');
title('Filtered sinogram')
xlabel('Position');
ylabel('Angle');

figure(14)
imagesc(r, angle/pi *180, sinogramfilt2'); colormap('gray'); axis('image')  
title('Filtered sinogram 2')
xlabel('Position')
ylabel('Angle')


%% Backproject filtered sinograms

bpf_recon2 = zeros(nx,ny);
for ia = 1:na
    fprintf('angle %g of %g', ia, na);
    fprintf('\n');

    % bakcproject through theta = 0;
    backprojectfilt2 = repmat(sinogramfilt2(:,ia)',nr,1);
    % rotation
    temp = imrotate(backprojectfilt2, angle(ia)*(180/pi), 'bicubic','crop'); 
    bpf_recon2 = bpf_recon2+temp;
end

figure(15)
imagesc(x, y, max(bpf_recon2,0)); colormap('gray'); axis('image')  
title('Filtered Backprojection Image 2')
xlabel('x position')
ylabel('y position')

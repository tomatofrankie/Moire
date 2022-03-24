function varargout=immax(I,nmax)

% Periodic image extension
K=20;
I=[I(:,(end-K+2):end) I I(:,1:(K-1))];
I=[I((end-K+2):end,:);I;I(1:(K-1),:)];

[M,N]=size(I);

s=sqrt(M^2+N^2)/4000;
I=imagefilter(I,@(wx,wy)exp(-s^2*(wx.^2+wy.^2)/2));
bnd=@(x,x0)mod(x-1,x0)+1;
[Px,Py]=improm(abs(I));
P=min(Px,Py);
[~,j]=sort(I(:).*(P(:)>=max(P(:))/1000),'descend');
k=mod(j-1,M)+1;l=floor((j-1)/M)+1;

I=I(K:(M-K+1),K:(N-K+1));
k=k-K+1;
l=l-K+1;
[M,N]=size(I);
n=find(k>=1&k<=M&l>=1&l<=N);
k=k(n);l=l(n);

varargout{1}=[k(1:nmax),l(1:nmax)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Px,Py]=improm(I)
[~,Px]=islocalmax(I);
[~,Py]=islocalmax(I');
Py=Py';
if nargout==1
    Px=max(Px,Py);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ia=imagefilter(I,freqresp,extension);
% Syntax: Ia=imagefilter(I,freqresp,extension);
% freqresp is a handle to a function describing the Frequency response of
% the filter. For instance, a Gaussian filtering would be called like:
% sigma=1;Ia=imagefilter(I,@(wx,wy)exp(-sigma^2*(wx.^2+wy.^2)/2));
% If unspecified, the extension is equal to 2 (whole-point mirror image
% extension). Use extension=1 for periodic boundary extension.
%
% Example: Gaussian filtering
% sigma=2;Ia=imagefilter(I,@(wx,wy)exp(-sigma^2*(wx.^2+wy.^2)/2));
if nargin~=3
    extension=2;
end
switch extension
    case 1
        % Periodic image extension
    case 2
        % Whole-point mirror image extension
        K=20;
        I=[I(:,K:-1:2) I I(:,(end-1):-1:(end-K+1))];
        I=[I(K:-1:2,:);I;I((end-1):-1:(end-K+1),:)];
end

[s1,s2]=size(I);
II=fft2(I);

k1=(1:s1)-1;
k2=(1:s2)-1;

% Wrapping frequencies from [0,N-1] to [-N/2,N/2]
N=s1;k1=k1.*(k1<N/2)+(k1-N).*(k1>=N/2);
k1=k1'*ones(1,s2);
N=s2;k2=k2.*(k2<N/2)+(k2-N).*(k2>=N/2);
k2=ones(s1,1)*k2;

% Building the frequency response of the filter
filter=freqresp(2*pi*k1/s1,2*pi*k2/s2);
filter(isinf(filter))=max(abs(filter(~isinf(filter))));

% Testing whether the output should be real-valued (if the input is)
Hermit=freqresp(2*pi*k1/s1,2*pi*k2/s2)-conj(freqresp(-2*pi*k1/s1,-2*pi*k2/s2));
realoutput=(max(abs(Hermit(:)))<=1e-10)&isreal(I);

% Filtering in the frequency domain
II=filter.*II;

% Retrieval of the spatial image
Ia=ifft2(II);
if realoutput&isreal(I)
    Ia=real(Ia);
end

switch extension
    case 1
        % Periodic image extension
    case 2
        % Whole-point mirror image extension
        Ia=Ia(K:(s1-K+1),K:(s2-K+1));
end
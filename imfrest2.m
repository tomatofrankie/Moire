function [Imin,param]=imfrest2(varargin)
% [Imin,param]=imfrest2(I,K,[options])
% Fits the complex image I using a model of the form
% \sum_{k=1}^K a(k)*\exp(i*w(1,k)*x+i*w(1,k)*y)
%
% EXAMPLE
% K=20;load lena256,B=I;[M,N]=size(B);B=randn(M,N);
% w0=2*pi*(rand(2,K)-0.5);a0=randn(K,1)+i*randn(K,1);[x,y]=coordinates([M,N]);I0=zeros(M,N);I0(:)=exp(i*x(:)*w0(1,:)+i*y(:)*w0(2,:))*a0;
% P=0;c=10^(P/20)/norm(I0(:))*norm(B(:));a0=a0*c;I0=I0*c;I=I0+B;sigma=norm(I(:)-I0(:))/sqrt(M*N);s=rng;
% rng(s),tic,[Imin,param]=imfrest2(I,K);toc,disp([param.RMSE,sigma]),semilogy(param.err)

I=varargin{1};
K=varargin{2};
[M,N]=size(I);

if nargin<=2
    options=struct;
else
    options=varargin{3};
end

if isfield(options,'numiter')
    numiter=options.numiter;
else
    numiter=10;
end

if isfield(options,'winit')
    w=options.winit;
else
    kl=immax(abs(fft2(I)),K);
    w(1,:)=2*pi*(kl(:,1)'-1)/M;
    w(2,:)=2*pi*(kl(:,2)'-1)/N;
    options.winit=w;
end

x=1:M;x=x'*ones(1,N);
y=1:N;y=ones(M,1)*y;
x0=(M+1)/2;
y0=(N+1)/2;
x=x-x0;
y=y-y0;

n=(1:M*N)';
Imin=I;
prevImin=I;
err=zeros(1,numiter);
for j=1:numiter
    J=exp(i*x(:)*w(1,:)+i*y(:)*w(2,:));
    a=J\I(:);
    Imin(:)=J*a;
    dw=0*w;
    for k=1:K
        A=[ones(M*N,1),i*x(:),i*y(:)].*J(:,k);
        da=A\(I(:)-Imin(:));
        dw(:,k)=da(2:3)/(a(k)+da(1));
    end
    err(j)=norm(prevImin(:)-Imin(:))/sqrt(M*N);
    w=w+real(dw);
    prevImin=Imin;
end

x=x+x0;
y=y+y0;
J=exp(i*x(:)*w(1,:)+i*y(:)*w(2,:));
a=J\I(:);
Imin(:)=J*a;
RMSE=norm(I(:)-Imin(:))/sqrt(M*N);
%e = exp(i*x(:)*param.w(1,:)+i*y(:)*param.w(2,:));
%param=struct('w',w,'a',a,'RMSE',RMSE,'model',@(x,y,param)exp(i*x(:)*param.w(1,:)+i*y(:)*param.w(2,:))*(param.a'*[1 0 0;0 2 0;0 0 2])','winit',options.winit,'err',err);
param=struct('w',w,'a',a,'RMSE',RMSE,'model',@(x,y,param)real(param.a(1))+real(2*param.a(2)*exp(1i*(param.w(1,3))*x+1i*param.w(2,3)*y)),'winit',options.winit,'err',err);



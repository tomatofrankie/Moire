clc;    
close all;  
clear;  
workspace;  
format long g;
format compact;
fontSize = 24;

%%
H = 50;K=50;
[y,x]=meshgrid(1:K,1:H);
moireImage0= imread('moire2.jpg');
[H,K,~] = size(moireImage0)
subplot(311),imshow(moireImage0)
%%
[R,G,B] = imsplit(moireImage0);
%moireImage0=(im2gray(moireImage0));
subplot(334),imshow(R),title('Red Channel')
subplot(335),imshow(G),title('Red Channel')
subplot(336),imshow(B),title('Blue Channel')
%moireImage0=double(moireImage0);

%%
X = double(reshape(moireImage0, H*K, 3));
coeff = pca(X);

Itransformed = X * coeff;

pca1Image = reshape(Itransformed(:,1), H,K);
pca2Image = reshape(Itransformed(:,2), H,K);
pca3Image = reshape(Itransformed(:,3), H,K);
subplot(337),imshow(pca1Image,[]),title('First principal component')
subplot(338),imshow(pca2Image,[]),title('Second principal component')
subplot(339),imshow(pca3Image,[]),title('Third principal component')


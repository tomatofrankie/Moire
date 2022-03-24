H=512;K=512;

x=1:H;
x=repmat(x',[1,K]);
imshow(x);
box
axis on
y=1:K;
y=repmat(y,[H,1]);
imshow(y)
axis on
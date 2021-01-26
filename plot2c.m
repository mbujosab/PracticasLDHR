function plot2c(A)
%PLOT2D	Two dimensional plot.
%	X is matrix with two rows and any number of columns.
%	PLOT2D(X) plots these columns as points in the plane
%	and connects them with lines.
%	The scale is set to [-10,10] in both directions.
%
%	For example, the statements
%           hand
%	or
%	    house
%	create sample X's representing common figures.
%	Then, for various 2 by 2 matrices A,
%	    plot2d(A*X)
%	demonstrates the effect of multiplication by A.
%	   
S=sum(A')';
X=kron(A,[0,1]);
x = X(1,:)';
y = X(2,:)';
plot(x,y,'-',x,y,'*');
%axis([-10 10 -10 10])
%axis('square')
axis("equal")
#box("off")
grid("on")
hold on

pause
SS=kron(S,[0,1]);
xs = SS(1,:)';
ys = SS(2,:)';
plot(xs,ys,'-2r',xs,ys,'*3r');

O=kron(A,[1 0])+kron(S,[0,1,0,1]);
xo = O(1,:)';
yo = O(2,:)';
plot(xo,yo,'--g');
%axis("auto")
m=max([X,S]')-min([X,S]');
a=[min([X,S]')-.1*m(1);max([X,S]')+.1*m(2)]
a=a(:)';
eje=kron(kron(a,[1 0;0 1]),[0 1]);
xe = eje(1,:)';
ye = eje(2,:)';
plot(xe,ye,':k');

axis([a])
hold off
function plot3c(A,b)

S=A*b;
X=kron([A,S],[0,1]);
x = X(1,:)';
y = X(2,:)';
z = X(3,:)';
plot3(x(1:end-2),y(1:end-2),z(1:end-2),'-',x(1:end-2),y(1:end-2),z(1:end-2),'*',x(end-1:end),y(end-1:end),z(end-1:end),'-2r',x(end-1:end),y(end-1:end),z(end-1:end),'*3r');
axis("equal");
view (110, 10);

grid("on")
hold on

m=max([X]')-min([X]');
a=[min([X]')-.1*m(1);max([X]')+.1*m(2)];

w=kron(eye(3),a(:)');
eje=kron(w(:,[1 2 9 10 17 18]),[0 1]);
xe = eje(1,:)';
ye = eje(2,:)';
ze = eje(3,:)';
plot(xe,ye,':k');
box off

hold off

function h = layertrans2(x,y,A,B)

% Transparency plot of two matrices. This particular code avoids annoyances
% that arise when trying to use surf with transparency options.

%% Basic assertions on input sizes
assert(all(size(A)==size(B)),'A and B incompatible');
assert(size(A,2)==numel(x)-1,'x mismatch');
assert(size(A,1)==numel(y)-1,'y mismatch');

%% Affine transformations to unit interval
A = affineunit(A);
B = affineunit(B);

%% Plot
hold on;
for j = 1:size(A,1)
    for k = 1:size(A,2)
        pA = patch([x(k),x(k+1),x(k+1),x(k)]+.5*(x(k+1)-x(k)),...
            [y(j),y(j),y(j+1),y(j+1)]+.5*(y(j+1)-y(j)),'b');
        set(pA,'FaceAlpha',A(j,k),'EdgeAlpha',0);
        pB = patch([x(k),x(k+1),x(k+1),x(k)]+.5*(x(k+1)-x(k)),...
            [y(j),y(j),y(j+1),y(j+1)]+.5*(y(j+1)-y(j)),'r');
        set(pB,'FaceAlpha',B(j,k),'EdgeAlpha',0);
    end
end
xlim(mean(x(1:2))+[min(x),max(x)]);
ylim(mean(y(1:2))+[min(y),max(y)]);

end

%% LOCAL FUNCTION

function A = affineunit(A)

minA = min(A(:));
maxA = max(A(:));
constA = (minA==maxA);  % 1 iff A is constant
A = (A-minA)/(maxA-minA+constA);

end

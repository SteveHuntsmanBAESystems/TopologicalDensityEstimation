% Script for plotting performance data for evaluating CV and TDE on the
% densities $f_j$ ($1 \le j \le 6$).

figure;

%% Normalization macro
normmacro = ['A = A./repmat(sum(A,1),[size(A,1),1]);',...
    'B = B./repmat(sum(B,1),[size(B,1),1]);'];

%% Plot macro
% "box on" gives inadequate results but retained for tickmarks; manually
% set axis and put in a box since "box on" doesn't give adequate results
plotmacro = ['layertrans2(0:5,[L,Inf],A,B);box on;',... 
    'set(gca,''XTick'',1:5,''XTickLabel'',',...
    '{''25'',''50'',''100'',''200'',''500''});',...
    'hold on;axis([0.5,5.5,lo,hi]);ax = axis;',...
    'line([ax(1),ax(1)]+1e-6,[ax(3),ax(4)],''Color'',''k'');',...
    'line([ax(2),ax(2)]-1e-6,[ax(3),ax(4)],''Color'',''k'');',...
    'line([ax(1),ax(2)],[ax(3),ax(3)]+1e-6,''Color'',''k'');',...
    'line([ax(1),ax(2)],[ax(4),ax(4)]-1e-6,''Color'',''k'');'];

%% c1
for j = 1:6
    % Normalized histogram
    lo = min([min(min(pre_c1_top(:,j,:))),min(min(pre_c1_CV(:,j,:)))]);
    hi = max([max(max(pre_c1_top(:,j,:))),max(max(pre_c1_CV(:,j,:)))]);
    L = linspace(lo,hi,50);
    for k = 1:5 % n = [25,50,100,200,500]
        A(:,k) = histc(squeeze(pre_c1_top(:,j,k)),L,1); 
        B(:,k) = histc(squeeze(pre_c1_CV(:,j,k)),L,1); 
    end
    eval(normmacro);
    % Plot
    subplot(3,6,j);  
    eval(plotmacro);
    title(['$f_',num2str(j),'$'],'Interpreter','latex');
end
subplot(3,6,1);
ylabel('$\hat h - h_{opt}$','Interpreter','latex');

%% c2, c3
for j = 1:6
    % Normalized histogram
    lo = 0;
    hi = max([max(max(pre_c23_top(:,j,:))),max(max(pre_c23_CV(:,j,:)))]);
    L = linspace(lo,hi,50);
    for k = 1:5 % n = [25,50,100,200,500]
        A(:,k) = histc(squeeze(pre_c23_top(:,j,k)),L,1); 
        B(:,k) = histc(squeeze(pre_c23_CV(:,j,k)),L,1); 
    end
    eval(normmacro);
    % Plot
    subplot(3,6,j+6);  
    eval(plotmacro);
end
subplot(3,6,7);
ylabel('ISE$(\hat h)$','Interpreter','latex');

%% c4, c5
for j = 1:6
    % Form histograms (note logspace). NB. The only difference between L1
    % and L2 versions is rescaling    
    lo = -4;
    hi = 0.25;
    L = linspace(lo,hi,50);
    % abs turns out to not have any effect, but included for correctness
    for k = 1:5 % n = [25,50,100,200,500]
        A(:,k) = histc(squeeze(log10(abs(pre_c45_top(:,j,k)))),L,1); 
        B(:,k) = histc(squeeze(log10(abs(pre_c45_CV(:,j,k)))),L,1); 
    end
    eval(normmacro);    
    % Plot
    subplot(3,6,j+12);  
    eval(plotmacro);
    xlabel('$n$','Interpreter','latex');
end
subplot(3,6,13);
ylabel('$c_{45}$','Interpreter','latex');

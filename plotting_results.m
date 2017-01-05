% Plots of results 

% Plotting of the load factor vs slenderness, (:,1 to 4,:) gives 4 combinations (or 12 comb for (1 to 3,1 to 4,:))

figure ('name','Load factor vs diameter')

% diameter range
h=[300 500 700 900];

% corresponding minimum load factor
a= Y(:,:,1);

plot(h,a,'o');

% attempt to show the percentage
h1=[300];
h2=[500];
h3=[700];
h4=[900];

a1= Y(1,1,1);
a2= Y(1,2,1);
a3= Y(1,3,1);
a4= Y(1,4,1);

percentage = ['  L',num2str(ceil(classification{1,1,1}(1))), ' D',num2str(ceil(classification{1,1,1}(2))), ' G',num2str(ceil(classification{1,1,1}(3)))];
text(h1,a1,percentage,'HorizontalAlignment','left');

percentage = ['  L',num2str(ceil(classification{1,2,1}(1))), ' D',num2str(ceil(classification{1,2,1}(2))), ' G',num2str(ceil(classification{1,2,1}(3)))];
text(h2,a2,percentage,'VerticalAlignment','bottom');

percentage = ['  L',num2str(ceil(classification{1,3,1}(1))), ' D',num2str(ceil(classification{1,3,1}(2))), ' G',num2str(ceil(classification{1,3,1}(3)))];
text(h3,a3,percentage,'VerticalAlignment','bottom');

percentage = ['  L',num2str(ceil(classification{1,4,1}(1))), ' D',num2str(ceil(classification{1,4,1}(2))), ' G',num2str(ceil(classification{1,4,1}(3)))];
text(h4,a4,percentage,'VerticalAlignment','bottom');

% Plotting of the load factor vs diameter, (:,:,1 to 9) gives 9 combinations (or 27 comb for (1 to 3,:,1 to 9)

figure ('name','Load factor vs slenderness')

% % slenderness range
z= [70 80 90 100 110 120 130 140 150];

% % corresponding minimum load factor
q= Y(3,1,:);

plot(z,squeeze(q),'o');

g1= [70];
g2= [80];
g3= [90];
g4= [100];

b1= Y(1,1,1);
b2= Y(1,1,2);
b3= Y(1,1,3);
b4= Y(1,1,4);


% attempt to show the percentage


percentage = ['  L',num2str(ceil(classification{1,1,1}(1))), ' D',num2str(ceil(classification{1,1,1}(2))), ' G',num2str(ceil(classification{1,1,1}(3)))];
text(g1,b1,percentage,'HorizontalAlignment','right');

percentage = ['  L',num2str(ceil(classification{1,1,2}(1))), ' D',num2str(ceil(classification{1,1,2}(2))), ' G',num2str(ceil(classification{1,1,2}(3)))];
text(g2,b2,percentage,'HorizontalAlignment','left');

percentage = ['  L',num2str(ceil(classification{1,1,3}(1))), ' D',num2str(ceil(classification{1,1,3}(2))), ' G',num2str(ceil(classification{1,1,3}(3)))];
text(g3,b3,percentage,'HorizontalAlignment','right');

percentage = ['  L',num2str(ceil(classification{1,1,4}(1))), ' D',num2str(ceil(classification{1,1,4}(2))), ' G',num2str(ceil(classification{1,1,4}(3)))];
text(g4,b4,percentage,'HorizontalAlignment','left');


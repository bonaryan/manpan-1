% Plots of results 

% Plotting of the load factor vs diameter, (:,:,1 to 10) gives 10 combinations (or 30 comb for (1 to 3,:,1 to 10))

figure ('name','Load factor vs diameter')

% diameter range
h=[300 500 700 900];

% corresponding minimum load factor
a= Y(1,:,1);

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

percentage = ['lokal = ',num2str(ceil(classification{1,1,1}(1))),'dist = ',num2str(ceil(classification{1,1,1}(2))), 'global = ',num2str(ceil(classification{1,1,1}(3)))];
text(h1,a1,percentage,'HorizontalAlignment','left');

percentage = ['lokal = ',num2str(ceil(classification{1,2,1}(1))), 'dist = ',num2str(ceil(classification{1,2,1}(2))), 'global = ',num2str(ceil(classification{1,2,1}(3)))];
text(h2,a2,percentage,'VerticalAlignment','bottom');

percentage = ['lokal = ',num2str(ceil(classification{1,3,1}(1))), 'dist = ',num2str(ceil(classification{1,3,1}(2))), 'global = ',num2str(ceil(classification{1,3,1}(3)))];
text(h3,a3,percentage,'VerticalAlignment','bottom');

percentage = ['lokal = ',num2str(ceil(classification{1,4,1}(1))), 'dist = ',num2str(ceil(classification{1,4,1}(2))), 'global = ',num2str(ceil(classification{1,4,1}(3)))];
text(h4,a4,percentage,'VerticalAlignment','bottom');


% Plotting of the load factor vs slenderness, (:,1 to 4,:) gives 4 combinations (or 12 comb for (1 to 3,1 to 4,:))

figure ('name','Load factor vs slenderness')

% % slenderness range
z= [10 70 80 90 100 110 120 130 140 150];

% % corresponding minimum load factor
b= Y(1,1,:);

plot(z,squeeze(b),'o');


% Hur kan vi plotta kombinationen ovan samt ange tillhorande interaktion (classification varldet)


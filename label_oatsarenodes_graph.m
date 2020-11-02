

% Script to read a movie, and label the G = G(E,N) formed by Physarum, where
% edges E are plasmodium, nodes N are oat flakess

% files:

% path to movie (a folder of pics)
MOVIE_PATH = '/Volumes/Asa Physarum Backup/Spanning Trees/DoubleRingExperimentsData/Physarum on middle/2020_07_13_2cm_DoubleRings_1/Movie/';
% path to mask; draw black ellipses around oats (define nodes)
MASK_PATH = '/Volumes/Asa Physarum Backup/Spanning Trees/DoubleRingExperimentsData/Physarum on middle/2020_07_13_2cm_DoubleRings_1/2020_09_06_19OatsWithPhOnThem_WithLid_MASK.tif';
% path to background mask; mask out area outside of experiment
MASK_PATH_BG = '/Volumes/Asa Physarum Backup/Spanning Trees/DoubleRingExperimentsData/Physarum on middle/2020_07_13_2cm_DoubleRings_1/2020_09_06_19OatsWithPhOnThem_WithLid_MASK_BG.tif';

% time points to read
FRAME_RANGE = [80 200];
STEP = 2;  % step size thru FRAME_RANGE
% 
THRESHOLD = 1;

STREL_SIZE = 31;
SIGMA_1 = 0.5;
MIN_LENGTH = 20;
MAX_DIST_TO_OAT = 35;

WINDOW = 1;
percent = 99;
% THRESH_COEFF_2 = 2;
% SIGMA_2 = 1;

AREA_THRESH = 10;

TIME = 0;

%% first, label nodes


% lets label the oats
MASK = logical(rgb2gray(imread(MASK_PATH)));
oat_cc = bwconncomp(~MASK); % assumes oats were labeled with 0s
all_centroids = regionprops(oat_cc,'Centroid');

% some noise might be in background still, so throw away small stuff (oats
% are wayyy bigger than noise)
oat_centroids = zeros([19,2]);  % row i holds (x,y) coords for centroid of oat i
oat_pixels = cell([19,1]);      % holds linear indices of pixels constituting each oat

label = 0;
for i=1:oat_cc.NumObjects
    if length(oat_cc.PixelIdxList{i}) < 3000  % noise
        continue
    else   % real oat
        label = label + 1;
        oat_centroids(label,:) = all_centroids(i).Centroid;
        oat_pixels{label} = oat_cc.PixelIdxList{i};
    end
end

% visualize nodes
figure(); imshow(~MASK);
hold on; scatter(oat_centroids(:,1),oat_centroids(:,2),100,'filled','MarkerFaceColor','red')
for i=1:19
    text(oat_centroids(i,1),oat_centroids(i,2),num2str(i),'FontSize',25);
end


%% Read movie
tic;
images = dir([MOVIE_PATH, '*.png']);
FRAMES = (FRAME_RANGE(2)-FRAME_RANGE(1) + 1)/STEP;  
[X,Y,~] = size(imread([MOVIE_PATH, images(1).name]));

movie = zeros([X,Y,FRAMES]);
for i=FRAME_RANGE(1):STEP:FRAME_RANGE(2)
    movie(:,:,i/STEP-FRAME_RANGE(1)/STEP+1) = double(rgb2gray(imread([MOVIE_PATH, images(i).name]))); %.* (MASK_BG&MASK);
end
disp(['------------ read: ', num2str(toc),'s ------------'])
TIME = TIME + toc;


%% Process movie (binarize & skeletonize)
tic;

movie_bw = zeros(size(movie)); % hold binary pictures
movie_skel = zeros(size(movie));  % hold skeletonized pictures

for i=1:FRAMES
    
    [this_bw,this_skel] = process_im(movie(:,:,i), STREL_SIZE, THRESHOLD);
    movie_bw(:,:,i) = this_bw;  
    movie_skel(:,:,i) = this_skel;  

    if mod(i,floor(FRAMES/10)) == 0
        disp([num2str(i) ' of ' num2str(FRAMES) ' processed']);
    end
end

% save
movie_stats.movie_bw = movie_bw;
movie_stats.movie_skel = movie_skel;

disp(['------------ processed: ', num2str(toc),'s ------------'])
TIME = TIME + toc;

% look at sample frame 
figure(); imagesc(movie_bw(:,:,10) + ~MASK*2);
title('green = nodes. light blue = binary plasmodium')
figure(); imagesc(imdilate(movie_skel(:,:,10),strel('disk',2)) + ~MASK*2);
title('green = nodes. light blue = skeleton')

%% Compute graph
tic;
movie_graphs = cell([FRAMES, 1]);
for i=1:FRAMES
    
    if sum(movie_skel(:,:,i),'all') < MIN_LENGTH
        % skip empty frames
        continue
    end
    % get nodes/links from Skel2Graph3D
    [~, nodes, links] = Skel2Graph3D(movie_skel(:,:,i), MIN_LENGTH);
    
    movie_graphs{i}.nodes = nodes;
    movie_graphs{i}.links = links;
    if mod(i,floor(FRAMES/10)) == 0
        disp([num2str(i) ' of ' num2str(FRAMES) ' processed']);
    end
end

% save 
movie_stats.movie_graphs = movie_graphs;
movie_stats.FRAME_RANGE=FRAME_RANGE;
movie_stats.THRESHOLD = THRESHOLD;
movie_stats.SIGMA_1 = SIGMA_1;
movie_stats.MIN_LENGTH = MIN_LENGTH;
movie_stats.WINDOW = WINDOW;
movie_stats.percent = percent;
movie_stats.AREA_THRESH = AREA_THRESH;
movie_stats.TIME = TIME;

disp(['------------ graphs computed: ', num2str(toc),'s ------------'])
TIME = TIME + toc;

%% Label graph
% now we construct an adjacency matrix for our desired graph, using the
% graph computed by Skel2Graph3D (which gets nodes/links from plasmodium)
tic;
newlinks = cell([FRAMES, 1]);
oat_connections = cell([FRAMES, 19]);
for i=1:FRAMES   % look at every frame
    
    % first, lets collapse each group of links into one of "label" links
    links = movie_graphs{i}.links;
    NUM_LABELS = double(max([links.label]));
    
    newlinks{i} = cell([NUM_LABELS, 1]);
    
    %lets assign each new link its pixels
    for l=1:length(links)
        newlinks{i}{links(l).label} = [newlinks{i}{links(l).label} links(l).point];
    end
    
    % now for this frame, we have the new links. 
    % test connectivity to each oat flake:
    
    for o=1:19
        [o1,o2,o3] = ind2sub([1200 1600], oat_pixels{o});
        ocoords = [o1,o2,o3];
        
        for l=1:NUM_LABELS
            if isempty(newlinks{i}{l})
                continue
            end
            [i1,i2,i3] = ind2sub([1200 1600], newlinks{i}{l});
            distances = pdist2([i1;i2;i3]',ocoords);
            
            if min(distances(:)) < MAX_DIST_TO_OAT  % oat 'o' is close enough to link 'l' for them to be connected
                oat_connections{i,o} = [oat_connections{i,o} l];
            end
        end
    end
    
    if mod(i,10) == 0
        disp(['labeled frame ' num2str(i)])
    end
end


% now we build adj matrix
adj_mat = zeros([FRAMES, 19, 19]);

for i=1:FRAMES
    for o=1:19  
        for o2=1:19
            if o==o2 
                continue
            else
                conn = oat_connections{i,o}(ismember(oat_connections{i,o},oat_connections{i,o2}));
                if isempty(conn)
                    continue
                else
                    for j=1:length(conn)
                        adj_mat(i,o,o2) = adj_mat(i,o,o2) + length(newlinks{i}{conn(j)});
                    end
                end
            end
        end
    end
end

s = cell([FRAMES, 1]);
t = cell([FRAMES, 1]);
w = cell([FRAMES, 1]);


for i=1:FRAMES
    track = zeros([19 19]);
    for o=1:19
        for o2=1:19
            if o==o2
                continue
            else
                if adj_mat(i,o,o2) == 0
                    continue
                else
                    if track(o,o2) == 1 || track(o2,o) == 1
                        continue
                    else
                        s{i} = [s{i} o];
                        t{i} = [t{i} o2];
                        w{i} = [w{i} adj_mat(i,o,o2)];
                        track(o,o2) = 1;
                        track(o2,o) = 1;
                    end
                end
            end
        end
    end
end

TIME = TIME + toc;

movie_stats.s = s;
movie_stats.t = t;
movie_stats.w = w;
movie_stats.adj_mat = adj_mat;
movie_stats.oat_connections = oat_connections;
movie_stats.newlinks = newlinks;
movie_stats.oat_centroids = oat_centroids;
movie_stats.oat_pixels = oat_pixels;
movie_stats.TIME = TIME;
movie_stats.MOVIE_PATH = MOVIE_PATH;

%% look at an example graph 
figure();
ex = 10; % which timepoint to grab
g = graph(s{ex},t{ex},w{ex});
h = plot(g,'XData',oat_centroids(1:max([s{ex} t{ex}]),1),'YData',oat_centroids(1:max([s{ex} t{ex}]),2),'LineWidth',4);
view(0,90); xlim([0 1600]); ylim([0 1200]);

%% compute network stats as function of time

totlength = zeros([FRAMES, 1]);         % total length of network at each time
totnumedges = zeros([FRAMES, 1]);       % total # edges at each time
edgeavg = zeros([FRAMES, 1]);           % avg edge weight at each time
totnodes = zeros([FRAMES, 1]);          % tot # connected oat flakes at each tme
nodedegree_avg = zeros([FRAMES,1]);         % avg node degree at each time
nodedegree_std = zeros([FRAMES,1]);         % std node degree at each time

graphs = cell([FRAMES, 1]);
for i=1:FRAMES
    totlength(i) = sum(w{i});
    totnumedges(i) = length(s{i});
    edgeavg(i) = mean(w{i});
    totnodes(i) = length(unique([s{i} t{i}]));
    
    graphs{i} = graph(s{i},t{i},w{i});
    these_degrees = degree(graphs{i});
    nodedegree_avg(i) = mean(these_degrees(these_degrees>0));
    nodedegree_std(i) = std(these_degrees(these_degrees>0));
end

movie_stats.totlength = totlength;
movie_stats.totnumedges = totnumedges;
movie_stats.edgeavg = edgeavg;
movie_stats.totnodes = totnodes;
movie_stats.nodedegree_avg = nodedegree_avg;
movie_stats.nodedegree_std = nodedegree_std;

%% compute convex hull
convex_areas = zeros([FRAMES, 1]);
convex_coords = cell([FRAMES, 1]);
for i=1:FRAMES
    g = graph(s{i},t{i},w{i});
    on_nodes = unique([s{i} t{i}]);
    activecoords = [];
    for j=1:length(on_nodes)
        thisnode = on_nodes(j);
        activecoords = [activecoords; oat_centroids(thisnode,:)];
    end
    if size(activecoords,1) < 3  % cant compute convex hull for line
        continue
    else
        [k,convex_areas(i)] = convhull(activecoords);
        convex_coords{i} = [activecoords(k,1), activecoords(k,2)];
    end
end
% 
movie_stats.convex_areas = convex_areas;
movie_stats.convex_coords = convex_coords;

%% I used this code to make gifs 
% gifname = '2020_09_06_19OatsWithPhOnThem_WithLid_graph.gif';
% figure();
% 1 = 1;
% resize = [600 800];
% DelayTime = 0.1;
% imshow(imadjust(uint8(movie(:,:,start))));
% hold on; plot(convex_coords{start}(:,1),convex_coords{start}(:,2),...
%     'LineWidth',2,'Color','red');
% set(gcf,'position',[1 1 600 800]);
% f = getframe(gcf);
% [A,map] = rgb2ind(f.cdata,256);
% imwrite(A,map,['/Users/asakalish/Desktop/',gifname],'gif','LoopCount',Inf,'DelayTime',DelayTime);
% 
% for i=start+1:FRAMES
%     imshow(imadjust(uint8(movie(:,:,i))));
%     if isempty(convex_coords{i}) || convex_areas(i) < 1
%         continue
%     else
%         hold on; 
%         g = graph(s{i},t{i},w{i});
%         h = plot(g,'XData',oat_centroids(1:max([s{i} t{i}]),1),'YData',oat_centroids(1:max([s{i} t{i}]),2),...
%             'LineWidth',6,'EdgeColor','red','NodeColor','red');
%         plot(convex_coords{i}(:,1),convex_coords{i}(:,2),'LineWidth',1,'Color','yellow','LineStyle','--');
%         title(num2str(i))
%     end
%     drawnow;
%     set(gcf,'position',[1 1 600 800]);
%     f = getframe(gcf);
%     [A,map] = rgb2ind(f.cdata,256);
%     if isempty(map)
%         continue
%     else
%         imwrite(A,map,['/Users/asakalish/Desktop/',gifname],'gif','WriteMode','append','DelayTime',DelayTime);
%     end
% end


%% I used this code to plot stuff
% timex = (FRAME_RANGE(1):FRAME_RANGE(2))*10/60;  % in hours; skip first 29 frames
% % timex = FRAME_RANGE(1):FRAME_RANGE(2);
% figure('Name', '2020_09_19_OatInMiddle_1');
% 
% totnodes(totnodes==0) = NaN;
% idx = ~isnan(totnodes);
% 
% subplot(3,2,1); 
% plot(timex(idx), totlength(idx)/totlength(36),'Color','red','LineWidth',1); hold on;
% title('Total network length'); ylabel('Network length / initial network length'); xlabel('Time (hours)');
% ylim([0 20]); 
% yline(1);
% 
% 
% subplot(3,2,3); 
% plot(timex(idx), totnumedges(idx), 'Color','red', 'LineWidth',1); hold on;
% title('Number of edges'); ylabel('# edges / initial # edges'); xlabel('Time (hours)');
% ylim([0 20]);
% 
% subplot(3,2,5); 
% plot(timex(idx), edgeavg(idx)/edgeavg(36), 'Color','red', 'LineWidth',1); hold on;
% title('Avg edge length'); ylabel('avg edge length / initial avg edge length'); xlabel('Time (hours)');
% ylim([0 5]); yline(1);
% 
% subplot(3,2,2); 
% plot(timex(idx), totnodes(idx), 'Color','blue', 'LineWidth',1); hold on;
% title('Total nodes'); ylabel('# nodes / initial # nodes'); xlabel('Time (hours)');
% ylim([0 20]);
% 
% subplot(3,2,4); 
% plot(timex(idx), nodedegree_avg(idx), 'Color','blue', 'LineWidth',1); hold on;
% title('Avg node degree'); ylabel('avg node degree / initial avg node degree'); xlabel('Time (hours)');
% ylim([0 5]);
% 
% idx2 = convex_areas == 0;
% convex_areas_update = convex_areas(~idx2);
% subplot(3,2,6); 
% plot(timex(idx), convex_areas(idx)/convex_areas_update(1), 'Color','black', 'LineWidth',1); hold on;
% title('Convex area'); ylabel('Convex area'); xlabel('Time (hours)');
% % ylim([0 2]); 
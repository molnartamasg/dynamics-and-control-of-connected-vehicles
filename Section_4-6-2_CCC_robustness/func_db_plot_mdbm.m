

function graphandle=db_plot_mdbm(varargin)
%it polots the output of the db_mdbm function
% db_plot_mdbm(outputsol)
%
% db_plot_mdbm(outputsol,plotcolor)
% plotcolor='r'
% plotcolor='b'
%
% db_plot_mdbm(outputsol,plotcolor,dimensionsorder)
%% e.g.:dimensionsorder=[3,1,2];
%
% db_plot_mdbm(outputsol,plotcolor,dimensionsorder,plotobjdim)
%
% db_plot_mdbm(outputsol,plotcolor,dimensionsorder,plotobjdim,gradplot)
% gradplot=true (1)
% gradplot=false (0)
% varargin{1} - order of plotting
%
% varargin{2}=plotobjdim; it defines the dimension of the plotting object
%   0 - points
%   1- lines
%   2 surfaces


outputsol=varargin{1};


Ndim=size(outputsol.posinterp,1);
Ncodim=length(outputsol.Hval{1});


plotcolor=[];%default plot color and marker
if length(varargin)>1
    if ~isempty(varargin{2})
        plotcolor=varargin{2};
    end
    
end

dimensionsorder=1:(min(3,Ndim));
if length(varargin)>2
    if and(~isempty(varargin{3}),length(varargin{3})<=Ndim)
        dimensionsorder=varargin{3};
        outputsol.posinterp=outputsol.posinterp(dimensionsorder,:);%reorder the points
        outputsol.gradient=outputsol.gradient(:,dimensionsorder,:);
        outputsol.ax=outputsol.ax(dimensionsorder);
        %         Ndim=length(dimensionsorder);
    end
    
end

plotobjdim=min([length(dimensionsorder),length(outputsol.DTalphashapecell),2]);
if length(varargin)>3
    plotobjdim=min([plotobjdim,varargin{4}]);
end
gradplot=0;
if length(varargin)>4
    gradplot=varargin{5};
end


switch Ndim%length(dimensionsorder)
    case 1
        graphandle=plot(outputsol.posinterp,outputsol.posinterp*0,[plotcolor,'.']);
        if gradplot
            hold on
            colorsgrad='krgbymc';
            for k=1:Ncodim
                gradvect=permute(outputsol.gradient(k,:,:),[2,3,1]);
                graphandle(2)=quiver(outputsol.posinterp(1,:),outputsol.posinterp(1,:)*0,...
                    gradvect(1,:),gradvect(1,:)*0,colorsgrad(mod(k-1,7)+1));
                
            end
        end
        axis([min(outputsol.ax(1).val),max(outputsol.ax(1).val),-1,1]) 
        xlabel([num2str(dimensionsorder(1)),'. parameter'])
    case 2
        
        switch plotobjdim
            case 0%points
                graphandle=plot(outputsol.posinterp(1,:),outputsol.posinterp(2,:),[plotcolor,'.']);
            case 1   %lines
                if isempty(plotcolor);
                    graphandle=trimesh(outputsol.DTalphashapecell{1},outputsol.posinterp(1,:),outputsol.posinterp(2,:),outputsol.posinterp(2,:)*0);
                else
                    graphandle=trimesh(outputsol.DTalphashapecell{1},outputsol.posinterp(1,:),outputsol.posinterp(2,:),outputsol.posinterp(2,:)*0,'EdgeColor',plotcolor);
                end
        end
        
        xlabel([num2str(dimensionsorder(1)),'. parameter'])
        ylabel([num2str(dimensionsorder(2)),'. parameter'])
        axis([min(outputsol.ax(1).val),max(outputsol.ax(1).val),min(outputsol.ax(2).val),max(outputsol.ax(2).val)])
        
        if gradplot
            hold on
            colorsgrad='krgbymc';
            for k=1:Ncodim
                gradvect=permute(outputsol.gradient(k,:,:),[2,3,1]);
                graphandle(end+1)=quiver(outputsol.posinterp(1,:),outputsol.posinterp(2,:),...
                    gradvect(1,:),gradvect(2,:),colorsgrad(mod(k-1,7)+1));
                
            end
        end
    case 3
        switch plotobjdim
            case 0 %points
                graphandle=plot3(outputsol.posinterp(1,:),outputsol.posinterp(2,:),outputsol.posinterp(3,:),[plotcolor,'.']);
            case 1 %lines
                if isempty(plotcolor);
                    graphandle=trimesh(outputsol.DTalphashapecell{1},outputsol.posinterp(1,:),outputsol.posinterp(2,:),outputsol.posinterp(3,:));
                else
                    graphandle=trimesh(outputsol.DTalphashapecell{1},outputsol.posinterp(1,:),outputsol.posinterp(2,:),outputsol.posinterp(3,:),'EdgeColor',plotcolor);
                end
            case 2 %surface
                if isempty(plotcolor);
                    graphandle=trisurf(outputsol.DTalphashapecell{2},outputsol.posinterp(1,:),outputsol.posinterp(2,:),outputsol.posinterp(3,:));
                else
                    graphandle=trisurf(outputsol.DTalphashapecell{2},outputsol.posinterp(1,:),outputsol.posinterp(2,:),outputsol.posinterp(3,:),'EdgeColor','none','FaceColor',plotcolor);
                end
        end
        
        

        xlabel([num2str(dimensionsorder(1)),'. parameter'])
        ylabel([num2str(dimensionsorder(2)),'. parameter'])
        zlabel([num2str(dimensionsorder(3)),'. parameter'])
        axis([min(outputsol.ax(1).val),max(outputsol.ax(1).val),min(outputsol.ax(2).val),max(outputsol.ax(2).val),min(outputsol.ax(3).val),max(outputsol.ax(3).val)])
        
        if gradplot
            hold on
            colorsgrad='krgbymc';
            for k=1:Ncodim
                gradvect=permute(outputsol.gradient(k,:,:),[2,3,1]);
                graphandle(end+1)=quiver3(outputsol.posinterp(1,:),outputsol.posinterp(2,:),outputsol.posinterp(3,:),...
                    gradvect(1,:),gradvect(2,:),gradvect(3,:),colorsgrad(mod(k-1,7)+1));
                
            end
        end
        
    otherwise %(4,5,6...) we cannot plot in higer dimensions (or it is pointles)
        switch plotobjdim
            case 0 %points
                graphandle=plot3(outputsol.posinterp(1,:),outputsol.posinterp(2,:),outputsol.posinterp(3,:),[plotcolor,'.']);
            case 1 %lines
                if isempty(plotcolor);
                    graphandle=trimesh(outputsol.DTalphashapecell{1},outputsol.posinterp(1,:),outputsol.posinterp(2,:),outputsol.posinterp(3,:));
                else
                    graphandle=trimesh(outputsol.DTalphashapecell{1},outputsol.posinterp(1,:),outputsol.posinterp(2,:),outputsol.posinterp(3,:),'EdgeColor',plotcolor);
                end
            case 2 %surface
                if isempty(plotcolor);
                    graphandle=trisurf(outputsol.DTalphashapecell{2},outputsol.posinterp(1,:),outputsol.posinterp(2,:),outputsol.posinterp(3,:),outputsol.posinterp(4,:));
                else
                    graphandle=trisurf(outputsol.DTalphashapecell{2},outputsol.posinterp(1,:),outputsol.posinterp(2,:),outputsol.posinterp(3,:),'EdgeColor','none','FaceColor',plotcolor);
                end
                %             graphandle=tetramesh(outputsol.DTalphashapecell{3},outputsol.posinterp(1:4,:)');
        end

              
        
        xlabel([num2str(dimensionsorder(1)),'. parameter'])
        ylabel([num2str(dimensionsorder(2)),'. parameter'])
        zlabel([num2str(dimensionsorder(3)),'. parameter'])
        axis([min(outputsol.ax(1).val),max(outputsol.ax(1).val),min(outputsol.ax(2).val),max(outputsol.ax(2).val),min(outputsol.ax(3).val),max(outputsol.ax(3).val)])
        
        if gradplot
            hold on
            colorsgrad='krgbymc';
            for k=1:Ncodim
                gradvect=permute(outputsol.gradient(k,:,:),[2,3,1]);
                graphandle(end+1)=quiver3(outputsol.posinterp(1,:),outputsol.posinterp(2,:),outputsol.posinterp(3,:),...
                    gradvect(1,:),gradvect(2,:),gradvect(3,:),colorsgrad(mod(k-1,7)+1));
                
            end
        end
        
end

grid on


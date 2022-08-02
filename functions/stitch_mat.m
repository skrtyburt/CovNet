function stitch_mat(Lmat,Umat,cmap,boxColor,tri_labels,node_labels,clims)

narginchk(2,8)

[r1,c1] = size(Lmat);
[r2,c2] = size(Umat);

if r1~=c1 || r2~=c2
    fprintf(2,'Asymmetrics adjacency matrix inputs! Exiting...\n');
    return
end
if r1~=r2
    fprintf(2,'Unequal size matrices! Exiting...\n');
    return
end
Lmat_symm = (Lmat+Lmat')/2;
Umat_symm = (Umat+Umat')/2;
if ~isequal(Lmat,Lmat_symm) || ~isequal(Umat,Umat_symm)
    fprintf(2,'One or both of the input matrices not symmetric!\n')
    Lres = max(max(Lmat_symm - Lmat));
    Ures = max(max(Umat_symm - Umat));
    if Lres > .01 || Ures > .01
        fprintf(2,'Asymmetry seems greater than rounding difference.\n')
        fprintf(2,'Check that matrices are correct/symmetric. Exiting...\n')
        return
    else
        fprintf('Asymmetrry <0.01; Lileky rounding difference.\n')
        fprintf('Continuting with forced symmetric matrices.\n')
        Lmat = Lmat_symm;
        Umat = Umat_symm;
    end
end
if ~exist('cmap','var')
    cmap = turbo(64);
end
if ~exist('boxColor','var')
    boxColor = [0 0 0];
end
if length(node_labels)~=r1
    fprintf(2,'node_labels different length than matrix.\n')
    clear node_labels
end


Ne = r1 + .5;
Ze = 0.5;
bx = [Ze Ze; Ze Ne; Ne Ne; Ze Ze; Ne Ze; Ne Ne];

stitchmat = tril(Lmat,-1) + triu(Umat,1);
if ~exist('clims','var')
    mn = min(stitchmat(:));
    mx = max(stitchmat(:));
else
    mn = clims(1);
    mx = clims(2);
end

f = figure(...
        'units','inches',...
        'position',[1 1 7.5 6],...
        'paperpositionmode','auto');

imagesc(stitchmat,[mn mx]); axis square
colormap(cmap); colorbar
box off
a=gca;
a.XAxis.TickLength = [0 0];
a.YAxis.TickLength = [0 0];
hold on
line(bx(:,1),bx(:,2),'Color',boxColor,'LineWidth',2)


if exist('node_labels','var')
    xticks(1:1:r1)
    xticklabels(node_labels)
    xtickangle(45)
    yticks(1:1:r1)
    yticklabels(node_labels)
end
if exist('tri_labels','var')
    ylabel(tri_labels{1})
    yyaxis('right')
    ylabel(tri_labels{2},'Rotation',270,'Position',[Ne+1.5 .5 0],'Color','k')
    yticks([])
end





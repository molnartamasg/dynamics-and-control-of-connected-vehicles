function point=p_remesh(point,new_degree,new_mesh)

% function rm_point=p_remesh(point,new_degree,new_mesh)
% INPUT:
%	point periodic solution point
%	new_degree new degree for new representation
%	new_mesh new mesh or new number of intervals
% OUTPUT:
%	rm_point interpolated point on new mesh

% (c) DDE-BIFTOOL v. 2.00, 23/11/2001

switch point.kind,
    case 'hcli',
        if length(new_mesh)>1
            point.profile=hcli_eva(point.profile,point.mesh,new_mesh,point.degree);
            point.mesh=new_mesh;
            point.degree=new_degree;
            l=(length(new_mesh)-1)/new_degree;
            if l~=floor(l)
                error('P_REMESH: length of new mesh %d and new degree %d do not match.',...
                    length(new_mesh),new_degree);
            end;
        else
            if isempty(point.mesh)
                mesh=0:1/(size(point.profile,2)-1):1;
                t_new=psol_msh(mesh,point.degree,point.profile,new_mesh,new_degree);
                point.profile=hcli_eva(point.profile,mesh,t_new,point.degree);
            else
                t_new=psol_msh(point.mesh,point.degree,point.profile,new_mesh,new_degree);
                point.profile=hcli_eva(point.profile,point.mesh,t_new,point.degree);
            end;
            point.degree=new_degree;
            point.mesh=t_new;
        end;
    case 'psol',
        if length(new_mesh)>1
            point.profile=psol_eva(point.profile,point.mesh,new_mesh,point.degree);
            point.mesh=new_mesh;
            point.degree=new_degree;
            l=(length(new_mesh)-1)/new_degree;
            if l~=floor(l)
                error('P_REMESH: length of new mesh %d and new degree %d do not match.',...
                    length(new_mesh),new_degree);
            end;
        else
            if isempty(point.mesh)
                mesh=0:1/(size(point.profile,2)-1):1;
                t_new=psol_msh(mesh,point.degree,point.profile,new_mesh,new_degree);
                point.profile=psol_eva(point.profile,mesh,t_new,point.degree);
            else
                t_new=psol_msh(point.mesh,point.degree,point.profile,new_mesh,new_degree);
                point.profile=psol_eva(point.profile,point.mesh,t_new,point.degree);
            end;
            point.degree=new_degree;
            point.mesh=t_new;
        end;
        if isfield(point,'stability');
            point.stability=[];
        end;
    otherwise,
        error('P_REMESH: no periodic solution point but %s.',point.kind);
end;
end

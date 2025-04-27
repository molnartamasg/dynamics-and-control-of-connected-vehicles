function s=p_measur(p,m)

% function sc_measure=p_measur(p,measure)
% INPUT:
%	p point
%	measure measure struct
% OUTPUT:
%	sc_measure scalar measure of point

% (c) DDE-BIFTOOL v. 1.00, 06/05/2000

row=m.row;
col=m.col;

if length(row)==3
  row(4)=' ';
end;
if length(col)==3
  col(4)=' ';
end;

f=getfield(p,m.field);

if isempty(f)
  s=[];
  return;
end;

if ~isempty(m.subfield)
  f=getfield(f,m.subfield);
  if isempty(f)
    s=[];
    return;
  end;
end;

if col=='max '
  s=max(f(row,:));
elseif col=='min '
  s=min(f(row,:));
elseif col=='mean'
  s=sum(f(row,:))/length(f(row,:));
elseif col=='ampl'
  s=max(f(row,:))-min(f(row,:));
elseif row=='max '
  s=max(f(:,col));
elseif row=='min '
  s=min(f(:,col));
elseif row=='mean'
  s=sum(f(:,col))/length(f(:,col));
elseif row=='ampl'
  s=max(f(:,col))-min(f(:,col));
elseif row=='all '
  if col=='all '
    s=f(:,:);
  else
    s=f(:,col);
  end;
elseif col=='all '
  s=f(row,:);
else
  s=f(row,col);
end;

if ~isempty(m.func)
  s=feval(m.func,s);
end;

return;

%% se 1 solo input, normalizza tra minimo e massimo, se 3, normalizza tra p1-esimo percentile e il p2-esimo

function v2=norma(v1,p1,p2)
if nargin==1
    v2=(v1-repmat(nanmin(v1),size(v1,1),1))./...
    (repmat(nanmax(v1),size(v1,1),1)-repmat(nanmin(v1),size(v1,1),1));
elseif nargin==3
     v2=(v1-repmat(prctile(v1,p1),size(v1,1),1))./...
    (repmat(prctile(v1,p2),size(v1,1),1)-repmat(prctile(v1,p1),size(v1,1),1));
    v2(v2<0)=NaN;
    v2(v2>1)=NaN;
end


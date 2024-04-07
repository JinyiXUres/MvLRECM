function [ weight ] = mvlrecm_weight( data, mvU, mveudis2, S, alpha, delta2, eta )
view_n = length(data);
data_n=size(data{1},1);
c=sum(S,2);
c_rep=repmat((c(1:end-1,:).^alpha),1,data_n);
Phi=zeros(1,view_n); % ¡ø
for m=1:view_n
    value1=c_rep .* (mvU{m}(1:end-1,:).^2) .* mveudis2{m}(1:end-1,:);
    value2=delta2*sum(mvU{m}(end,:).^2);
    Phi(m)=(sum(sum(value1))+value2);
end
t_zoom=exp(sum(floor(log(Phi)))/length(Phi));% Phi is too large to calculate e^Phi
Phi_zoom=Phi./t_zoom;
eta_zoom=eta/t_zoom;
weight=zeros(1,view_n);
e=exp(1);
for m=1:view_n
    value1=e^(-(eta_zoom+Phi_zoom(m))/eta);
    value2=sum(e.^(-(eta_zoom+Phi_zoom)./eta));
    weight(m)=value1/value2;
end

end


function out = Zero_Intercept_Times(data, timeVector)
% this function assumes a relitavely well behaved perodic function like a
% sine curve, monotonic data near zero crossings.  It returns the linearly
% interpolated time of the zero crossing.

if size(timeVector,2)>size(timeVector,1)
    timeVector = timeVector';
end

dd = (data >0);
ee = dd(2:end) - dd(1:end-1);
gg = (ee ~= 0);
cross_indx = find(gg == 1);

pair_ind = [cross_indx cross_indx+1];

out = (-data(pair_ind(:,1)).*(timeVector(pair_ind(:,2)) - timeVector(pair_ind(:,1))))./ ...
      (data(pair_ind(:,2)) -data(pair_ind(:,1))) + timeVector(pair_ind(:,1));
  
  if data(pair_ind(1,1)) < 0
      out(1) = [];
  end

end

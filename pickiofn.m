function chosen = pickiofn(i,n)
% chosen = pickiofn(i,n);
% Constructs a random logical n-vector with exactly i truths, and 
% all combinations equally likely.

reverse = i > n/2;
if (reverse)
    i = n-i;
end

chosen = false(1,n);

for p = 1:i
   m = n-p+1;     % number of possibilities on this roll of the dice
   ithis = ceil(m*rand);
   
   % mark the next one not taken yet
   % This loop must terminate because there cannot be more than p-1
   % true values in the last (n-m+1) = p elements, and ithis <= m.
   
   while chosen(ithis)
       ithis = ithis + 1;
   end
   chosen(ithis) = true;
end

if reverse, chosen = ~chosen; end

end

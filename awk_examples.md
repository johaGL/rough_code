add +1 to third column : 
```
cd ldh_multi/mulitplex/coex/
awk -v s=1 '{print $1, $2, $3+s}' OFS='\t' edM10M_old_zeroesin.tsv > ed10M_sgControlPOS.tsv
```



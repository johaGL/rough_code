Add +1 to third column : 
```
cd ldh_multi/mulitplex/coex/
awk -v s=1 '{print $1, $2, $3+s}' OFS='\t' edM10M_old_zeroesin.tsv > ed10M_sgControlPOS.tsv
```


Multiply by -1 sixth column, sparing header (two steps, not very smart but works):
note: input is tab separated (if FS not set, espace is also took as separator ! )
```
mkdir test

head -n 1 atoy > test/atoy

awk -v FS="\t" 'NR>1 {print $1,$2,$3,$4,$5, $6 * -1, $7,$8,$9,$10}' OFS="\t" atoy >> test/atoy

```
I used it like this : cd $tableslocation; for i in *.tsv; do head -n 1 $i > test/$i; awk -v FS="\t" 'NR>1 {print $1,$2,$3,$4,$5, $6 * -1, $7,$8,$9,$10}' OFS="\t" $i >> test/$i ;done   . That yielded files named identically but in test subfolder, then replaced parent content:  mv test/\* . ;rm test.


More serious examples with explanations : [here](https://www.theunixschool.com/2012/11/awk-examples-insert-remove-update-fields.html)


#### simple but important
Rights for accession and execution, all my group:
chmod g+r /home/jgalvis
chmod g+rx /home/jgalvis

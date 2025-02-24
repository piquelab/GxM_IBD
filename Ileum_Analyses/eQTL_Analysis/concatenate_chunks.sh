for i in $(seq 1 30); do
for j in $(seq 1 30); do
     zcat  output/PC1-$i.permutations.chunk$j.txt.gz
done | gzip -c >  output/PC1-$i.permutations.eQTL.txt.gz;
done

for j in $(seq 1 30); do
     zcat  output/PC1-0.permutations.chunk$j.txt.gz
done | gzip -c >  output/PC1-0.permutations.eQTL.txt.gz
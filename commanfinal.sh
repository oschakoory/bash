#command sur subset


#parallel 'bzcat {1}.fastq.bz2 | split -a 3 -l '400000' /dev/stdin splitted.{1.} ; parallel -I %% "bzip2 -f %%" ::: splitted.{1.}???' ::: *.fastq.bz2

parallel 'cat {1} | statsFastQ --label={1}' ::: *.subset.fastq.bz2 > statisticsFastq.raw.tsv

|parallel 'cat {1}| fastqTagCleaner -64 | fastqCutLongAtN | fastqPolyATendTrim --repeatLimit=6 | fastqComplexity --complexity=30 | fastqSelectLength --min=30 | bzip2 -c > clean.{1}' ::: *.subset.fastq.bz2

|parallel 'bzcat {1} | statsFastQ --label={1}' ::: clean.* > statisticsFastq.clean.tsv

|parallel 'bwa samse All_reference.fasta <(bwa aln -t 8 All_reference.fasta <(bzcat {1}) ) <(bzcat {1})| samtools view -S -f 4 /dev/stdin | awk "{print \"@\"\$1\"\n\"\$10\"\n+\n\"\$11}"| bzip2 -c > decont.{1}' ::: clean.*

|parallel 'bzcat {1} | statsFastQ --label={1}' ::: decont.* > statisticsFastq.decont.tsv


|parallel --xapply 'fastqPairEndSynchronizer -s -p1={1} -p2={2}' ::: decont.*_1.subset.fastq.bz2 ::: decont.*_2.subset.fastq.bz2

|parallel 'bzcat {1} | statsFastQ --label={1}' ::: paired.decont.clean.* > statisticsFastq.paired.tsv

|parallel 'bwa mem popGenome.73013transcripts.fasta \
                                     <(bzcat {1}) \
                                     <(bzcat {2}) \
                                | samtools view -bS /dev/stdin \
                                | samtools sort -m 8G /dev/stdin \
                                > mapped.{1}.bam' \
        ::: paired.decont.clean.*_1.subset.fastq.bz2 \
        ::: paired.decont.clean.*_2.subset.fastq.bz2

#|ls paired.decont.{1}_1.fastq???.bz2 | sed -r 's:\.[^.]+$::;s:\.[^.]+$::;' | sort | uniq > targetSplittedSamples

#|parallel --xapply 'samtools merge -f -r total.mapped.{1}.bam $(ls mapped.{1}.fastq???.bz2.bam)' ::: targetSplittedSamples

#|rm mapped.{1}.bam

|parallel --xapply 'samtools index {1}' ::: mapped.paired.decont.clean.*.bam

|parallel 'samtools idxstats {1} > resultsCounts.{= s/^\.[^.]//; =}.txt' ::: mapped.paired.decont.clean.*.bam


|paste  <(echo -e "TranscriptID\tTranscriptBp" ; cut -f1,2 *.bam.txt) <(cat   <(for NAME in *.bam.txt ; do NAME=$(basename seq) ; echo -ne "${NAME%_*}.map\t${NAME%_*}.bkg\t" ; done ; echo ) <( parallel -j8 -k 'cut -f3- {} > {}.cut' ::: *.bam.txt ; paste *.bam.txt.cut )) > Score.csv


########################################################################################################################################

#command sur subset et split

parallel 'bzcat {1} | head -n1600000 | bzip2 -c > {1}.subset.fastq.bz2' ::: *.fastq.bz2

|parallel 'bzcat {1} | statsFastQ --label={1}' ::: splitted.* > statisticsFastq.Unsplitfinal.tsv

|parallel 'bzcat {1}| split -a 3 -l '400000' /dev/stdin splitted.{1.} ; parallel -I %% "bzip2 -f %%" ::: splitted.{1.}???' ::: *.subset.fastq.bz2

|parallel 'bzcat {1} | statsFastQ --label={1}' ::: splitted.* > statisticsFastq.raw.tsv

|parallel 'bzcat {1}| fastqTagCleaner -64 | fastqCutLongAtN | fastqPolyATendTrim --repeatLimit=6 | fastqComplexity --complexity=30 | fastqSelectLength --min=30 | bzip2 -c > clean.{1}' ::: splitted.*

|parallel 'bzcat {1} | statsFastQ --label={1}' ::: clean.* > statisticsFastq.clean.tsv

|parallel 'bwa samse All_reference.fasta <(bwa aln -t 8 All_reference.fasta <(bzcat {1}) ) <(bzcat {1})| samtools view -S -f 4 /dev/stdin | awk "{print \"@\"\$1\"\n\"\$10\"\n+\n\"\$11}"| bzip2 -c > decont.{1}' ::: clean.*

|rm clean.*

|parallel 'bzcat {1} | statsFastQ --label={1}' ::: decont.* > statisticsFastq.decont.tsv


#|parallel --xapply 'fastqPairEndSynchronizer -s -p1={1} -p2={2}' ::: decont.*_1.fastq.bz2.subset.fastq???.bz2 ::: decont.*_2.fastq.bz2.subset.fastq???.bz2
|parallel --xapply 'fastqPairEndSynchronizer -s -p1={1} -p2={2}' ::: decont.*_1.* ::: decont.*_2*
|rm decont.*

|parallel 'bzcat {1} | statsFastQ --label={1}' ::: paired.decont.clean.* > statisticsFastq.paired.tsv

|parallel 'bwa mem popGenome.73013transcripts.fasta \
                                     <(bzcat {1}) \
                                     <(bzcat {2}) \
                                | samtools view -bS /dev/stdin \
                                | samtools sort -m 8G /dev/stdin \
                                > mapped.{1}.bam' \
        ::: paired.decont.clean.*_1* \
        ::: paired.decont.clean.*_2*

|rm paired.decont.*

|toto=$(ls mapped.*_1.* | sed -r 's:\.[^.]+$:: ; s:\.[^.]+$:: ; s:\.[^.]+$:: ;' | sort | uniq)

|parallel 'samtools merge -f -r total.{1}.bam $(ls *.fastq???.bz2.bam)' ::: $toto


|parallel 'samtools index {1}' ::: total.*

|parallel 'samtools idxstats {1} > resultsCounts.{= s/^\.[^.]//; =}.txt' ::: total.*.bam

#error present
|paste  <(echo -e "TranscriptID\tTranscriptBp" ; cut -f1,2 *.bam.txt) <(cat   <(for NAME in *.bam.txt ; do NAME=$(basename seq) ; echo -ne "${NAME%_*}.map\t${NAME%_*}.bkg\t" ; done ; echo ) <( parallel -j8 -k 'cut -f3- {} > {}.cut' ::: *.bam.txt ; paste *.bam.txt.cut )) > Score.csv

|paste  <(echo -e "TranscriptID\tTranscriptBp" ; cut -f1,2 *.bam.txt) <(cat   <(while NAME in *.bam.txt ; do NAME=$(basename seq) ; echo -ne "${NAME%_*}.map\t${NAME%_*}.bkg\t" ; done ; echo ) <( parallel -j8 -k 'cut -f3- {} > {}.cut' ::: *.bam.txt ; paste *.bam.txt.cut )) > Score.csv

###########################################################################################################################################"


#command sur tous les fichiers

parallel 'bzcat {1} | statsFastQ --label={1}' ::: *.fastq.bz2 > statisticsFastq.Unsplitfinal.tsv

|parallel 'bzcat {1} | split -a 3 -l '16 000 000' /dev/stdin splitted.{1.} ; parallel -I %% "bzip2 -f %%" ::: splitted.{1.}???' ::: *.fastq.bz2

| parallel 'bzcat {1} | statsFastQ --label={1}' ::: splitted.* > statisticsFastq.Finalraw.tsv

| parallel 'bzcat {1} | fastqTagCleaner -64 | fastqCutLongAtN | fastqPolyATendTrim --repeatLimit=6 | fastqComplexity --complexity=30 | fastqSelectLength --min=30 | bzip2 -c > clean.{1}' ::: splitted.*

| parallel 'bzcat {1} | statsFastQ --label={1}' ::: clean.* > statisticsFastq.Finalclean.tsv

| parallel 'bwa samse All_reference.fasta \ 		
                            <(bwa aln -t 8 All_reference.fasta <(bzcat {1}) ) \ 
                            <(bzcat {1}) \				
                        | samtools view -S -f 4 /dev/stdin \
                        | awk "{print \"@\"\$1\"\n\"\$10\"\n+\n\"\$11}" \
                        | bzip2 -c \
                        > decont.{1}' \
     ::: clean.*

| rm clean.*

| parallel 'bzcat {1} | statsFastQ --label={1}' ::: decont.* > statisticsFastq.Finaldecont.tsv

| parallel --xapply 'fastqPairEndSynchronizer -s -p1={1} -p2={2}' :::decont.*_1.* :::decont.*_2.*

| rm decont.*

| parallel 'bzcat {1} | statsFastQ --label={1}' ::: paired.* > statisticsFastq.paired.tsv

| parallel 'bwa mem popGenome.73013transcripts.fasta \
                                     <(bzcat {1}) \
                                     <(bzcat {2}) \
                                | samtools view -bS /dev/stdin \
                                | samtools sort -m 8G /dev/stdin \
                                > mapped.{1}.bam' \
        ::: paired.decont.*_1.* \
        ::: paired.decont.*_2.*

|rm paired.*

| scriptfinal=$(ls mapped.* | sed -r 's:\.[^.]+$::;s:\.[^.]+$::;' | sort | uniq) 

|parallel 'samtools merge -f -r total.{1}.bam $(ls *.fastq???.bz2.bam)' ::: ${scriptfinal}

rm mapped.*.bam

| parallel 'samtools index total.mapped.{1}.bam' ::: ${scriptfinal}

| parallel 'samtools idxstats total.mapped.{1}.bam > resultsCounts.{= s/^\.[^.]//; =}.txt' ::: ${scriptfinal}

| exit 







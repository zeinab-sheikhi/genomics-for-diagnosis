Reference indexing is the process of creating search indexes for a reference genome to enable fast sequence alignment.

Output: 
.bwt - Main search index for finding read locations
.sa - Maps index positions back to genome coordinates
.pac - Provides the actual reference sequence
.ann - Tells BWA chromosome names and lengths
.amb - Handles ambiguous bases (N's) properly

QC observation: "Unusual GC content detected in R1"
Hypothesis: "May indicate structural variants affecting GC distribution"

INDEXES="ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/alignment.index "
INDEXES="$INDEXES ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/sequence.index"

for i in $INDEXES ; do
	wget $i
done

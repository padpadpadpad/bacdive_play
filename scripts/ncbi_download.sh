# download genome assemblies from NCBI

conda env list

# create ncbigenomedownload
conda create -n ncbi-download -c bioconda ncbi-genome-download

conda activate ncbi-download

ncbi-genome-download -help

# run one for genbank
ncbi-genome-download -F "assembly-stats" -A data/bacdive_ncbi_accession.txt -o data/ncbi_download -s "genbank" bacteria

# run one for refseq
ncbi-genome-download -F " assembly-stats" -A data/bacdive_ncbi_accession.txt -o data/ncbi_download -s "refseq" bacteria
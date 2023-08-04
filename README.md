# Analysis for publication "High-throughput parallel DNA RNA sequencing of single nuclei"


The analysis for the manuscript was performed in the following steps:

1) Data preprocessing with cellranger pipelines

For HIPSD-seq and sciHIPSD-seq  we used standard cellranger ATAC pipeline and for RNA from HIPSD&R-seq we use standard cellranger RNA pipeline

2. Cells definition
    - RNA part

        The cells were defined by cellranger pipeline.
    - DNA part

        The stadard cellranger ATAC pipeline can not define cells from HIPSD-seq output as it bases the cell selection based on peak enrichment. Therefore, selection based on total counts is used. Example of this step is given bellow:
        ```python
        filename = "atac_fragments.tsv.gz"
        read_counter = dict()
        chunksize = 10 ** 6
        with pd.read_csv(filename, chunksize=chunksize, comment = "#", sep = "\t", header = None) as reader:
            for chunk in tqdm(reader):
                for tup in chunk.itertuples():
                    if tup[4] not in read_counter:
                        read_counter[tup[4]] = tup[5]
                    else:
                        read_counter[tup[4]] = read_counter[tup[4]] + tup[5]
        counts = np.array(list(read_counter.values()))
        hist, bins = np.histogram(counts, bins=np.arange(1, np.max(counts)+2))
        cumulative_hist = np.cumsum(hist[::-1])[::-1]
        
        # Create the elbow plot
        sns.scatterplot(y=bins[:-1], x = cumulative_hist)
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel("Barcodes")
        plt.ylabel("Counts")
        plt.legend()
        plt.show()

        # select cut-off based on the elbow plot
        cutoff = 10000
        valid_barcodes = [x for x in read_counter if read_counter[x]>=cutoff]

        # save the results
        bfile = data_path + "barcodes.tsv"
        with open(bfile, "w") as out:
            for cell in valid_barcodes:
                out.write(cell + "\n")


        ``` 

2) Counts creation at 100kb resolution

The code snippet above can be used to have a list of cell barcodes that is needed to extract cell BAM files from the BAM file produced by cellranger. For more details on extraction of individual BAM files, please refer to the original cellranger [software](https://github.com/10XGenomics/subset-bam). Each cell BAM file should also should be filetered, sorted and indexed. The code snippet below can be used for this purpose:

```bash
samtools view -F 3844 -q 30 "$CELL_BAM" | samtools view -bS - | samtools sort -o "$CLEAN_SORT_BAM" -
samtools index "$CLEAN_SORT_BAM"
```
To extract counts from each cell BAM file, we used [hmmcopy_utils](https://github.com/shahcompbio/hmmcopy_utils). E.g:
```bash
    readCounter -q 30 -w 100000 -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY $CLEAN_SORT_RH > $CELL_COUNT_FILE

```     
Please note that for step 5, one would also need genome mappability and gc content files. These files can also be generated using [hmmcopy_utils](https://github.com/shahcompbio/hmmcopy_utils)

3) Metacells creation

Metacells are produced based on a greedy algorithm. Please see the full notebook `metacell_creation.ipynb` for details. To execute the notebook, one would need to have the following files:

- readcounts for each cell (output of `readCounter` command)

- `create_metacells.py` script

The result of the metacelling is a directory with `.txt` files. In each file there is a path to readcount files that are assigned to the metacell.

4) Copy number calling

Copy number colling is performed using [hmmcopy](https://bioconductor.org/packages/release/bioc/html/HMMcopy.html) R package. 

The original documentation can be used to reproduce the results. We also developed a wrapper script that can be used to reprocduce the analysis precisly. The script is called `run_hmmcopy.R` and it can be used as follows:

```bash
for file in metacell_files/*.txt; do  Rscript run_hmmcopy.R -f $file -o hmmcopy_results/ -r ../reference/hmmcopy_references/cellranger_xy_ -l -b  -m -3.5; done
```
To view all available options, please run `Rscript run_hmmcopy.R -h`.
The script produces `.bed`  files for each cell at 100kb resolution.




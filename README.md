# Simulating WGS cfDNA with scWGS
A python package that simulates cfDNA data by mixing scWGS data. 

## Getting started 

1. Create a working directory:
   ```
   mkdir -p path/to/project
   cd path/to/project
   ```
2. Clone the repository through git:
    ```
    git clone --depth 1 https://github.com/matteolepur/cfsim.git
    ```
3. Install
    ```
    cd cfsim 
    pip install .
    ```
4. Verify the installation worked:
    ```
    cfsim
    ```
--------
## Simulating cfDNA 

```bash
simulate \
  -d <hapclone-data-file> \
  -r <hapclone-results-file> \
  -s <snp-file> \
  -o <out-file> \
  --clone-prevalence-file <clone-prevalence-file> \
  -c <coverage> \
  -p <clone-prevalence-prior> \
  -t <tumour-content> \
  --read-length <read-length> \
  --seed <seed>
```


| Flag | Description |
|---|---|
| `--hapclone-data-file` | Path to single cell counts data file h5 file. |
| `--hapclone-results-file` | Path to the single cell clone cn results tsv file. |
| `--snp-file` | Path to the SNP file. |
| `--clone-prevalence-file` | Path to the file specifying clone prevalence. |
| `--coverage` | Target sequencing coverage for the simulation. |
| `--clone-prevalence-prior` | Prior value for clone prevalence. |
| `--tumour-content` | Fraction of tumour content to simulate. |
| `--read-length` | Length of simulated reads. |
| `--seed` | Random seed for reproducibility. |
| `--out-file` | Destination path where the simulation output is written. |

**Example inputs and output**

`--hapclone-data-file` (see [DataSet](src/cfsim/dataset.py) class).

```
data = DataSet.from_file(arg.hapclone_data_file)

data.normal_reads
# array([ 12,  15,   9, ...,  11,  14,  10])                                                # shape (num_bins,)

data.reads
# array([[[12], [15], [ 9], ..., [11], [14], [10]],
#        [[ 8], [10], [ 7], ..., [ 9], [12], [ 8]],
#        ...,
#        [[14], [16], [11], ..., [13], [15], [12]]])                                        # shape (num_cells, num_bins, 1)

data.cells
# array(['cell_0001', 'cell_0002', ..., 'cell_0499', 'cell_0500'])                          # shape (num_cells,)

data.bins
# array(['chr1:0:500000', 'chr1:500000:1000000', ..., 'chrY:56500000:57000000'])            # shape (num_bins,)
```

`--hapclone-results-file`
```
cell_id	cluster_id	chrom	beg	end	a	b	a_cor	b_cor	rdr	rdr_scaled	reads	cn_A	cn_B	cn_A_cell	cn_B_cell	cell_seg	seg
TFRIPAIR4_FL_A98167_R04-C08	0	chr1	0	500000	0	0	0	0			8	2	0	1	1	1	1
TFRIPAIR4_FL_A98167_R04-C08	0	chr1	500000	1000000	0	0	0	0			20	2	0	1	1	1	1
TFRIPAIR4_FL_A98167_R04-C08	0	chr1	1000000	1500000	1	0	1	0	0.5900735864607976	1.204670758098118	39	2	0	1	1	1	1
TFRIPAIR4_FL_A98167_R04-C08	0	chr1	1500000	2000000	3	2	3	2	0.5900735864607976	1.204670758098118	41	0	1	1	1	1	2
TFRIPAIR4_FL_A98167_R04-C08	0	chr1	2000000	2500000	0	2	0	2	0.7175930920078096	1.4650095072717724	60	0	1	0	2	2	2
...
```

`--out-file`
```
chrom	start	end	reads	a	b	gc	map	valid	gc_correction	rdr
chr1	1000000	1500000	1565	72	66	0.605014	0.952257	True	1619.0143652063894	0.9666375009591075
chr1	1500000	2000000	1465	127	126	0.540692	0.905416	True	1653.072475793659	0.8862285359246798
chr1	2000000	2500000	1532	78	90	0.585454	0.985942	True	1628.8321623323736	0.940551172446327
chr1	3000000	3500000	1549	57	70	0.573834	0.988141	True	1634.887704895853	0.9474656854787931
chr1	3500000	4000000	1601	135	160	0.555584	0.988141	True	1644.734001033254	0.9734096814404157
chr1	4000000	4500000	1555	160	211	0.479274	0.975042	True	1690.348083044838	0.9199288688510627
...
```
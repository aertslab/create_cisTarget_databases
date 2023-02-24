# Create cisTarget databases

> :information_source: **The SCENIC+ motif collection PWMs will be made available upon publication.**
>
> You can find precomputed databases for human, mouse and fly at https://resources.aertslab.org/cistarget/.
> If you want early access to the PWMs, you can send us an email to carmen.bravogonzalezblas@kuleuven.be,
> seppe.dewinter@kuleuven.be and/or stein.aerts@kuleuven.be with your motivation to request the collection
> and we may decide to share it privately with you.


## Installation


### Clone `create_cisTarget_databases` source code

```bash
# Clone git repo.
git clone https://github.com/aertslab/create_cisTarget_databases

cd create_cisTarget_databases

# Display to which value ${create_cistarget_databases_dir} variable should be set.
echo "create_cistarget_databases_dir='""${PWD}""'"
```


### Create conda environment

```bash
# Create conda environment.
conda create -n create_cistarget_databases \
    'python=3.10' \
    'numpy=1.21' \
    'pandas>=1.4.1' \
    'pyarrow>=7.0.0' \
    'numba>=0.55.1' \
    'python-flatbuffers'
```

### Install Cluster-Buster

Install [Cluster-Buster](https://github.com/ghuls/cluster-buster/) for scoring regulatory regions with motifs.

#### Install precompiled binary

```bash
# Activate conda environment.
conda activate create_cistarget_databases

cd "${CONDA_PREFIX}/bin"

# Download precompiled Cluster-Buster binary.
wget https://resources.aertslab.org/cistarget/programs/cbust

# Make downloaded binary executable.
chmod a+x cbust
```

#### Compile from source

```bash
# Clone Cluster-Buster repo.
#git clone https://github.com/weng-lab/cluster-buster/
git clone -b change_f4_output https://github.com/ghuls/cluster-buster/

cd cluster-buster

# Compile Cluster-Buster.
make cbust

# Compile Cluster-Buster with AMD Math Library (LibM).
#
# This can be much faster (like twice as fast) than using the
# glibc math library on older distributions.
#
#   https://developer.amd.com/amd-aocl/amd-math-library-libm/
#
#     mkdir ../aocl-libm
#     cd ../aocl-libm
#
#     AMD_LIB_VERSION=3.1.0
#
#     tar xzf aocl-libm-linux-aocc-${AMD_LIB_VERSION}.tar.gz
#     mv amd-libm amd-libm-aocc
make cbust_amd_libm_aocc

# Activate conda environment.
conda activate create_cistarget_databases

# Copy CLuster-Buster binary of your choice in conda environment.
cp -a cbust "${CONDA_PREFIX}/bin/cbust"
cp -a cbust_amd_libm_aocc "${CONDA_PREFIX}/bin/cbust"
```


### Install UCSC tools

Install some UCSC tools:
  - `liftOver`: Move regulatory regions from one assembly to another. Needed when generating regulatory regions in
    other species when creating cross-species regulatory regions.
  - `bigWigAverageOverBed`: Compute average/max score of bigWig files for each regulatory region. Used for scoring TF
    ChIP-seq bigWig files.

```bash
# Activate conda environment.
conda activate create_cistarget_databases

cd "${CONDA_PREFIX}/bin"

# Download liftOver.
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver

# Download bigWigAverageOverBed.
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed

# Make downloaded binaries executable.
chmod a+x liftOver bigWigAverageOverBed
```


### Activate environment

Before running any of the scripts, load `create_cistarget_databases` conda environment
and set the `create_cistarget_databases_dir` variable to the dir that contains the cloned repo.

```bash
# Activate conda environment.
conda activate create_cistarget_databases

# Set ${create_cistarget_databases_dir} variable to path where the repo was cloned to.
create_cistarget_databases_dir=""
```


## Memory requirements

Creating cisTarget databases can be a very memory intensive job as it needs to create/store a 2D matrix with dimensions
(number of motifs/tracks vs number of regions/genes) or vice versa. Besides this, it needs (relatively little) memory
to store motif/tracks and regions/genes names.

Memory size of **cisTarget scores database** when loaded in memory:
  - `4 bytes x number of regions/genes x number of motifs/tracks`
  -  memory needed to store region/genes names and motifs/tracks names

Memory size of **cisTarget rankings database** when loaded in memory:
  - 32768 regions/genes or less:
      - `2 bytes x number of regions/genes x number of motifs/tracks`
      - memory needed to store region/genes names and motifs/tracks names
  - more than 32768 regions/genes:
      - `4 bytes x number of regions/genes x number of motifs/tracks`
      -  memory needed to store region/genes names and motifs/tracks names

*Examples:*

| cisTarget database type | byte size (1 element) | #genes/regions | #motifs/tracks | RAM requirement for 2D matrix |
| ---: | :---: | :---: | :---: | :---: |
| *cisTarget scores database*   | 4 bytes | 20000 genes     | 10000 motifs | `4 bytes x 20000 x 10000 = 0.745 GB` |
| *cisTarget rankings database* | 2 bytes | 20000 genes     | 10000 motifs | `2 bytes x 20000 x 10000 = 0.373 GB` |
| *cisTarget scores database*   | 4 bytes | 1000000 regions | 10000 motifs | `4 bytes x 1000000 x 10000 = 37.253 GB` |
| *cisTarget rankings database* | 4 bytes | 1000000 regions | 10000 motifs | `4 bytes x 1000000 x 10000 = 37.253 GB` |
| *cisTarget scores database*   | 4 bytes | 1000000 regions | 1000 tracks  | `4 bytes x 1000000 x 1000 = 3.725 GB` |
| *cisTarget rankings database* | 4 bytes | 1000000 regions | 1000 tracks  | `4 bytes x 1000000 x 1000 = 3.725 GB` |


When running the scripts in this repo, you might need around 3 times the amount of RAM of the actual database.



## Scripts overview

|                                                                                                                                                                   script | description                                                                                                                                                                                                                                                |
|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|                                                                                             [`create_cistarget_motif_databases.py`](#create_cistarget_motif_databasespy) | Create cisTarget motif databases.                                                                                                                                                                                                                          |
|                                                                                             [`create_cistarget_track_databases.py`](#create_cistarget_track_databasespy) | Create cisTarget track databases.                                                                                                                                                                                                                          |
|           [`combine_partial_motifs_or_tracks_vs_regions_or_genes_scores_cistarget_dbs.py`](#combine_partial_motifs_or_tracks_vs_regions_or_genes_scores_cistarget_dbspy) | Combine partial cisTarget motifs or tracks vs regions or genes scores databases to: **1)** a complete cisTarget motifs or tracks regions or genes scores database and **2)** a complete cisTarget regions or genes vs motifs or tracks scores database.    |
|           [`combine_partial_regions_or_genes_vs_motifs_or_tracks_scores_cistarget_dbs.py`](#combine_partial_regions_or_genes_vs_motifs_or_tracks_scores_cistarget_dbspy) | Combine partial cisTarget regions or genes vs motifs or tracks scores databases to: **1)** a complete cisTarget regions or genes vs motifs or tracks scores database and **2)** a complete cisTarget motifs or tracks vs regions or genes scores database. |
| [`convert_motifs_or_tracks_vs_ regions_or_genes_scores_to_ rankings_cistarget_dbs.py`](#convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbspy) | Convert cisTarget motifs or tracks vs regions or genes scores database to cisTarget rankings database.                                                                                                                                                     |
|                                                                               [`create_cross_species_motifs_rankings_db.py`](#create_cross_species_motifs_rankings_dbpy) | Create cisTarget cross-species motifs rankings databases.                                                                                                                                                                                                  |


### Usage


#### create_cistarget_motif_databases.py

```bash
❯ ${create_cistarget_databases_dir}/create_cistarget_motif_databases.py --help
usage: create_cistarget_motif_databases.py [-h] -f FASTA_FILENAME [-F ORIGINAL_SPECIES_FASTA_FILENAME]
                                           -M MOTIFS_DIR -m MOTIFS_LIST_FILENAME
                                           [-5 MOTIF_MD5_TO_MOTIF_ID_FILENAME] -o DB_PREFIX
                                           [-c CLUSTER_BUSTER_PATH] [-t NBR_THREADS]
                                           [-p CURRENT_PART NBR_TOTAL_PARTS]
                                           [-g EXTRACT_GENE_ID_FROM_REGION_ID_REGEX_REPLACE]
                                           [-b BG_PADDING] [--min MIN_NBR_MOTIFS] [--max MAX_NBR_MOTIFS]
                                           [-l] [-s SEED] [-r SSH_COMMAND]

Create cisTarget motif databases.

options:
  -h, --help            show this help message and exit
  -f FASTA_FILENAME, --fasta FASTA_FILENAME
                        FASTA filename which contains the regions/genes to score with Cluster-Buster for
                        each motif. When creating a cisTarget species database from regions/genes lifted
                        over from a different species, provide the original FASTA file for that species
                        to -F.
  -F ORIGINAL_SPECIES_FASTA_FILENAME, --fasta-original-species ORIGINAL_SPECIES_FASTA_FILENAME
                        FASTA filename which contains all the regions/genes of the original species. The
                        fasta file provided to -f can contain less regions (not all regions could be
                        lifted over) than the one provided to -F, but to create a cisTarget cross-
                        species database later, all individual cisTarget species databases need to
                        contain the same amount of regions/genes.
  -M MOTIFS_DIR, --motifs_dir MOTIFS_DIR
                        Path to directory with Cluster-Buster motifs.
  -m MOTIFS_LIST_FILENAME, --motifs MOTIFS_LIST_FILENAME
                        Filename with list of motif IDs or motif MD5 names to be scored from directory
                        specified by "--motifs_dir".
  -5 MOTIF_MD5_TO_MOTIF_ID_FILENAME, --md5 MOTIF_MD5_TO_MOTIF_ID_FILENAME
                        Filename with motif MD5 to motif ID mappings to map Cluster-Buster motif MD5
                        filenames to motif IDs.
  -o DB_PREFIX, --output DB_PREFIX
                        Feather database prefix output filename.
  -c CLUSTER_BUSTER_PATH, --cbust CLUSTER_BUSTER_PATH
                        Path to Cluster-Buster (https://github.com/weng-lab/cluster-buster/). Default:
                        "cbust".
  -t NBR_THREADS, --threads NBR_THREADS
                        Number of threads to use when scoring motifs. Default: 1.
  -p CURRENT_PART NBR_TOTAL_PARTS, --partial CURRENT_PART NBR_TOTAL_PARTS
                        Divide the motif list in a number of total parts (of similar size) and score
                        only the part defined by current_part. This allows creating partial databases on
                        machines which do not have enough RAM to score all motifs in one iteration. This
                        will only create a partial regions/genes vs motifs scoring database ({db_prefix}
                        .part_000{current_part}_of_000{nbr_total_parts}.regions_vs_motifs.scores.feather
                        or {db_prefix}.part_000{current_part}_of_000{nbr_total_parts}.genes_vs_motifs.sc
                        ores.feather).
  -g EXTRACT_GENE_ID_FROM_REGION_ID_REGEX_REPLACE, --genes EXTRACT_GENE_ID_FROM_REGION_ID_REGEX_REPLACE
                        Take top CRM score for a gene by taking the maximum CRM score of multiple
                        regions for that gene. Define a regex which will remove the non-gene part of the
                        region ID, so only the gene ID remains. Examples: "gene_id#some_number":
                        "#[0-9]+$" or "region_id@@gene_id": "^.+@@".
  -b BG_PADDING, --bgpadding BG_PADDING
                        Background padding in bp that was added for each sequence in FASTA file.
                        Default: 0.
  --min MIN_NBR_MOTIFS  Minimum number of motifs needed per Cluster-Buster motif file to be considered
                        for scoring (filters motifs list). Default: 1.
  --max MAX_NBR_MOTIFS  Maximum number of motifs needed per Cluster-Buster motif file to be considered
                        for scoring (filters motifs list). Default: None.
  -l, --mask            Consider masked (lowercase) nucleotides as Ns.
  -s SEED, --seed SEED  Random seed used for breaking ties when creating rankings for a range of tied
                        scores. When setting this seed to a specific value and running this script with
                        the same input, will result in the same rankings databases as output.
  -r SSH_COMMAND, --ssh SSH_COMMAND
                        If defined, run Cluster-Buster over ssh by running the provided command to make
                        the connection before running Cluster-Buster itself. Example: 'ssh -o
                        ControlMaster=auto -o ControlPath=/tmp/ssh-control-path-%l-%h-%p-%r -o
                        ControlPersist=600 <hostname>'
```


#### create_cistarget_track_databases.py

```bash
❯ ${create_cistarget_databases_dir}/create_cistarget_track_databases.py --help
usage: create_cistarget_track_databases.py [-h] -b BED_FILENAME -T TRACKS_DIR -d TRACKS_LIST_FILENAME
                                           -o DB_PREFIX [-a BIGWIG_AVERAGE_OVER_BED_PATH]
                                           [-t NBR_THREADS] [-p CURRENT_PART NBR_TOTAL_PARTS]
                                           [-g EXTRACT_GENE_ID_FROM_REGION_ID_REGEX_REPLACE] [-s SEED]
                                           [-r SSH_COMMAND]

Create cisTarget track databases.

options:
  -h, --help            show this help message and exit
  -b BED_FILENAME, --bed BED_FILENAME
                        BED filename which contains the regions/genes to score with
                        bigWigAverageOverBed for each bigwig track (ChIP-seq) files.
  -T TRACKS_DIR, --tracks_dir TRACKS_DIR
                        Path to directory with bigwig track (ChIP-seq) files.
  -d TRACKS_LIST_FILENAME, --tracks TRACKS_LIST_FILENAME
                        Filename with list of track IDs to be scored from directory specified by "--
                        tracks_dir".
  -o DB_PREFIX, --output DB_PREFIX
                        Feather database prefix output filename.
  -a BIGWIG_AVERAGE_OVER_BED_PATH, --bwaob BIGWIG_AVERAGE_OVER_BED_PATH
                        Path to bigWigAverageOverBed
                        (http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed).
                        Default: "bigWigAverageOverBed".
  -t NBR_THREADS, --threads NBR_THREADS
                        Number of threads to use when scoring tracks. Default: 1.
  -p CURRENT_PART NBR_TOTAL_PARTS, --partial CURRENT_PART NBR_TOTAL_PARTS
                        Divide the tracks list in a number of total parts (of similar size) and score
                        only the part defined by current_part. This allows creating partial databases
                        on machines which do not have enough RAM to score all tracks in one iteration.
                        This will only create a partial regions/genes vs tracks scoring database ({db_p
                        refix}.part_000{current_part}_of_000{nbr_total_parts}.regions_vs_tracks.scores.
                        feather or {db_prefix}.part_000{current_part}_of_000{nbr_total_parts}.genes_vs_
                        tracks.scores.feather).
  -g EXTRACT_GENE_ID_FROM_REGION_ID_REGEX_REPLACE, --genes EXTRACT_GENE_ID_FROM_REGION_ID_REGEX_REPLACE
                        Take top score for a gene by taking the maximum score of multiple regions for
                        that gene. Define a regex which will remove the non-gene part of the region ID,
                        so only the gene ID remains. Examples: "gene_id#some_number": "#[0-9]+$" or
                        "region_id@@gene_id": "^.+@@".
  -s SEED, --seed SEED  Random seed used for breaking ties when creating rankings for a range of tied
                        scores. When setting this seed to a specific value and running this script with
                        the same input, will result in the same rankings databases as output.
  -r SSH_COMMAND, --ssh SSH_COMMAND
                        If defined, run bigWigAverageOverBed over ssh by running the provided command
                        to make the connection before running bigWigAverageOverBed itself. Example:
                        'ssh -o ControlMaster=auto -o ControlPath=/tmp/ssh-control-path-%l-%h-%p-%r -o
                        ControlPersist=600 <hostname>'
```


#### combine_partial_motifs_or_tracks_vs_regions_or_genes_scores_cistarget_dbs.py

```bash
❯ ${create_cistarget_databases_dir}/combine_partial_motifs_or_tracks_vs_regions_or_genes_scores_cistarget_dbs.py --help
usage: combine_partial_motifs_or_tracks_vs_regions_or_genes_scores_cistarget_dbs.py
       [-h] -i INPUT -o OUTPUT_DIR

Combine partial cisTarget motifs or tracks vs regions or genes scores databases to: 1) a complete
cisTarget motifs or tracks vs regions or genes scores database and2) a complete cisTarget regions or
genes vs motifs or tracks scores database.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input directory or database prefix with partial cisTarget motif or track vs
                        regions or genes scores database Feather files.
  -o OUTPUT_DIR, --output OUTPUT_DIR
                        Output directory to which the 1) complete cisTarget motifs or tracks vs
                        regions or genes scores database Feather files and 2) complete cisTarget
                        regions or genes vs motif or track scores database Feather files will be
                        written.
```

#### combine_partial_regions_or_genes_vs_motifs_or_tracks_scores_cistarget_dbs.py

```bash
❯ ${create_cistarget_databases_dir}/combine_partial_regions_or_genes_vs_motifs_or_tracks_scores_cistarget_dbs.py --help
usage: combine_partial_regions_or_genes_vs_motifs_or_tracks_scores_cistarget_dbs.py
       [-h] -i INPUT -o OUTPUT_DIR

Combine partial cisTarget regions or genes vs motifs or tracks scores databases to: 1) a complete
cisTarget regions or genes vs motifs or tracks scores database and 2) a complete cisTarget motifs or
tracks vs regions or genes scores database.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input directory or database prefix with partial cisTarget regions or genes vs
                        motif or track scores database Feather files.
  -o OUTPUT_DIR, --output OUTPUT_DIR
                        Output directory to which the 1) complete cisTarget regions or genes vs motif
                        or track scores database Feather files and 2) complete cisTarget motifs or
                        tracks vs regions or genes scores database Feather files will be written.
```


#### convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py

```bash
❯ ${create_cistarget_databases_dir}/convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py --help
usage: convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py
       [-h] -i CT_SCORES_DB_MOTIFS_OR_TRACKS_VS_REGIONS_OR_GENES_FILENAME [-s SEED]

Convert cisTarget motifs or tracks vs regions or genes scores database to cisTarget rankings
database.

optional arguments:
  -h, --help            show this help message and exit
  -i CT_SCORES_DB_MOTIFS_OR_TRACKS_VS_REGIONS_OR_GENES_FILENAME, --db CT_SCORES_DB_MOTIFS_OR_TRACKS_VS_REGIONS_OR_GENES_FILENAME
                        cisTarget motifs or tracks vs regions or genes scores database filename.
                        The cisTarget rankings database Feather file will be written to the same
                        directory.
  -s SEED, --seed SEED  Random seed used for breaking ties when creating rankings for a range of tied
                        scores. When setting this seed to a specific value and running this script with
                        the same input, will result in the same cisTarget rankings databases as output.

```


#### create_cross_species_motifs_rankings_db.py

```bash
❯ ${create_cistarget_databases_dir}/create_cross_species_motifs_rankings_db.py -h
usage: create_cross_species_motifs_rankings_db.py [-h] -i INPUT -o OUTPUT_DIR

Create cisTarget cross-species motifs rankings databases.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input directory or database prefix with cisTarget motifs vs regions or genes
                        rankings databases per species.
  -o OUTPUT_DIR, --output OUTPUT_DIR
                        Output directory to which the cisTarget cross-species motifs rankings database
                        files will be written.
```



# Creating cisTarget databases: details

To create cisTarget databases:
  - FASTA file with regulatory regions:
      - gene-based
      - region-based
  - motifs or TF ChIP-seq tracks:
      - motifs: in Cluster-Buster format
      - tracks: bigWig files of TF ChIP-seq data


## Creating cisTarget motif databases

cisTarget motif databases can be created in 2 ways:
  - [Score all motifs at once and create rankings](#score-all-motifs-at-once-and-create-rankings)
  - [Score motifs in different parts and generate rankings in a separate step](#score-motifs-in-different-parts-and-generate-rankings-in-a-separate-step)


### Score all motifs at once and create rankings

Create cisTarget motif databases:
  - [`create_cistarget_motif_databases.py`](#create_cistarget_motif_databasespy)
      - for each motif score all regulatory regions and create a cisTarget motifs vs regions/genes scores db:
          - `*.motifs_vs_regions.scores.feather`
          - `*.motifs_vs_genes.scores.feather`
      - transpose cisTarget motifs vs regions/genes scores db to cisTarget regions/genes vs motifs scores db:
          - `*.regions_vs_motifs.scores.feather`
          - `*.genes_vs_motifs.scores.feather`
      - creating a ranking for each regulatory region per motif based on the CRM score of the motif for that region and
        create a cisTarget motifs vs regions/genes rankings db:
          - `*.motifs_vs_regions.rankings.feather`
          - `*.motifs_vs_genes.rankings.feather`
      - transpose cisTarget motifs vs regions/genes rankings db to cisTarget regions/genes vs motifs rankings db:
          - `*.regions_vs_motifs.rankings.feather`
          - `*.genes_vs_motifs.rankings.feather`


```bash
# FASTA file with sequences per region IDs / gene IDs.
fasta_filename=
# Directory with motifs in Cluster-Buster format.
motifs_dir=
# File with motif IDs (base name of motif file in ${motifs_dir}).
motifs_list_filename=
# cisTarget motif database output prefix.
db_prefix=

nbr_threads=22


"${create_cistarget_databases_dir}/create_cistarget_motif_databases.py" \
    -f "${fasta_filename}" \
    -M "${motifs_dir}" \
    -m "${motifs_list_filename}" \
    -o "${db_prefix}" \
    -t "${nbr_threads}"
```


### Score all tracks at once and create rankings

Create cisTarget tracks databases:
  - [`create_cistarget_track_databases.py`](#create_cistarget_track_databasespy)
      - for each track score all regulatory regions and create a cisTarget tracks vs regions/genes scores db:
          - `*.tracks_vs_regions.scores.feather`
          - `*.tracks_vs_genes.scores.feather`
      - transpose cisTarget track vs regions/genes scores db to cisTarget regions/genes vs tracks scores db:
          - `*.regions_vs_tracks.scores.feather`
          - `*.genes_vs_tracks.scores.feather`
      - creating a ranking for each regulatory region per motif based on the track score of the track for that region
        and create a cisTarget track vs regions/genes rankings db:
          - `*.tracks_vs_regions.rankings.feather`
          - `*.tracks_vs_genes.rankings.feather`
      - transpose cisTarget tracks vs regions/genes rankings db to cisTarget regions/genes vs tracks rankings db:
          - `*.regions_vs_tracks.rankings.feather`
          - `*.genes_vs_tracks.rankings.feather`


```bash
# BED file with regions to score.
regions_bed_filename=
# Directory with bigWig tracks of TF-ChIP-seq data.
tracks_dir=
# File with track IDs (base names of bigWig files in ${tracks_dir}).
tracks_list_filename=
# cisTarget track database output prefix.
db_prefix=

nbr_threads=22


"${create_cistarget_databases_dir}/create_cistarget_track_databases.py" \
    -b "${regions_bed_filename}" \
    -T "${tracks_dir}" \
    -d "${tracks_list_filename}" \
    -o "${db_prefix}" \
    -t "${nbr_threads}"
```



### Score motifs in different parts and generate rankings in a separate step

Create cisTarget motif databases:
  1) [`create_cistarget_motif_databases.py`](#create_cistarget_motif_databasespy):
      - score the whole list of motifs in several parts by running `create_cistarget_motif_databases.py` with the
        `-p ${current_part} ${nbr_total_parts}`  option with `${current_part}` set from `1` to `${nbr_total_parts}`.
        Each run will create motif scores (for the current subset of motifs) for all regulatory regions and create a partial
        cisTarget motifs vs regions/genes scores db:
          - `*.part_${current_part}_of_${nbr_total_parts}.motifs_vs_regions.scores.feather`
          - `*.part_${current_part}_of_${nbr_total_parts}.motifs_vs_genes.scores.feather`
  2) [`combine_partial_motifs_or_tracks_vs_regions_or_genes_scores_cistarget_dbs.py`](#combine_partial_motifs_or_tracks_vs_regions_or_genes_scores_cistarget_dbspy):
      - Combine partial cisTarget motifs vs regions/genes scores db to:
          - a complete cisTarget motifs vs regions/genes scores database:
              - `*.motifs_vs_regions.scores.feather`
              - `*.motifs_vs_genes.scores.feather`
          - a complete cisTarget regions/genes vs motifs scores database:
              - `*.regions_vs_motifs.scores.feather`
              - `*.genes_vs_motifs.scores.feather`
          - partial cisTarget scores databases can be deleted afterwards.
  3) [`convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py`](#convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbspy)
      - creating a ranking for each regulatory region per motif based on the CRM score of the motif for that region and
        create a cisTarget motifs vs regions/genes rankings db:
          - `*.motifs_vs_regions.rankings.feather`
          - `*.motifs_vs_genes.rankings.feather`
      - transpose cisTarget motifs vs regions/genes rankings db to cisTarget regions/genes vs motifs rankings db:
          - `*.regions_vs_motifs.rankings.feather`
          - `*.genes_vs_motifs.rankings.feather`


#### Step 1

Using `-p` or `--partial` of [`create_cistarget_motif_databases.py` (help)](#create_cistarget_motif_databasespy) will
divide the motif list in a number of total parts (`${nbr_total_parts}`) (of similar size) and score only the part
defined by`${current_part}`.

This allows creating partial databases on machines which do not have enough RAM to score all motifs in one iteration
and/or running the motif scoring on multiple nodes (where each node runs the motif scoring with a different value for
`${current_part}`) in parallel. This is quite useful if the number of region/genes is quite high.

This will only create a partial cisTarget motifs vs regions/genes scores database files:
   - `${output_db_prefix}.part_0*${current_part}_of_0*${nbr_total_parts}.motifs_vs_regions.scores.feather`
   - `${output_db_prefix}.part_0*${current_part}_of_0*${nbr_total_parts}.motifs_vs_genes.scores.feather`


```bash
# FASTA file with sequences per region IDs / gene IDs.
fasta_filename=
# Directory with motifs in Cluster-Buster format.
motifs_dir=
# File with motif IDs (base name of motif file in ${motifs_dir}).
motifs_list_filename=
# cisTarget motif database output prefix.
db_prefix=

nbr_threads=22
nbr_parts=10

# Create a partial directory, so partial cisTarget database files can be deleted easily afterwards.
mkdir partial

# Each invocation of the for loop (with different ${current_part}) can also be submitted to a different node to speedup
# the motif scoring.

for current_part in $(seq 1 ${nbr_total_parts}) ; do
    "${create_cistarget_databases_dir}/create_cistarget_motif_databases.py" \
         -f "${fasta_filename}" \
         -M "${motifs_dir}" \
         -m "${motifs_list_filename}" \
         -p "${current_part}" "${nbr_total_parts}" \
         -o "partial/${db_prefix}" \
         -t "${nbr_threads}"
done
```


When creating cross-species databases, motif scoring should be done with regulatory regions after liftover
(`fasta_filename`), while the original regulatory regions fasta file (`original_species_fasta_filename`) should also be
provided. The latter is only used to make sure that all regions/genes of the original regulatory regions are in the
generated cisTarget database as some regions might get lost after liftover. cisTarget motif score databases need to be
generated for each species of interest.

```bash
# FASTA file with sequences per region IDs / gene IDs.
fasta_filename=
# FASTA file with sequences per region IDs / gene IDs of the original species.
original_species_fasta_filename=
# Directory with motifs in Cluster-Buster format.
motifs_dir=
# File with motif IDs (base name of motif file in ${motifs_dir}).
motifs_list_filename=
# cisTarget motif database output prefix.
db_prefix=

nbr_threads=22


"${create_cistarget_databases_dir}/create_cistarget_motif_databases.py" \
    -f "${fasta_filename}" \
    -F "${original_species_fasta_filename}" \
    -M "${motifs_dir}" \
    -m "${motifs_list_filename}" \
    -o "${db_prefix}" \
    -t "${nbr_threads}"
```

#### Step 2

> *See [Memory requirements](#memory-requirements) to
have a rough guess about the amount of memory needed in case you have problems running this step.*


When all partial cisTarget regions/genes vs motifs scores database files are created:
  - `${db_prefix}.part_0*of_0*${nbr_total_parts}.motifs_vs_regions.scores.feather`
  - `${db_prefix}.part_0*_of_0*${nbr_total_parts}.motifs_vs_genes.scores.feather`

they can be combined with
[`combine_partial_motifs_or_tracks_vs_regions_or_genes_scores_cistarget_dbs.py` (help)](#combine_partial_rmotifs_or_tracks_vs_regions_or_genes_scores_cistarget_dbspy)
to:
  - a complete cisTarget regions/genes vs motifs scores database:
      - `*.motifs_vs_regions.scores.feather`
      - `*.motifs_vs_genes.scores.feather`
  - a complete cisTarget motifs vs regions/genes scores database:
      - `*.regions_vs_motifs.scores.feather`
      - `*.genes_vs_motifs.scores.feather`

Partial cisTarget scores databases can be deleted afterwards.


```bash
"${create_cistarget_databases_dir}/combine_partial_motifs_or_tracks_vs_regions_or_genes_cistarget_dbs.py \
    -i partial/ \
    -o .

# Partial cisTarget databases can be removed.
#rm -r partial
```


#### Step 3

> *See [Memory requirements](#memory-requirements) to
have a rough guess about the amount of memory needed in case you have problems running this step.*


Create rankings from a complete cisTarget regions/genes vs motifs scores database with
[`convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py` (help)](#convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbspy):
  - creating a ranking for each regulatory region per motif based on the CRM score of the motif for that region and
    create a cisTarget motifs vs regions/genes rankings db:
      - `*.motifs_vs_regions.rankings.feather`
      - `*.motifs_vs_genes.rankings.feather`
  - transpose cisTarget motifs vs regions/genes rankings db to cisTarget regions/genes vs motifs rankings db:
      - `*.regions_vs_motifs.rankings.feather`
      - `*.genes_vs_motifs.rankings.feather`


```bash
# cisTarget database filename:
#   - *.motifs_vs_regions.rankings.feather
#   - *.motifs_vs_genes.rankings.feather
db_filename=

"${create_cistarget_databases_dir}/convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py \
    -i "${db_filename}"
```



## Create cisTarget cross-species motifs rankings database

To create cisTarget cross-species motifs rankings database:
  - liftover regulatory regions of main species to each species of interest.
  - [`Create cisTarget motif databases`](#creating-cistarget-motif-databases) for each lifted over regulatory regions
    FASTA file.convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py
  - Create cisTarget cross-species motifs rankings from individual (per species) cisTarget motifs rankings databases
    with [`create_cross_species_motifs_rankings_db.py`](#create_cross_species_motifs_rankings_dbpy)

```bash
# cisTarget database prefix which matches the common part of all cisTarget rankings databases (or just the directory
# that contains them).
db_prefix=

# Output directory.
output_dir=


"${create_cistarget_databases_dir}/create_cross_species_motifs_rankings_db.py \
    -i "${db_prefix}" \
    -o "${output_dir}"
```


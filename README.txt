
## Installation and Usage

### Step 1: Create a Conda virtual environment

Create a Conda virtual environment named `OHRC-Sentinel` with Python 3.10:

conda create --name OHRC-Sentinel python=3.10
conda activate OHRC-Sentinel

Step 2: Configure Conda channels
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

Step 3: Install required dependencies

conda install abricate=1.0.1
conda install padloc=2.0.0
conda install snippy
conda install mob_suite
conda install r-fastbaps
pip install pandas numpy openpyxl scikit-learn joblib biopython


Step 4: Run OHRC-Sentinel

python create_summary_tab.py "./*.fna" ./output/ metadata.tsv

PS: The `metadata.tsv` file should be tab-separated and formatted as follows, with `Sample_ID`, `ST`, and `Source` as the required column names:

Sample_ID    ST      Source
sample_001   ST11    Human
sample_002   ST11    Animal
sample_003   ST11    Food
sample_004   ST11    Environment
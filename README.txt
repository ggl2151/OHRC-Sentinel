Step 1: Create a conda virtual environment named “PPRP”
conda create --name PPRP python=3.10
conda activate PPRP

Step 2: Install required dependenciesconda install abricate
conda install abricate  (version 1.0.1)
conda install padloc   (version 2.0.0)


Step 3: python create_summary_tab.py "./*.fna" ./output/
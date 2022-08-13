# Variant Calling

## Introduction:
This tool was designed to call SNPs from Nanopore data using a reference genome and produce a phylogenetic tree for epidemiology purposes. The Variant calling is done using the [PEPPER-Margin-DeepVariant](https://github.com/kishwarshafin/pepper) pipeline. RAXmL is used to produced the final tree.


## Installation
### Docker
The PEPPER-Margin-DeepVariant pipeline runs through docker. Please refer to official documentation to install docker:
* https://docs.docker.com/engine/install/ubuntu/

You also need to run the post-installation steps to be able to run docker without sudo:
* https://docs.docker.com/engine/install/linux-postinstall/

The pipeline has been designed and tested with `docker.io=20.10.7-0ubuntu5~20.04.2`.

### Virtual environment
It is highly recommended installing this pipeline in a virtual environment to manage all the dependencies and avoid conflicts with other tools on your system. I recommend using [miniconda](https://docs.conda.io/en/latest/miniconda.html#linux-installers) to manage the virtual environment and installing dependencies. I also find `mamba` to be more efficient than `conda` to install packages. See [here](https://github.com/mamba-org/mamba) for more info.
```
# Create conda environment (replace "conda" by "mamba" if you installed it)
conda create -n pepper -c bioconda -y bcftools=1.9 samtools=1.9 \
    minimap2=2.24  minimap2=2.24 vcftools=0.1.16 raxml=8.2.12 git

# Activate environment (don't use "mamba" instead of "conda" here)
conda activate pepper

# Clone repository
git clone https://github.com/duceppemo/nanopore_variant_calling_pepper

# Test pipeline. The help file should print on screen without any error messages.
cd nanopore_variant_calling_pepper
python Pepper_ont.py -h 
```
## Usage
```
usage: python pepper_ont.py [-h] -i /path_to_input_fastq/ -r reference.fasta -o /path_to_output_folder/ [-t 4] [-p PARALLEL]

VCF and GVCF file generation

optional arguments:
  -h, --help            show this help message and exit
  -i /path_to_input_fastq/, --input /path_to_input_fastq/
                        Input folder containing nanopore fastq files.
  -r reference.fasta, --reference reference.fasta
                        Reference fasta file.
  -o /path_to_output_folder/, --output /path_to_output_folder/
                        Output folder.
  -t 4, --threads 4     Number of threads.
  -p PARALLEL, --parallel PARALLEL
                        Number of samples to run in parallel.
```

# CFIA
# Instructions for how to install all package dependencies and execute the script

# Instructions on how to install Docker with Image
$ sudo apt install docker.io=20.10.7-0ubuntu5~20.04.2 \
$ sudo docker run hello-world will download a test container to confirm that the docker has been \
successfully downloaded \
The command to download the image is available in the run_pepper_iter function. Should only be used if the \
image container has not already been downloaded. 

# Install bcftools, samtools, gatk and minimap2
conda install -c bioconda bcftools=1.9 \
conda install -c bioconda samtools=1.9 \
conda install -c bioconda gatk4=4.2.5.0 \
conda install -c bioconda minimap2=2.24 

# Instructions to run docker without the sudo command
$ sudo groupadd docker \
$ sudo gpasswd -a $USER docker \
$ newgrp docker \
$ docker run hello-world to check if docker can run without sudo 

# How to run this script in terminal:
Choose directory in which script is located and then run the following command in terminal: \
Create a conda environment for Pepper (run in terminal: conda create --name Pepper) \
Conda activate Pepper \
Run sudo setfacl -m user:($USER):rw /var/run/docker.sock to be able to execute docker and then run the \
following command \
python Pepper_ont.py -i [folder with the query fastq files] -r [reference fasta file] -o [output folder]  [-t 4] 
[-p PARALLEL]

-i folder containing the query fastq files \
-r reference fasta file  \
-o output folder \
-t number of threads \
-p number of samples to be run in parallel

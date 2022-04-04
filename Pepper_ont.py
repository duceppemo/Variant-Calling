import subprocess
import os
import argparse
from concurrent import futures
from pathlib import Path
import logging


# Instructions on how to run Docker with Image
# $ sudo apt-get install docker-ce=<VERSION_STRING> docker-ce-cli=<VERSION_STRING> containerd.io
# (version e.g 5:18.09.1~3-0~ubuntu-xenial)
# $ sudo docker run hello-world will download a test container to confirm that the docker has been
# successfully downloaded
# The command to download the image is available in the run_pepper_iter function. Should only be used if the
# image container has not already been downloaded.

# Instructions to run docker without the sudo command
# $ sudo groupadd docker
# $ sudo gpasswd -a $USER docker
# $ newgrp docker
# $ docker run hello-world to check if docker can run without sudo

# How to run this script in Pycharm:
# Choose directory in which script is located and then run the following command in terminal:
# Create a conda environment for Pepper (run in terminal: conda create --name Pepper)
# Conda activate Pepper
# Run sudo setfacl -m user:($USER):rw /var/run/docker.sock to be able to execute docker and then run the
# following command
# python Pepper_ont.py -i [folder with the query fastq files] -r [reference fasta file] -o [output folder]  [-t 4]
# [-p PARALLEL]


class Pepper:
    def __init__(self, args):
        # Logging information of all the tools used in thie script
        logging.basicConfig(level=logging.INFO)
        info = 'bcftools v=1.9, samtools v=1.9, minimap2 v=2.24-r1122, gatk v=4.2.5.0, docker v=20.10.7'
        logging.info(info)

        # Arguments from the command line
        self.fastq_folder = args.input
        self.reference = args.reference
        self.output_folder = args.output
        self.cpu = args.threads
        self.parallel = args.parallel

        # Output folders needed for various steps
        self.bam_folder = self.output_folder + '/bam/'
        self.pepper_folder = self.output_folder + '/pepper/'

        # List all the input fastq files
        self.fastq_list = Methods_Calling.list_file(self.fastq_folder, (".fastq", ".fastq.gz", ".fq", ".fq.gz"))

        # Run
        Pepper.run_all(self)

    def run_all(self):

        # Test docker installation
        # Only use this command to verify whether docker is working after which it can be removed
        Methods_Calling.docker_test()

        # Creates the output folder
        Methods_Calling.folder_create(self.output_folder)
        Methods_Calling.folder_create(self.bam_folder)
        Methods_Calling.folder_create(self.pepper_folder)

        # Creates the bam file for each query file in the list
        Methods_Calling.map_fastq_parallel(self.reference, self.fastq_list, self.bam_folder, self.cpu, self.parallel)

        # Produces the GVCF and VCF files using variant calling
        Methods_Calling.run_pepper_iter(self.reference, self.bam_folder, self.pepper_folder, self.cpu)

        # Merges all GVCF files into a single VCF file
        g_vcf_list = Methods_Calling.list_file(self.pepper_folder, ".g.vcf.gz")
        Methods_Calling.gatk_merge_gvcf(self.reference, g_vcf_list, self.pepper_folder)


class Methods_Calling:
    @staticmethod
    def docker_test():
        docker_test_cmd = ['docker', 'run', 'hello-world']
        if subprocess.run(docker_test_cmd):
            pass

    @staticmethod
    def list_file(folder, extension):
        file_list = list()
        for root, dirs, files in os.walk(folder):
            for file in files:
                if file.endswith(extension):
                    file_list.append(os.path.join(root, file))
        return file_list

    @staticmethod
    def folder_create(folder):
        # Create output folder
        Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def map_fastq(ref_file, fastq_file, output_folder, cpu):
        sample_name = (os.path.basename(fastq_file).split(".")[0]).split("_")[0]
        bam_out = output_folder + sample_name + '.bam'
        # Creating the bam file using minimap2 and samtools
        minimap2_cmd = ['minimap2',
                        '-x', 'map-ont',
                        '-a',
                        '-t', str(cpu),
                        ref_file, fastq_file]

        samtools_cmd = ['samtools', 'sort',
                        '-o', bam_out, '-@', str(cpu)]

        # Pipe processes to avaoid wrtiting intermediate SAM file to disk
        p1 = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(samtools_cmd, stdin=p1.stdout)
        p1.stdout.close()
        p2.communicate()

        # Index the BAM file
        index_cmd = ['samtools', 'index', bam_out]
        subprocess.run(index_cmd)

    @staticmethod
    def map_fastq_parallel(ref_file, fastq_list, bam_folder, cpu, parallel):
        with futures.ThreadPoolExecutor(max_workers=parallel) as executor:
            args = ((ref_file, fastq_file, bam_folder, int(cpu / parallel))
                    for fastq_file in fastq_list)
            for results in executor.map(lambda x: Methods_Calling.map_fastq(*x), args):
                pass

    @staticmethod
    def run_pepper(ref_file, bam_file, pepper_folder, cpu):
        ref_folder = os.path.dirname(ref_file)
        ref_name = os.path.basename(ref_file)
        bam_folder = os.path.dirname(bam_file)
        bam_name = os.path.basename(bam_file)
        sample_name = (bam_name.split(".")[0]).split("_")[0]
        # Generate the gvcf and vcf files with pepper + deepvariant
        pepper_cmd = ['docker', 'run',
                      '-v', '{}:/ref'.format(ref_folder),
                      '-v', '{}:/bam'.format(bam_folder),
                      '-v', '{}:/output'.format(pepper_folder),
                      'kishwars/pepper_deepvariant:r0.7',
                      'run_pepper_margin_deepvariant', 'call_variant',
                      '-b', "/bam/{}".format(bam_name),
                      '-f', "/ref/{}".format(ref_name),
                      '-o', "/output/",
                      '-p', sample_name,
                      '-t', str(cpu),
                      '--ont_r9_guppy5_sup',
                      '--gvcf']
        subprocess.run(pepper_cmd)

    @staticmethod
    def run_pepper_iter(ref_file, bam_folder, pepper_folder, cpu):
        # Pull docker if image not already downloaded
        docker_pull_cmd = ['docker', 'pull', 'kishwars/pepper_deepvariant:r0.7']
        subprocess.run(docker_pull_cmd)
        bam_list = Methods_Calling.list_file(bam_folder, '.bam')
        print(bam_list)

        for bam_file in bam_list:
            Methods_Calling.run_pepper(ref_file, bam_file, pepper_folder, cpu)

    @staticmethod
    def gatk_merge_gvcf(ref_file, g_vcf_list, pepper_folder):
        merged_gvcf_file = pepper_folder + "mycobacterium_sp.g.vcf"
        merged_vcf_file = pepper_folder + 'mycobacterium_sp.vcf'
        gvcf_create_sequence_dict_cmd = ['gatk', 'CreateSequenceDictionary', '-R', ref_file]

        vcf_create_cmd = ['bcftools', 'convert',
                          '--gvcf2vcf', merged_gvcf_file, '-f', ref_file, '-o', merged_vcf_file]

        gvcf_cmd = ["gatk", 'CombineGVCFs']
        for file in g_vcf_list:
            gvcf_cmd = gvcf_cmd + ["-V", file]
        gvcf_cmd = gvcf_cmd + ['-R', ref_file, '-O', merged_gvcf_file]

        subprocess.run(gvcf_create_sequence_dict_cmd)
        subprocess.run(gvcf_cmd)
        subprocess.run(vcf_create_cmd)


if __name__ == "__main__":
    # Required arguments
    parser = argparse.ArgumentParser(description="VCF and GVCF file generation")
    parser.add_argument('-i', '--input', metavar='/path_to_input_fastq/',
                        required=True, type=str,
                        help="Input folder containing nanopore fastq files.")
    parser.add_argument('-r', "--reference", metavar='reference.fasta',
                        required=True, type=str,
                        help="Reference fasta file.")
    parser.add_argument('-o', "--output", metavar='/path_to_output_folder/',
                        required=True, type=str,
                        help="Output folder.")
    parser.add_argument('-t', "--threads", metavar='4',
                        required=False, type=int, default=4,
                        help="Number of threads.")
    parser.add_argument('-p', "--parallel",
                        required=False, type=int, default=1,
                        help="Number of samples to run in parallel.")
    arg = parser.parse_args()
    Pepper(arg)

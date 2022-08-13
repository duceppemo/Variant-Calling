import subprocess
import os
import argparse
from concurrent import futures
from pathlib import Path
import logging
import pandas as pd


# Instructions on how to install Docker with Image
# $ sudo apt install docker.io=20.10.7-0ubuntu5~20.04.2
# $ sudo docker run hello-world will download a test container to confirm that the docker has been
# successfully downloaded
# The command to download the image is available in the run_pepper_iter function. Should only be used if the
# image container has not already been downloaded.

# Instructions to run docker without the sudo command
# $ sudo groupadd docker
# $ sudo gpasswd -a $USER docker
# $ newgrp docker
# $ docker run hello-world to check if docker can run without sudo

# How to run this script in terminal:
# Choose directory in which script is located and then run the following command in terminal:
# Create a conda environment for Pepper (run in terminal: conda create --name Pepper)
# Conda activate Pepper
# Run sudo setfacl -m user:($USER):rw /var/run/docker.sock to be able to execute docker and then run the
# following command
# python pepper_ont.py -i [folder with the query fastq files] -r [reference fasta file] -o [output folder]  [-t 4]
# [-p PARALLEL]

class Pepper:
    def __init__(self, args):
        # Logging information for all the required tools needs to run the script
        logging.basicConfig(level=logging.INFO)
        info = 'Required tools: bcftools==1.9, samtools==1.9, minimap2==2.24, gatk==4.2.5.0, docker==20.10.7, ' \
               'vcftools==0.1.16'
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
        self.fastq_list = MethodsCalling.list_file(self.fastq_folder, (".fastq", ".fastq.gz", ".fq", ".fq.gz"))

        # Run
        Pepper.run_all(self)

    def run_all(self):
        # Checking version number for gatk, minimap2, samtools, bcftools and docker
        MethodsCalling.check_version()

        # Creates the output folder
        MethodsCalling.create_folder(self.output_folder)
        MethodsCalling.create_folder(self.bam_folder)
        MethodsCalling.create_folder(self.pepper_folder)

        # Creates the bam file with read group tag for each query file in the list
        MethodsCalling.map_fastq_parallel(self.reference, self.fastq_list, self.bam_folder, self.cpu, self.parallel)

        # Produces the gVCF and VCF files using variant calling and renames samples with its accession number
        MethodsCalling.run_pepper_iter(self.reference, self.bam_folder, self.pepper_folder, self.cpu)

        # create a list of all the gVCF files
        g_vcf_list = MethodsCalling.list_file(self.pepper_folder, ".g.vcf.gz")

        # Merges all the new gVCF files to create a multi-sample gVCF and VCF
        MethodsCalling.merge_gvcf(self.reference, g_vcf_list, self.cpu, self.output_folder)

        # Reformat the multi-sample VCF file
        vcf_in = self.output_folder + '/multi_sample.vcf'
        vcf_out = self.output_folder + '/multi_sample_v1.vcf'
        MethodsCalling.fix_merged_vcf(vcf_in, vcf_out)

        # Filters through the new multi-sample vcf and selects sites that are bi-allelic, contain the
        # desired quality score etc. All sites that do not pass these requirements are removed
        # Adjust the filter settings accordingly
        vcf_filtered = self.output_folder + '/multi_sample_final'
        MethodsCalling.filter_vcf(vcf_out, vcf_filtered)

        # Converts the multi-sample VCF into an aligned fasta file
        final_vcf = self.output_folder + '/multi_sample_final.recode.vcf'
        fasta_file = self.output_folder + '/final.fasta'
        MethodsCalling.vcf2fasta(final_vcf, fasta_file)

        # Creates a phylogenetic tree from the aligned fasta file
        MethodsCalling.make_tree_raxml(fasta_file, self.output_folder, self.cpu)


class MethodsCalling:
    @staticmethod
    def check_version():
        # Checks the version number for all the tools to see if they have been installed
        docker_vers_check_cmd = ['docker', '--version']
        bcftools_vers_check_cmd = ['bcftools', '--version']
        samtools_vers_check_cmd = ['samtools', '--version']
        minimap2_vers_check_cmd = ['minimap2', '--version']
        vcftools_vers_check_cmd = ['vcftools', '-h']
        gatk_vers_check_cmd = ['gatk', '--version']

        # TODO: stop running if dependency is missing
        # TODO: log versions in log file for QA
        subprocess.run(docker_vers_check_cmd)
        subprocess.run(bcftools_vers_check_cmd)
        subprocess.run(samtools_vers_check_cmd)
        subprocess.run(minimap2_vers_check_cmd)
        subprocess.run(gatk_vers_check_cmd)
        subprocess.run(vcftools_vers_check_cmd)

    @staticmethod
    def list_file(folder, extension):
        # Creates a list of all the query files
        file_list = list()
        for root, dirs, files in os.walk(folder):
            for file in files:
                if file.endswith(extension):
                    file_list.append(os.path.join(root, file))

        # Returns this query folder is empty if the there are no query files present
        if file_list == list():
            print("{} is empty".format(folder))
            exit()
        else:
            return file_list

    @staticmethod
    def create_folder(folder):
        # Create output folder. Will create parent folder is not exist. Won't warn is folder already exists.
        Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def map_fastq(ref_file, fastq_file, output_folder, cpu):
        # The sample name/accession number of the query file
        sample_name = (os.path.basename(fastq_file).split(".")[0]).split("_")[0]

        # The bam file name for the respected query file
        bam_out = output_folder + sample_name + '.bam'

        # Stops the program if ref file is missing
        if os.listdir(os.path.dirname(ref_file)) == list():
            print('ref file not found')

        # Creating the bam file using minimap2 and samtools
        minimap2_cmd = ['minimap2',
                        '-x', 'map-ont',
                        '-a',
                        '-R', "@RG\\tID:{}\\tSM:{}".format(sample_name, sample_name),
                        '-t', str(cpu),
                        ref_file, fastq_file]

        samtools_cmd = ['samtools', 'sort',
                        '-@', str(cpu),
                        '-o', bam_out]

        # Pipe processes to avoid writing intermediate SAM file to disk
        p1 = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(samtools_cmd, stdin=p1.stdout)
        p1.stdout.close()
        p2.communicate()

        # Index the created bam file
        index_cmd = ['samtools', 'index', bam_out]
        subprocess.run(index_cmd)

    @staticmethod
    def map_fastq_parallel(ref_file, fastq_list, bam_folder, cpu, parallel):
        # Maps the fastq files to the reference in parallel
        with futures.ThreadPoolExecutor(max_workers=parallel) as executor:
            args = ((ref_file, fastq_file, bam_folder, int(cpu / parallel))
                    for fastq_file in fastq_list)
            for results in executor.map(lambda x: MethodsCalling.map_fastq(*x), args):
                pass

    @staticmethod
    def run_pepper(ref_file, bam_file, pepper_folder, cpu):
        # The directory of the reference file
        ref_folder = os.path.dirname(ref_file)

        # Name of the reference file
        ref_name = os.path.basename(ref_file)

        # Directory of the bam file
        bam_folder = os.path.dirname(bam_file)

        # Name of the bam file
        bam_name = os.path.basename(bam_file)

        # Unique sample name/accession number
        sample_name = bam_name.split(".")[0].split("_")[0]

        # Generate the gVCF and vcf files with pepper_deepvariant
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
                      '-s', sample_name,
                      '--ont_r9_guppy5_sup',
                      '--gvcf']
        subprocess.run(pepper_cmd)

    @staticmethod
    def run_pepper_iter(ref_file, bam_folder, pepper_folder, cpu):
        # Pull docker if image not already downloaded
        # docker_pull_cmd = ['docker', 'pull', 'kishwars/pepper_deepvariant:r0.7']
        # subprocess.run(docker_pull_cmd)
        bam_list = MethodsCalling.list_file(bam_folder, '.bam')

        # Stops the code if the bam files have not been generated
        if bam_list == list():
            print("bam folder is empty")
            exit()
        # Uses the run_pepper function to create gcvf and vcf files for each bam file
        else:
            for bam_file in bam_list:
                MethodsCalling.run_pepper(ref_file, bam_file, pepper_folder, cpu)

    @staticmethod
    def merge_gvcf(ref_file, g_v_list, cpu, output):
        # Merges the gVCF files into a single multi-sample gVCF and VCF
        merged_gvcf_file = output + '/multi_sample.g.vcf'
        merged_vcf_file = output + '/multi_sample.vcf'
        gvcf_create_cmd = ['bcftools', 'merge',
                           '-g', ref_file,
                           '-O', "v",
                           '-o', merged_gvcf_file,
                           '--threads', str(cpu),
                           ] + g_v_list
        vcf_create_cmd = ['bcftools',
                          'convert',
                          '--gvcf2vcf',
                          '-o', merged_vcf_file,
                          '-O', 'v',
                          '--threads', str(cpu),
                          '-f', ref_file,
                          merged_gvcf_file]
        subprocess.run(gvcf_create_cmd)
        subprocess.run(vcf_create_cmd)

    @staticmethod
    def fix_merged_vcf(input_vcf, output_vcf):
        # Creates a new file with the positions that contain variants and removes the lines that have <*> as
        # these positions contain no variants throughout all species. Creates version 1 of the modified file
        with open(output_vcf, 'w') as out_file:
            with open(input_vcf, 'r') as in_file:
                for line in in_file:
                    if '\t<*>\t' in line:
                        continue
                    else:
                        line = line.replace('<*>,', '')
                        line = line.replace(',<*>', '')
                        line = line.replace('2/2', '1/1')
                        line = line.replace('0/2', '0/1')
                    out_file.write(line)

    @staticmethod
    def filter_vcf(vcf_in, vcf_out):
        # Creates the final version of the multi-sample vcf by filtering through the file
        filter_cmd = ['vcftools',
                      '--vcf', vcf_in,
                      '--out', vcf_out,
                      '--remove-indels',
                      # '--maf', str(0.05),
                      '--min-alleles', str(2),
                      '--max-alleles', str(2),
                      '--minDP', str(10),
                      # '--max-missing', str(0.2),
                      '--minQ', str(10),
                      '--remove-filtered-all',
                      '--recode']
        subprocess.run(filter_cmd)

    @staticmethod
    def vcf2fasta(input_vcf, output_fasta):
        iupac_dict = {'AG': 'R',
                      'CT': 'Y',
                      'GC': 'S',
                      'AT': 'W',
                      'GT': 'K',
                      'AC': 'M'}

        with open(output_fasta, 'w') as out_fh:
            with open(input_vcf, 'r') as in_fh:
                for line in in_fh:
                    if line.startswith('##'):
                        # out_fh.write(line)
                        continue
                    elif line.startswith('#CHROM'):
                        sample_list = line.rstrip().split('\t', 9)[9].split('\t')
                        df = pd.DataFrame(columns=sample_list)
                    else:
                        # Split data lines into 10 fields, the last one is the samples info
                        field_list = line.rstrip().split('\t', 9)
                        # Dictionary to convert 0/0 and 1/1 geno to REF or ALT call
                        ref = field_list[3]
                        alt = field_list[4]
                        try:
                            iupac = iupac_dict[ref + alt]
                        except KeyError:
                            iupac = iupac_dict[alt + ref]

                        geno_dict = {'0/0': ref,
                                     '1/1': alt,
                                     '0/1': iupac,
                                     '1/0': iupac,
                                     './.': '-'}
                        # Split the last field (sample info) and only keep the genotype
                        df.loc[len(df)] = [geno_dict[x.split(':')[0]] for x in field_list[9].split('\t')]

                # Convert dataframe to fasta output
                fasta_dict = df.to_dict(orient='list')
                for sample, seq_list in fasta_dict.items():
                    out_fh.write('>{}\n{}\n'.format(sample, ''.join(seq_list)))

    @staticmethod
    def make_tree_raxml(aligned_fasta, output_folder, cpu):

        cmd = ['raxmlHPC-PTHREADS-AVX2',
               '-s', aligned_fasta,
               '-w', output_folder,
               '-n', '{}.tree'.format('.'.join(os.path.basename(aligned_fasta).split('.')[:-1])),
               '-m', 'GTRCAT',
               '-N', str(100),
               '-d', '-f', 'a', '-T', str(cpu),
               '-x', str(1234), '-p', str(123)]

        subprocess.run(cmd)

    @staticmethod
    def run_raxml(fasta):
        # Generates the Phylogenetic tree
        phylo_tree_cmd = ['raxmlHPC-PTHREADS-AVX2',
                          '-s', fasta,
                          '-n', 'multi-sample_tree',
                          '-m', 'GTRCAT',
                          '-f', 'a',
                          '-x', '123',
                          '-N', 'autoMRE',
                          '-p', '456']
        return subprocess.run(phylo_tree_cmd)


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

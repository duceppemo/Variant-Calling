import subprocess
import os
import argparse
from concurrent import futures
from pathlib import Path
import logging
import pandas as pd
import io


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
# python Pepper_ont.py -i [folder with the query fastq files] -r [reference fasta file] -o [output folder]  [-t 4]
# [-p PARALLEL]

```python
class Pepper:
    def __init__(self, args):
        # Logging information for all the required tools needs to run the script
        logging.basicConfig(level=logging.INFO)
        info = 'Required tools: bcftools v=1.9, samtools v=1.9, minimap2 v=2.24, gatk v=4.2.5.0, docker v=20.10.7, ' \
               'vcftools v=0.1.16'
        logging.info(info)
```

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
        # Checking version number for gatk,minimap2,samtools,bcftools and docker
        # Use this check to see if all the required tools in the log info have been installed
        Methods_Calling.check_version()

        # Creates the output folder
        Methods_Calling.folder_create(self.output_folder)
        Methods_Calling.folder_create(self.bam_folder)
        Methods_Calling.folder_create(self.pepper_folder)

        # Creates the bam file with read group tag for each query file in the list
        Methods_Calling.map_fastq_parallel(self.reference, self.fastq_list, self.bam_folder, self.cpu, self.parallel)

        # Produces the GVCF and VCF files using variant calling
        Methods_Calling.run_pepper_iter(self.reference, self.bam_folder, self.pepper_folder, self.cpu)

        # create a list of all the gvcf files
        g_vcf_list = Methods_Calling.list_file(self.pepper_folder, ".g.vcf.gz")

        # Renames each gcvf file with its unique accession number
        Methods_Calling.rename_sample(self.reference, g_vcf_list)

        # Merges all the new gcvf files to create a multi-sample gvcf and vcf
        g_vcf_rg_list = Methods_Calling.list_file(self.pepper_folder, "_rg.g.vcf.gz")
        Methods_Calling.merge_gvcf(self.reference, g_vcf_rg_list, self.cpu)

        # Reformat the multi-sample vcf file
        Methods_Calling.fix_merged_vcf('multi_sample.vcf', 'multi_sample_v1.vcf')

        # Filters through the new multi-sample vcf and selects sites that are bi-allelic, contain the
        # desired quality score etc. All sites that do not pass these requirements are removed
        # Adjust the filter settings accordingly
        Methods_Calling.filter_vcf('multi_sample_v1.vcf', 'multi_sample_final')

        # Converts the multi-sample vcf into an aligned fasta file
        Methods_Calling.Fasta(g_vcf_rg_list, 'multi_sample_v2.recode.vcf')

        # Creates a phylogenetic tree from the aligned fasta file
        Methods_Calling.Create_Phylo_tree('final.fasta')


class Methods_Calling:
    @staticmethod
    def check_version():
        # Checks the version number for all the tools to see if they have been installed
        docker_vers_check_cmd = ['docker', '--version']
        bcftools_vers_check_cmd = ['bcftools', '--version']
        samtools_vers_check_cmd = ['samtools', '--version']
        minimap2_vers_check_cmd = ['minimap2', '--version']
        vcftools_vers_check_cmd = ['vcftools', '-h']
        gatk_vers_check_cmd = ['gatk', '--version']
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
            print("folder is empty")
            exit()
        else:
            return file_list

    @staticmethod
    def folder_create(folder):
        # Create output folder
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
            for results in executor.map(lambda x: Methods_Calling.map_fastq(*x), args):
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
        sample_name = ((bam_name.split(".")[0]).split("_")[0])

        # Name of new gvcf file
        gvcf_file = sample_name + 'g.vcf'

        # Generate the gvcf and vcf files with pepper_deepvariant
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
        # docker_pull_cmd = ['docker', 'pull', 'kishwars/pepper_deepvariant:r0.7']
        # subprocess.run(docker_pull_cmd)
        bam_list = Methods_Calling.list_file(bam_folder, '.bam')

        # Stops the code if the bam files have not been generated
        if bam_list == list():
            print("bam folder is empty")
            exit()
        # Uses the run_pepper function to create gcvf and vcf files for each bam file
        else:
            for bam_file in bam_list:
                Methods_Calling.run_pepper(ref_file, bam_file, pepper_folder, cpu)

    @staticmethod
    def rename_sample(ref_file, g_vcf_list):
        # Gets the directory of the reference fasta file
        ref_file_cp = os.path.dirname(ref_file) + '/ref_cp.fasta'

        # Sample name of the ref file
        ref_base = os.path.basename(ref_file)

        # Saves a copy of the original reference fasta file
        Copy_ref_file = ['cp', ref_file, ref_file_cp]

        # Produces the reference file in a bgzip form
        bgzip_cmd = ['bgzip', ref_file]

        # Indexes the bgzip file
        faidx_cmd = ['samtools', 'faidx', (ref_file + '.gz')]

        # Run command to execute above commands
        subprocess.run(Copy_ref_file)
        p1 = subprocess.Popen(bgzip_cmd, stdout=subprocess.PIPE)
        p1.communicate()
        subprocess.run(faidx_cmd)

        # Rename the copied file to the original name
        ref_original = os.path.dirname(ref_file_cp) + '/' + ref_base
        Modify_ref_file = ['mv', ref_file_cp, ref_original]
        subprocess.run(Modify_ref_file)

        # Gvcf files not found if these files were not generated
        if g_vcf_list == list():
            print("gcvf files not found")
            exit()

        for file in g_vcf_list:
            # Changes the sample name to that of the accession number and indexes the files
            file_directory = os.path.dirname(file)
            sample_name = os.path.basename(file).split(".")[0]
            new_file = file_directory + '/' + sample_name + '_rg.g.vcf.gz'
            Change_sample_ID = ['gatk', "RenameSampleInVcf",
                                '-I', file,
                                '--NEW_SAMPLE_NAME', sample_name,
                                '-O', new_file]
            Index_new_file = ['tabix', new_file]

            subprocess.run(Change_sample_ID)
            subprocess.run(Index_new_file)

    @staticmethod
    def merge_gvcf(ref_file, g_v_list, cpu):
        # Merges the gvcf files into a single multi-sample gvcf and vcf
        merged_gvcf_file = 'multi_sample.g.vcf'
        merged_vcf_file = 'multi_sample.vcf'
        gvcf_create_cmd = ['bcftools',
                           'merge'] + g_v_list + ['-g', ref_file, '-O', "v", '-o', merged_gvcf_file]
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
                    out_file.write(line)

    @staticmethod
    def filter_vcf(vcf_in, vcf_out):
        # Creates the final version of the multi-sample vcf by filtering through the file
        filter_cmd = ['vcftools',
                      '--vcf', vcf_in,
                      '--out', vcf_out,
                      '--remove-indels',
                      '--maf', str(0.05),
                      '--min-alleles', str(2),
                      '--max-alleles', str(2),
                      '--minDP', str(20),
                      '--max-missing', str(0.5),
                      '--minQ', str(20),
                      '--remove-filtered-all',
                      '--recode']
        subprocess.run(filter_cmd)
    @staticmethod
    def IUPAC_base(ref, alt):
        # Adds the appropriate symbol for a heterzygous genotype (0/1)
        if ref in 'AG' and alt in 'AG':
            return 'R'
        if ref in 'CT' and alt in 'CT':
            return 'Y'
        if ref in 'GC' and alt in 'GC':
            return 'S'
        if ref in 'AT' and alt in 'AT':
            return 'W'
        if ref in 'GT' and alt in 'GT':
            return 'K'
        if ref in 'AC' and alt in 'AC':
            return 'M'

    @staticmethod
    def iterate_pandas(df, sample):
        # Iterates through the pandas vcf data frame to replace the genotypes with codons/gaps
        string = '<' + sample + '\n'
        for i in range(len(df)):
            if df.loc[i, sample][0:3] == ('./.'):
                string = string + '-'
            if df.loc[i, sample][0:3] == '0/1':
                string = string + Methods_Calling.IUPAC_base(df.loc[i, 'REF'], df.loc[i, 'ALT'])
            if df.loc[i, sample][0:3] == '0/0':
                string = string + df.loc[i, 'REF']
            if df.loc[i, sample][0:3] == '1/1':
                string = string + df.loc[i, 'ALT']
        return string

    @staticmethod
    def vcf_pandas(vcf):
        # Creates a pandas data frame from the created multi-sampple vcf
        with open(vcf, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]
        df = pd.read_csv(io.StringIO(''.join(lines)), dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str,
                                                          'ALT': str, 'QUAL': str, 'FILTER': str, 'INFO': str},
                                                          sep='\t').rename(columns={'#CHROM': 'CHROM'})
        return df



    @staticmethod
    def Fasta(gv_list, vcf):
        # Converts the multi-sample vcf to a fasta file
        df = Methods_Calling.vcf_pandas(vcf)
        s = ''
        L = list()
        for i in gv_list:
            L.append((os.path.basename(i).split(".")[0]).split("_")[0])
        for d in L:
            s = s + Methods_Calling.iterate_pandas(df,d) + '\n\n'
        with open('multi_sample.fasta', 'w') as output1:
            output1.write(s)
        Methods_Calling.Fasta_format_correction('multi_sample.fasta')

    @staticmethod
    def Fasta_format_correction(file):
        # fixes the formatting of the fasta file by specifying a certain amount of characters per line
        with open('final.fasta', 'w') as output2:
            with open(file, 'r') as input1:
                for line in input1:
                    if ('-' or 'A' or 'G' or 'T' or 'C') in line:
                        line = '\n'.join(line[i:i + 50] for i in range(0, len(line), 50))
                    output2.write(line)
        cmd = ['rm', '-f', file]
        subprocess.run(cmd)

    @staticmethod
    def Create_Phylo_tree(fasta):
        # Generates the Phylogenetic tree
        phylo_tree_cmd = ['raxmlHPC', '-s',
                          fasta, '-n', 'multi-sample_tree',
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

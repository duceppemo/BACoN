import glob
import os
import sys
import warnings
from argparse import ArgumentParser
from multiprocessing import cpu_count
from psutil import virtual_memory
from bacon_methods import Methods


__author__ = 'duceppemo'
__version__ = '0.1'


class Bacon(object):
    def __init__(self, args):
        # I/O
        self.reference = os.path.abspath(args.reference)
        self.input = os.path.abspath(args.input)
        self.output_folder = os.path.abspath(args.output)
        self.ref_size = args.size
        self.keep_bam = args.keep_bam
        self.baiting = args.baiting_method
        self.assembler = args.assembly_method
        self.min_size = args.min_size
        self.kmer = args.kmer_size
        self.snp_method = args.snp_method

        # Performance
        self.cpu = args.threads
        self.parallel = args.parallel
        self.mem = args.memory

        # Data
        self.sample_dict = dict()

        # Run
        self.run()

    def run(self):
        print('Checking a few things...')

        # Check if number of CPU and memory requested are valid
        self.cpu, self.parallel = Methods.check_cpus(self.cpu, self.parallel)
        self.mem = Methods.check_mem(self.mem)

        # Check input file compatibility
        Methods.check_input(self.input)

        ############################################################

        # Step completion report files
        done_trimming = self.output_folder + '/done_trimming'
        done_filtering = self.output_folder + '/done_filtering'
        done_extracting = self.output_folder + '/done_extracting'
        done_assembling = self.output_folder + '/done_assembling'
        done_comparing = self.output_folder + '/done_comparing'

        # Output folders to create
        extracted_folder = self.output_folder + '/1_extracted/'
        trimmed_folder = self.output_folder + '/2_trimmed/'
        filtered_folder = self.output_folder + '/3_filtered/'
        assembled_folder = self.output_folder + '/4_assembled/'
        compared_folder = self.output_folder + '/5_compared/'

        # Create output folder
        Methods.make_folder(self.output_folder)

        # Get input files and place info in dictionary
        if os.path.isdir(self.input):
            self.sample_dict['raw'] = Methods.get_files(self.input)
        else:
            sample = os.path.basename(self.input).split('.')[0].replace('_pass', '')
            self.sample_dict['raw'] = {sample: os.path.realpath(self.input)}

        if len(self.sample_dict['raw']) == 1:
            print('Only one sample to process. All threads will be used for that sample.')
            self.parallel = 1  # Use all cpu for single sample.

        # Drop the unclassified sample
        # self.sample_dict['raw'].pop('unclassified_pass', None)

        print('\tAll good!')

        ##################
        #
        # 1- Bait reads
        #
        ##################

        if not os.path.exists(done_extracting):
            if self.baiting == 'minimap2':
                print('Extracting reads matching {} with Minimap2...'.format(os.path.basename(self.reference)))
                Methods.run_minimap2_parallel(extracted_folder, self.reference, self.sample_dict['raw'],
                                              self.cpu, self.parallel, self.keep_bam)
            else:  # if self.baiting == 'bbduk':
                print('Extracting reads matching {} with BBduk...'.format(os.path.basename(self.reference)))
                Methods.bait_bbduk_parallel(extracted_folder, self.reference, self.sample_dict['raw'],
                                            self.cpu, self.parallel, self.kmer, self.mem)
            Methods.flag_done(done_extracting)
        else:
            print('Skipping extracting. Already done.')

        # Update sample_dict after extracting
        self.sample_dict['extracted'] = Methods.get_files(extracted_folder)

        ##################
        #
        # 2- Trim reads
        #
        ##################

        if not os.path.exists(done_trimming):
            print('Removing Nanopore adapters with Porechop...')
            Methods.run_porechop_parallel(self.sample_dict['extracted'], trimmed_folder, self.cpu, self.parallel)
            Methods.flag_done(done_trimming)
        else:
            print('Skipping trimming. Already done.')

        # Update sample_dict after trimming
        self.sample_dict['trimmed'] = Methods.get_files(trimmed_folder)

        ##################
        #
        # 3- Filter reads
        #
        ##################

        # Get reference size
        if not self.ref_size:
            self.ref_size = Methods.fasta_length(self.reference)
        ref_name = '.'.join(os.path.basename(self.reference).split('.')[:-1])
        print('Reference ({}): {} bp'.format(ref_name, self.ref_size))

        if not os.path.exists(done_filtering):
            print('Filtering lower quality reads with Filtlong...')
            Methods.run_filtlong_parallel(self.sample_dict['trimmed'], filtered_folder,
                                          self.ref_size, self.parallel)
            Methods.flag_done(done_filtering)
        else:
            print('Skipping filtering. Already done.')

        # Update sample_dict after trimming
        self.sample_dict['filtered'] = Methods.get_files(filtered_folder)

        ##################
        #
        # 4- Assemble reads
        #
        ##################

        if not os.path.exists(done_assembling):
            if self.assembler == 'flye':
                print('Assembling extracted reads with Flye...')

                Methods.assemble_flye_parallel(self.sample_dict['filtered'], assembled_folder,
                                               self.ref_size, self.min_size, self.cpu, self.parallel)
                Methods.flye_assembly_stats(assembled_folder, self.output_folder)  # Get stats
            elif self.assembler == 'shasta':
                print('Assembling extracted reads with Shasta...')

                Methods.assemble_shasta_parallel(self.sample_dict['filtered'], assembled_folder,
                                                 self.min_size, self.cpu, self.parallel)
                Methods.shasta_assembly_stats(assembled_folder, self.output_folder)  # Get stats
            else:  # elif self.assembler == 'rebaler':
                print('Performing reference-guided assembly with Rebaler...')

                Methods.assemble_rebaler_parallel(self.reference, self.sample_dict['filtered'], assembled_folder,
                                                  self.cpu, self.parallel)
                # Methods.shasta_assembly_stats(assembled_folder, self.output_folder)  # Get stats
            # Completion flag
            Methods.flag_done(done_assembling)
        else:
            print('Skipping assembling organelle genomes. Already done.')

        # Update sample_dict after trimming
        self.sample_dict['assembled'] = Methods.get_files(assembled_folder)
        assembly_list = Methods.list_files_in_folder(assembled_folder + 'all_assemblies/', "*.fasta")

        ##################
        #
        # 5- Compare samples
        #
        ##################

        # Make stats and create SNP VCF file (populations)
        if not os.path.exists(done_comparing):
            if len(assembly_list) < 3:
                warnings.warn('Cannot build a tree with less than three samples!')
                sys.exit()

            if self.snp_method == 'parsnp':
                print('Making core-SNP tree with Parsnp...')

                # Create core-SNP tree
                Methods.run_parsnp(assembly_list, compared_folder, self.reference, self.cpu)

                # Check if Parsnp ran to completion. Sometimes when assembly sizes are too different from the ref
                # Parsnp won't run to completion
                if not os.path.exists(compared_folder + 'parsnp.xmfa'):
                    raise Exception('Parsnp could not run. Sample comparison not be performed.')

                # Plot tree to PDF, .SVG or .PNG
                Methods.plot_newick_tree(compared_folder + 'parsnp.tree', compared_folder + 'parsnp.pdf')

                # Convert xmfa to fasta to compute a new tree with bootstraps
                Methods.convert_xmfa_to_fastq(compared_folder + 'parsnp.xmfa', compared_folder + 'parsnp.fasta')

                # Create ML tree with bootstraps
                print('Making RAxML tree with 1,000 bootstraps...')
                Methods.make_tree_raxml(compared_folder + 'parsnp.fasta', compared_folder, self.cpu)
                # Plot tree to PDF, .SVG or .PNG
                # This command won't work for this tree...
                # Methods.plot_newick_tree(compared_folder + 'RAxML_bipartitionsBranchLabels.raxml.tree',
                #                          compared_folder + 'raxml.pdf')

                # Create FastTree
                print('Making FastTree tree with 1,000 peudo-bootstraps...')
                Methods.make_tree_fasttree(compared_folder + 'parsnp.fasta', compared_folder, self.cpu)
                Methods.plot_newick_tree(compared_folder + 'fasttree.tree',
                                         compared_folder + 'fasttree.pdf')
            else:  # elif self.snp_method == 'ksnp3':
                print('Making pan-SNP ans core-SNP trees with kSNP3...')

                # Create pan-SNP tree
                Methods.run_ksnp3(assembly_list, compared_folder, self.reference, self.cpu)

            Methods.flag_done(done_comparing)
        else:
            print('Skipping sample comparison. Already done.')

        print('DONE!')


if __name__ == "__main__":
    max_cpu = cpu_count()
    max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB

    parser = ArgumentParser(description='Extract, assemble and compare Nanopore reads matching a reference sequence.')
    parser.add_argument('-r', '--reference', metavar='/path/to/reference_organelle/genome.fasta',
                        required=True,
                        help='Reference genome for read mapping. Mandatory.')
    parser.add_argument('-i', '--input', metavar='/path/to/input/folder/ or /path/to/my_fastq.gz',
                        required=True,
                        help='Folder that contains the fastq files or individual fastq file. Mandatory.')
    parser.add_argument('-o', '--output', metavar='/path/to/output/folder/',
                        required=True,
                        help='Folder to hold the result files. Mandatory.')
    parser.add_argument('-b', '--baiting-method',
                        required=False, default='minimap2',
                        choices=['minimap2', 'bbduk'],
                        type=str,
                        help='Read baiting method. Default "minimap2". Optional.')
    parser.add_argument('-a', '--assembly-method',
                        required=False, default='flye',
                        choices=['flye', 'shasta', 'rebaler'],
                        type=str,
                        help='Assembly method. Default "flye". Optional.')
    parser.add_argument('--min-size', metavar='3000',
                        required=False, default='3000',
                        type=int,
                        help='Minimum read size for Shasta assembler or minimum read overlap for Flye. '
                             'Default 3000 for Shasta and auto for Flye. Optional.')
    parser.add_argument('-t', '--threads', metavar=str(max_cpu),
                        required=False,
                        type=int, default=max_cpu,
                        help='Number of threads. Default is maximum available({}). Optional.'.format(max_cpu))
    parser.add_argument('-p', '--parallel', metavar='2',
                        required=False,
                        type=int, default=2,
                        help='Number of samples to process in parallel. Default is 2. Optional.')
    parser.add_argument('-m', '--memory', metavar=str(max_mem),
                        required=False,
                        type=int, default=max_mem,
                        help='Memory in GB. Default is 85%% of total memory ({})'.format(max_mem))
    parser.add_argument('-s', '--size', metavar='150000',
                        required=False,
                        type=int,
                        help='Override automatically detected reference size. Optional.')
    parser.add_argument('-k', '--kmer-size', metavar='99',
                        required=False,
                        type=int, default=99,
                        help='Kmer size for baiting. Only used if "--baiting-method" is "bbduk". Optional.')
    parser.add_argument('--keep-bam',
                        required=False, action='store_true',
                        help='Do not delete BAM files. Only used if "--baiting-method" is "minimap2". Optional.')
    parser.add_argument('-snp', '--snp-method',
                        required=False, default='parsnp',
                        choices=['parsnp', 'ksnp3'],
                        type=str,
                        help='SNP calling method. Default "parsnp". Optional.')
    parser.add_argument('-v', '--version', action='version',
                        version=f'{os.path.basename(__file__)}: version {__version__}')

    # Get the arguments into an object
    arguments = parser.parse_args()

    Bacon(arguments)

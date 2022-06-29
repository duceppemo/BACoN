import glob
import os
import sys
import warnings
from argparse import ArgumentParser
from multiprocessing import cpu_count
from bacon_methods import Methods


__author__ = 'duceppemo'
__version__ = 'v0.1'


class GBS(object):
    def __init__(self, args):
        # I/O
        self.reference = os.path.abspath(args.reference)
        self.input_folder = os.path.abspath(args.input)
        self.output_folder = os.path.abspath(args.output)
        self.ref_size = args.size

        # Performance
        self.cpu = args.threads
        self.parallel = args.parallel

        # Data
        self.sample_dict = dict()

        # Run
        self.run()

    def run(self):
        print('Checking a few things...')

        # Check if number of CPU and memory requested are valid
        self.cpu, self.parallel = Methods.check_cpus(self.cpu, self.parallel)

        # Check input file compatibility
        error_message = Methods.check_input_folder(self.input_folder)
        if error_message:
            raise Exception(error_message)

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
        self.sample_dict['raw'] = Methods.get_files(self.input_folder)
        # Drop the unclassified sample
        # self.sample_dict['raw'].pop('unclassified_pass', None)

        ##################
        #
        # 1- Extract reads by mapping
        #
        ##################

        if not os.path.exists(done_extracting):
            print('Extracting organelle reads with Minimap2...')
            Methods.run_minimap2_parallel(extracted_folder, self.reference, self.sample_dict['raw'],
                                          self.cpu, self.parallel)
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
        print('Reference ({}): {}bp'.format(ref_name, self.ref_size))

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
            print('Assembling extracted reads with Flye...')

            Methods.assemble_flye_parallel(self.sample_dict['filtered'], assembled_folder,
                                           self.ref_size, self.cpu, self.parallel)
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
            print('Making core SNP tree with Parsnp...')

            # Create Cores SNP tree
            Methods.run_parsnp(assembly_list, compared_folder, self.reference, self.cpu)

            # Plot tree to PDF, .SVG or .PNG
            Methods.plot_newick_tree(compared_folder + 'parsnp.tree', compared_folder + 'parsnp.pdf')

            # Convert xmfa to fasta so we can compute a new tree with bootstraps
            Methods.convert_xmfa_to_fastq(compared_folder + 'parsnp.xmfa', compared_folder + 'parsnp.fasta')

            # Create ML tree with bootstraps
            print('Making tree with 1,000 bootstraps...')
            Methods.make_tree_raxml(compared_folder + 'parsnp.fasta', compared_folder, self.cpu)
            # Plot tree to PDF, .SVG or .PNG
            # This command wont work for this tree...
            # Methods.plot_newick_tree(compared_folder + 'RAxML_bipartitionsBranchLabels.raxml.tree',
            #                          compared_folder + 'raxml.pdf')

            # Create FastTree
            Methods.make_tree_fasttree(compared_folder + 'parsnp.fasta', compared_folder, self.cpu)
            Methods.plot_newick_tree(compared_folder + 'fasttree.tree',
                                     compared_folder + 'fasttree.pdf')

            Methods.flag_done(done_comparing)
        else:
            print('Skipping sample comparison. Already done.')


if __name__ == "__main__":
    max_cpu = cpu_count()
    # max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB

    parser = ArgumentParser(description='Organelle genome assembly from Nanopore genome skimming data.')
    parser.add_argument('-r', '--reference', metavar='/path/to/reference_organelle/genome.fasta',
                        required=True,
                        help='Reference genome for read mapping. Mandatory.')
    parser.add_argument('-i', '--input', metavar='/path/to/input/folder/',
                        required=True,
                        help='Folder that contains the read files. Mandatory.')
    parser.add_argument('-o', '--output', metavar='/path/to/output/folder/',
                        required=True,
                        help='Folder to hold the result files. Mandatory.')
    parser.add_argument('-t', '--threads', metavar=str(max_cpu),
                        required=False,
                        type=int, default=max_cpu,
                        help='Number of threads. Default is maximum available({}). Optional'.format(max_cpu))
    parser.add_argument('-p', '--parallel', metavar='2',
                        required=False,
                        type=int, default=2,
                        help='Number of samples to process in parallel. Default is 2. Optional')
    parser.add_argument('-s', '--size', metavar='150000',
                        required=False,
                        type=int,
                        help='Override automatically detected reference size. Optional')

    # Get the arguments into an object
    arguments = parser.parse_args()

    GBS(arguments)

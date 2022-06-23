import subprocess
import os
import sys
from concurrent import futures
import pathlib
from shutil import move
from multiprocessing import cpu_count
import gzip
from itertools import groupby
from ete3 import Tree, TreeStyle
from ete3.parser.newick import NewickError
from glob import glob
import pysam


class Methods(object):
    accepted_extensions = ['.fq', '.fq.gz',
                           '.fastq', '.fastq.gz',
                           '.fasta', '.fasta.gz',
                           '.fa', '.fa.gz',
                           '.fna', '.fna.gz',
                           '.bam']

    @staticmethod
    def check_input_folder(folder):
        error_message_list = ['Please provide a folder as input.',
                              'The input folder provided does not exist.',
                              'Make sure files in input folder end with {}'.format(Methods.accepted_extensions)
                              ]

        # Check if input folder exists and is a folder
        isdir_status = os.path.isdir(folder)
        exists_status = os.path.exists(folder)

        # List content of folder
        file_list = os.listdir(folder)

        # if folder is not empty and all files have the accepted extensions
        extension_status_list = [f.endswith(tuple(Methods.accepted_extensions)) for f in file_list]

        if not isdir_status:
            return error_message_list[0]
        elif not exists_status:
            return error_message_list[1]
        elif not any(extension_status_list):
            return error_message_list[2]
        else:
            return

    @staticmethod
    def check_cpus(requested_cpu, n_proc):
        total_cpu = cpu_count()

        if 1 > requested_cpu > total_cpu:
            requested_cpu = total_cpu
            sys.stderr.write("Number of threads was set to {}".format(requested_cpu))
        if 1 > n_proc > total_cpu:
            n_proc = 2
            sys.stderr.write("Number of threads was set to {}".format(2))

        return requested_cpu, n_proc

    @staticmethod
    def check_version(log_file):
        # Not being used right now because versions are captured in the requirements.txt file
        with open(log_file, 'w') as f:
            # Python
            p = subprocess.Popen(['python', '--version'])
            stderr, stdout = p.communicate()

            # Porechop
            p = subprocess.Popen(['porechop', '--version'])
            stderr, stdout = p.communicate()

            # bbmap suite
            p = subprocess.Popen(['bbduk.sh', '--version'])
            stderr, stdout = p.communicate()

            # Filtlong
            p = subprocess.Popen(['filtlong', '--version'])
            stderr, stdout = p.communicate()

            # SAMtools
            p = subprocess.Popen(['samtools', '--version'])
            stderr, stdout = p.communicate()

            # Minimap2
            p = subprocess.Popen(['minimap2', '--version'])
            stderr, stdout = p.communicate()

    @staticmethod
    def make_folder(folder):
        # Will create parent directories if don't exist and will not return error if already exists
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def get_files(in_folder):
        sample_dict = dict()

        # Look for input sequence files recursively
        for root, directories, filenames in os.walk(in_folder):
            for filename in filenames:
                if filename.endswith(tuple(Methods.accepted_extensions)):  # accept a tuple or string
                    file_path = os.path.join(root, filename)
                    sample = filename.split('.')[0].replace('_pass', '')
                    if filename.endswith('gz'):
                        sample = sample.split('.')[0]
                    sample_dict[sample] = file_path
        if not sample_dict:
            raise Exception('Sample dictionary empty!')

        return sample_dict

    @staticmethod
    def list_to_file(my_list, output_file):
        with open(output_file, 'wt') as f:
            for l in my_list:
                f.write('{}\n'.format(l))

    @staticmethod
    def list_files_in_folder(folder, extension):
        return glob(folder + '/*' + extension)

    @staticmethod
    def flag_done(flag_file):
        with open(flag_file, 'w') as f:
            pass

    @staticmethod
    def convert_bam_to_fastq(sample, bam_file, output_folder):
        """Converting the minimap2 BAM files directly to fasta with reformat.sh results in
        more than one sequence with the ID, so I had to write a custom function to extrac the mapped reads."""
        cmd = ['reformat.sh',
               'ow=t',
               'in={}'.format(bam_file),
               'out={}'.format(output_folder + sample + '.fastq.gz')]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def get_fastq_from_bam(sample, bam_file, fastq_file, output_folder):
        read_dict = dict()
        with pysam.AlignmentFile(bam_file, 'rb') as f:
            for read in f.fetch():
                read_name = read.qname
                read_dict[read_name] = ''  # push reads into dict to avoid duplicates
                # if read_name not in read_dict:
                #     read_dict[read_name] = ''  # push reads into dict to avoid duplicates

        extracted_fastq = output_folder + sample + '.fastq.gz'

        # Parse fastq file to dictionary
        line_counter = 0
        with gzip.open(extracted_fastq, 'wb') as out_fh:
            with gzip.open(fastq_file, 'rt') if fastq_file.endswith('.gz') else open(fastq_file, 'r') as f:
                fastq_entry = list()
                for line in f:
                    line_counter += 1
                    fastq_entry.append(line)
                    if line_counter == 4:  # last line of fastq entry

                        # Ditch the leading "@" and everything after 1st space
                        seq_id = fastq_entry[0].split()[0][1:]

                        # Write to
                        if seq_id in read_dict:
                            out_fh.write(''.join(fastq_entry).encode('ascii'))

                        # Prepare for new fastq entry
                        line_counter = 0
                        fastq_entry = list()

    @staticmethod
    def run_minimap2(sample, ref, fastq_file, cpu, output_folder):
        print('\t{}'.format(sample))

        output_bam = output_folder + sample + '.bam'

        minimap2_cmd = ['minimap2',
                        '-a',
                        '-x', 'map-ont',
                        '-t', str(cpu),
                        '--MD',
                        '--secondary=no',
                        ref,
                        fastq_file]
        samtools_view_cmd = ['samtools', 'view',
                             '-@', str(cpu),
                             '-F', '4', '-h',
                             '-T', ref,
                             '-']
        samtools_sort_cmd = ['samtools', 'sort',
                             '-@', str(cpu),
                             '--reference', ref,
                             '-']
        samtools_markdup_cmd = ['samtools', 'markdup',
                                '-r',
                                '-@', str(cpu),
                                '-',
                                output_bam]
        # samtools can only index chromosomes up to 512M bp.
        samtools_index_cmd = ['samtools', 'index',
                              output_bam]

        p1 = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2 = subprocess.Popen(samtools_view_cmd, stdin=p1.stdout, stdout=subprocess.PIPE)#, stderr=subprocess.DEVNULL)
        p1.stdout.close()
        p3 = subprocess.Popen(samtools_sort_cmd, stdin=p2.stdout, stdout=subprocess.PIPE)#, stderr=subprocess.DEVNULL)
        p2.stdout.close()
        p4 = subprocess.Popen(samtools_markdup_cmd, stdin=p3.stdout, stdout=subprocess.PIPE)#, stderr=subprocess.DEVNULL)
        p3.stdout.close()
        p4.communicate()

        # Index bam file
        subprocess.run(samtools_index_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        # Convert bam to fastq
        Methods.get_fastq_from_bam(sample, output_bam, fastq_file, output_folder)
        # Methods.convert_bam_to_fastq(sample, output_bam, output_folder)

        # Remove bam
        bam_list = glob(output_folder + '*.bam*')
        for bam in bam_list:
            os.remove(bam)

    @staticmethod
    def run_minimap2_parallel(output_folder, ref, sample_dict, cpu, parallel):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, ref, path, int(cpu / parallel), output_folder)
                    for sample, path in sample_dict.items())
            for results in executor.map(lambda x: Methods.run_minimap2(*x), args):
                pass

    @staticmethod
    def run_porechop(sample, input_fastq, trimmed_folder, cpu):
        cmd = ['porechop',
               '-i', input_fastq,
               '-o', trimmed_folder + '/' + sample + '.fastq.gz',
               '--threads', str(cpu),
               '--check_reads', str(1000)]

        print('\t{}'.format(sample))
        subprocess.run(cmd)#, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def run_porechop_parallel(sample_dict, output_folder, cpu, parallel):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, path, output_folder, int(cpu / parallel))
                    for sample, path in sample_dict.items())
            for results in executor.map(lambda x: Methods.run_porechop(*x), args):
                pass

    @staticmethod
    def run_filtlong(sample, input_fastq, filtered_folder, genome_size):
        print('\t{}'.format(sample))

        cmd = ['filtlong',
               '--keep_percent', str(95),  # drop bottom 5% reads
               '--target_bases', str(genome_size * 100),  # keep top 100X if more reads
               input_fastq]

        # Filtlong writes to stdout
        filtered_fastq = filtered_folder + sample + '.fastq.gz'
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        with gzip.open(filtered_fastq, 'wb') as f:
            f.write(p.communicate()[0])

    @staticmethod
    def run_filtlong_parallel(sample_dict, output_folder, genome_size, parallel):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, path, output_folder, genome_size)
                    for sample, path in sample_dict.items())
            for results in executor.map(lambda x: Methods.run_filtlong(*x), args):
                pass

    @staticmethod
    def fasta_length(input_fasta):
        with gzip.open(input_fasta, 'rt') if input_fasta.endswith('.gz') else open(input_fasta, 'r') as f:

            # Create iterator in case there are more than 1 contig in reference genome
            faiter = (x[1] for x in groupby(f, lambda line: line[0] == '>'))

            total_len = 0
            for header in faiter:
                # Join all sequence lines of fasta entry in one line, measure its length
                # and add it to length of any other previous sequences if present
                total_len += len(''.join(s.rstrip() for s in faiter.__next__()))

            return total_len

    @staticmethod
    def assemble_flye(sample, input_fastq, output_folder, genome_size, cpu):
        print('\t{}'.format(sample))

        # Create a subfolder for each sample
        output_subfolder = output_folder + sample + '/'
        Methods.make_folder(output_subfolder)

        cmd_flye = ['flye',
                    '--genome-size', str(genome_size),
                    '--nano-hq', input_fastq,
                    '--threads', str(cpu),
                    '--out-dir', output_subfolder,
                    '--iterations', str(3)]
        subprocess.run(cmd_flye)

        # Rename and move assembly file
        assemblies_folder = output_folder + 'all_assemblies/'
        Methods.make_folder(assemblies_folder)
        move(output_subfolder + 'assembly.fasta', assemblies_folder + sample + '.fasta')

        # Rename and move assembly graph file
        assembly_graph_folder = output_folder + 'assembly_graphs/'
        Methods.make_folder(assembly_graph_folder)
        move(output_subfolder + 'assembly_graph.gfa', assembly_graph_folder + sample + '_graph.gfa')

        # Assembly graph
        cmd_bandage = ['Bandage', 'image',
                       '{}'.format(assembly_graph_folder + sample + '_graph.gfa'),
                       '{}'.format(assembly_graph_folder + sample + '_graph.png')]
        subprocess.run(cmd_bandage)

    @staticmethod
    def assemble_flye_parallel(sample_dict, output_folder, genome_size, cpu, parallel):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, path, output_folder, genome_size, int(cpu / parallel))
                    for sample, path in sample_dict.items())
            for results in executor.map(lambda x: Methods.assemble_flye(*x), args):
                pass

    @staticmethod
    def run_parsnp(assembly_list, output_folder, ref, cpu):
        Methods.make_folder(output_folder)

        cmd = ['parsnp',
               '-r', ref,
               '-p', str(cpu),
               '-o', output_folder,
               '-d'] + assembly_list
        subprocess.run(cmd)

    @staticmethod
    def convert_xmfa_to_fastq(xmfa_input, fasta_output):
        cmd = ['harvesttools',
               '-x', xmfa_input,
               '-M', fasta_output]
        subprocess.run(cmd)

    @staticmethod
    def make_tree_raxml(aligned_fasta, output_folder, cpu):
        cmd = ['raxmlHPC-PTHREADS-AVX2',
               '-s', aligned_fasta,
               '-w', output_folder,
               '-n', 'raxml.tree',
               '-m', 'GTRCAT',
               '-N', str(1000),
               '-d', '-f', 'a', '-T', str(cpu),
               '-x', str(1234), '-p', str(123)]
        subprocess.run(cmd)

    @staticmethod
    def make_tree_fasttree(aligned_fasta, output_folder, cpu):
        cmd = ['FastTree',
               '-nt',
               '-gtr',
               aligned_fasta]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        with open(output_folder + 'fasttree.tree', 'wb') as f:
            f.write(p.communicate()[0])

    @staticmethod
    def plot_newick_tree(tree_file, output_file):
        try:
            t = Tree(tree_file)
            ts = TreeStyle()
            ts.show_leaf_name = True
            ts.show_branch_support = True
            t.render(output_file, w=183, units='mm', tree_style=ts)
        except NewickError:
            print('Could not convert tree to picture file.')

import subprocess
import os
import sys
from concurrent import futures
import pathlib
from shutil import move
from psutil import virtual_memory
from multiprocessing import cpu_count
import gzip
from itertools import groupby
from ete3 import Tree, TreeStyle
from ete3.parser.newick import NewickError
from glob import glob
import pysam
import warnings
import pandas as pd
import io


class Methods(object):
    accepted_extensions = ['.fq', '.fq.gz',
                           '.fastq', '.fastq.gz',
                           '.fasta', '.fasta.gz',
                           '.fa', '.fa.gz',
                           '.fna', '.fna.gz']

    @staticmethod
    def check_input(my_input):
        error_message_list = ['Please provide a folder as input.',
                              'The input folder provided does not exist.',
                              'Make sure files in input folder end with {}'.format(Methods.accepted_extensions),

                              ]

        if not os.path.exists(my_input):
            raise Exception('Please select an existing file or folder as input.')

        # Check if folder
        if os.path.isdir(my_input):
            file_list = os.listdir(my_input)  # List content of folder
        else:  # elif os.path.isfile(my_input):
            file_list = [my_input]

        # if folder is not empty and all files have the accepted extensions
        if not all([f.endswith(tuple(Methods.accepted_extensions)) for f in file_list]):
            raise Exception('Make sure files in input folder end with {}'.format(Methods.accepted_extensions))

    @staticmethod
    def check_cpus(requested_cpu, n_proc):
        total_cpu = cpu_count()

        if 1 > requested_cpu > total_cpu:
            requested_cpu = total_cpu
            sys.stderr.write("Number of threads was set to {}".format(requested_cpu))
        if 1 > n_proc > total_cpu:
            n_proc = total_cpu
            sys.stderr.write("Number of samples to parallel process was set to {}".format(total_cpu))

        return requested_cpu, n_proc

    @staticmethod
    def check_mem(requested_mem):
        max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB
        if requested_mem:
            if requested_mem > max_mem:
                requested_mem = max_mem
                sys.stderr.write("Requested memory was set higher than available system memory ({})".format(max_mem))
                sys.stderr.write("Memory was set to {}".format(requested_mem))
        else:
            requested_mem = max_mem

        return requested_mem

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
                    file_path = os.path.realpath(file_path)  # follow symbolic links
                    sample = filename.split('.')[0].replace('_pass', '').replace('_filtered', '')
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
    def gzipped_file_size(gzipped_file):
        with gzip.open(gzipped_file, 'rb') as f:
            return f.seek(0, whence=2)

    @staticmethod
    def run_minimap2(sample, ref, fastq_file, cpu, output_folder, keep_bam):
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
        p2 = subprocess.Popen(samtools_view_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p1.stdout.close()
        p3 = subprocess.Popen(samtools_sort_cmd, stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2.stdout.close()
        p4 = subprocess.Popen(samtools_markdup_cmd, stdin=p3.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p3.stdout.close()
        p4.communicate()

        # Index bam file
        if os.path.exists(output_bam):
            if os.stat(output_bam).st_size != 0:  # bam file exists and not empty
                subprocess.run(samtools_index_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

                # Convert bam to fastq
                Methods.get_fastq_from_bam(sample, output_bam, fastq_file, output_folder)

                # Remove bam
                if not keep_bam:
                    bam_list = glob(output_folder + '*.bam*')
                    for bam in bam_list:
                        os.remove(bam)
        else:
            warnings.warn('No reads were extracted for {}!'.format(sample))

    @staticmethod
    def run_minimap2_parallel(output_folder, ref, sample_dict, cpu, parallel, keep_bam):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, ref, path, int(cpu / parallel), output_folder, keep_bam)
                    for sample, path in sample_dict.items())
            for results in executor.map(lambda x: Methods.run_minimap2(*x), args):
                pass

    @staticmethod
    def bait_bbduk(sample, ref, fastq_file, cpu, output_folder, kmer, mem):
        print('\t{}'.format(sample))

        extracted_fastq = output_folder + sample + '.fastq.gz'

        cmd = ['bbduk.sh',
               '-Xmx{}g'.format(mem),
               '-eoom',
               'overwrite=true',
               'in={}'.format(fastq_file),
               'ref={}'.format(ref),
               'threads={}'.format(cpu),
               'k={}'.format(kmer),
               'maskmiddle=f',
               'qin=33',
               'hdist=2',
               'outm={}'.format(extracted_fastq)]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        # Check if gzipped fastq is not empty
        if Methods.gzipped_file_size(extracted_fastq) == 0:
            warnings.warn('No reads were extracted for {}!'.format(sample))

    @staticmethod
    def bait_bbduk_parallel(output_folder, ref, sample_dict, cpu, parallel, kmer, mem):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, ref, path, int(cpu / parallel), output_folder, kmer, int(mem / parallel))
                    for sample, path in sample_dict.items())
            for results in executor.map(lambda x: Methods.bait_bbduk(*x), args):
                pass

    @staticmethod
    def run_porechop(sample, input_fastq, trimmed_folder, cpu):
        cmd = ['porechop',
               '-i', input_fastq,
               '-o', trimmed_folder + '/' + sample + '.fastq.gz',
               '--threads', str(cpu),
               '--check_reads', str(1000)]

        print('\t{}'.format(sample))
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

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
               '--min_length', str(500),  #  remove rejected reads from targeted sequencing
               input_fastq]

        # Filtlong writes to stdout
        filtered_fastq = filtered_folder + sample + '.fastq.gz'
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
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
    def assemble_flye(sample, input_fastq, output_folder, genome_size, min_size, cpu):
        print('\t{}'.format(sample))

        # Create a subfolder for each sample
        output_subfolder = output_folder + sample + '/'
        Methods.make_folder(output_subfolder)

        if min_size:
            cmd_flye = ['flye',
                        '--genome-size', str(genome_size),
                        '--nano-hq', input_fastq,
                        '--threads', str(cpu),
                        '--out-dir', output_subfolder,
                        '--iterations', str(3),
                        '--min-overlap', str(min_size)]
        else:
            cmd_flye = ['flye',
                        '--genome-size', str(genome_size),
                        '--nano-hq', input_fastq,
                        '--threads', str(cpu),
                        '--out-dir', output_subfolder,
                        '--iterations', str(3)]
        subprocess.run(cmd_flye, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        # Rename and move assembly file
        assemblies_folder = output_folder + 'all_assemblies/'
        Methods.make_folder(assemblies_folder)
        if os.path.exists(output_subfolder + 'assembly.fasta'):
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
        else:
            warnings.warn('No assembly for {}!'.format(sample))

    @staticmethod
    def assemble_flye_parallel(sample_dict, output_folder, genome_size, min_size, cpu, parallel):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, path, output_folder, genome_size, min_size, int(cpu / parallel))
                    for sample, path in sample_dict.items())
            for results in executor.map(lambda x: Methods.assemble_flye(*x), args):
                pass

    @staticmethod
    def flye_assembly_stats(assembly_folder, output_folder):
        # Output file
        output_stats_file = output_folder + '/flye_stats.tsv'

        # Pandas data frame to save values
        df = pd.DataFrame(columns=['Sample', 'TotalBases', 'ReadsN50', 'AssemblyLength', 'Contigs', 'Coverage'])

        # Find log file(s) and parse info of interest
        log_list = glob(assembly_folder + '/**/flye.log', recursive=True)

        for log_file in log_list:
            sample = log_file.split('/')[-2]
            with open(log_file, 'r') as f:
                read_len = 0
                n50 = 0
                assembly_len = 0
                contigs = 0
                mean_cov = 0

                for line in f:
                    line = line.rstrip()
                    if 'Total read length:' in line:
                        read_len = line.split()[-1]
                    elif 'N50/N90' in line:
                        n50 = line.split(':')[-1].split('/')[0].strip()
                    elif 'Total length:' in line:
                        assembly_len = line.split('\t')[-1]
                    elif 'Fragments:' in line:
                        contigs = line.split('\t')[-1]
                    elif 'Mean coverage:' in line:
                        mean_cov = line.split('\t')[-1]

                data_dict = {'Sample': [sample],
                             'TotalBases': [read_len],
                             'ReadsN50': [n50],
                             'AssemblyLength': [assembly_len],
                             'Contigs': [contigs],
                             'Coverage':  [mean_cov]}

                df = pd.concat([df, pd.DataFrame.from_dict(data_dict)], axis='index', ignore_index=True)

        # Convert df to tsv file
        df.to_csv(output_stats_file, sep='\t', index=False)

    @staticmethod
    def assemble_shasta(sample, input_fastq, output_folder, min_size, cpu):
        print('\t{}'.format(sample))

        # Create a subfolder for each sample
        output_subfolder = output_folder + sample + '/'

        # Unzipped fastq file (needed for shasta)
        unzipped_fastq = output_folder + sample + '.fastq'

        cmd_ungzip = ['pigz', '-dkc', input_fastq]  # To stdout
        cmd_shasta = ['shasta',
                      '--config', 'Nanopore-Oct2021',
                      '--input', unzipped_fastq,
                      '--assemblyDirectory', output_subfolder,
                      '--command', 'assemble',
                      '--threads', str(cpu),
                      '--Reads.minReadLength', str(min_size)]
        cmd_shasta_clean = ['shasta',
                            '--assemblyDirectory', output_subfolder,
                            '--command', 'cleanupBinaryData']

        # Decompress fastq for shasta
        with open(unzipped_fastq, 'w') as f:
            subprocess.run(cmd_ungzip, stdout=f)

        # Run shasta assembler
        shasta_stdout = output_folder + sample + '_shasta_stdout.txt'
        with open(shasta_stdout, 'w') as f:
            subprocess.run(cmd_shasta, stdout=f, stderr=subprocess.DEVNULL)

        # Cleanup temporary files
        subprocess.run(cmd_shasta_clean, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        os.remove(unzipped_fastq)

        # Rename and move assembly file
        assemblies_folder = output_folder + 'all_assemblies/'
        Methods.make_folder(assemblies_folder)
        if os.path.exists(output_subfolder + 'Assembly.fasta'):
            move(output_subfolder + 'Assembly.fasta', assemblies_folder + sample + '.fasta')

            # Rename and move assembly graph file
            assembly_graph_folder = output_folder + 'assembly_graphs/'
            Methods.make_folder(assembly_graph_folder)
            move(output_subfolder + 'Assembly.gfa', assembly_graph_folder + sample + '.gfa')

            # Assembly graph
            cmd_bandage = ['Bandage', 'image',
                           '{}'.format(assembly_graph_folder + sample + '.gfa'),
                           '{}'.format(assembly_graph_folder + sample + '.png')]
            subprocess.run(cmd_bandage)
        else:
            warnings.warn('No assembly for {}!'.format(sample))

    @staticmethod
    def assemble_shasta_parallel(sample_dict, output_folder, min_size, cpu, parallel):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((sample, path, output_folder, min_size, int(cpu / parallel))
                    for sample, path in sample_dict.items())
            for results in executor.map(lambda x: Methods.assemble_shasta(*x), args):
                pass

    @staticmethod
    def shasta_assembly_stats(assembly_folder, output_folder):
        # Output file
        output_stats_file = output_folder + '/shasta_stats.tsv'

        # Pandas data frame to save values
        df = pd.DataFrame(columns=['Sample', 'ReadNumber', 'TotalBases', 'ReadsN50', 'AssemblyLength', 'Contigs'])

        # Find log file(s) and parse info of interest
        log_list = glob(assembly_folder + '/*_shasta_stdout.txt', recursive=False)

        for log_file in log_list:
            sample = os.path.basename(log_file).replace("_shasta_stdout.txt", '')
            with open(log_file, 'r') as f:
                read_cnt = 0
                read_len = 0
                n50 = 0
                assembly_len = 0
                contigs = 0

                for line in f:
                    line = line.rstrip()

                    if 'Total number of reads is' in line:
                        read_cnt = line.split()[-1].replace('.', '')
                    elif 'Total number of raw bases is' in line:
                        read_len = line.split()[-1].replace('.', '')
                    elif 'N50 for read length is' in line:
                        n50 = line.split()[-2]
                    elif 'Total length of assembled sequence is' in line:
                        assembly_len = line.split()[-1]
                    elif 'The assembly graph has' in line:
                        contigs = line.split()[-3]

                data_dict = {'Sample': [sample],
                             'ReadNumber': [read_cnt],
                             'TotalBases': [read_len],
                             'ReadsN50': [n50],
                             'AssemblyLength': [assembly_len],
                             'Contigs': [contigs]}

                df = pd.concat([df, pd.DataFrame.from_dict(data_dict)], axis='index', ignore_index=True)

        # Convert df to tsv file
        df.to_csv(output_stats_file, sep='\t', index=False)

        # Remove log files
        for log_file in log_list:
            os.remove(log_file)

    @staticmethod
    def assemble_rebaler(ref, reads, sample, output_folder, cpu):
        print('\t{}'.format(sample))
        cmd = ['rebaler',
               '--threads', str(cpu),
               ref,
               reads]

        assemblies_folder = output_folder + 'all_assemblies/'
        Methods.make_folder(assemblies_folder)
        output_assembly = assemblies_folder + sample + '.fasta'

        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

        # Change sequence title from reference to sample
        with open(output_assembly, 'w') as f:
            seq_counter = 0
            for line in io.TextIOWrapper(p.stdout, encoding="utf-8"):
                if line.startswith('>'):
                    seq_counter += 1
                    line = '>{}_{}\n'.format(sample, seq_counter)
                f.write(line)

    @staticmethod
    def assemble_rebaler_parallel(ref, sample_dict, output_folder, cpu, parallel):
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            args = ((ref, path, sample, output_folder, int(cpu / parallel))
                    for sample, path in sample_dict.items())
            for results in executor.map(lambda x: Methods.assemble_rebaler(*x), args):
                pass

    # @staticmethod
    # def run_ragtag(ref, assembly, out_folder, cpu):
    #     cmd = ['ragtag.py', 'scaffold',
    #            '-t', str(cpu),  # threads
    #            '-w',  # overwrite
    #            '-o', out_folder,  # output folder
    #            ref,  # reference
    #            assembly]  # query
    #
    #     subprocess.run(cmd)

    @staticmethod
    def run_parsnp(assembly_list, output_folder, ref, cpu):
        Methods.make_folder(output_folder)

        cmd = ['parsnp',
               '-r', ref,
               '-p', str(cpu),
               '-o', output_folder,
               '-d'] + assembly_list
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def run_ksnp3(assembly_list, output_folder, ref, cpu):
        Methods.make_folder(output_folder)

        # Add ref to assembly list
        assembly_list.append(ref)

        ass_file_ksnp3 = output_folder + 'assembly.list'
        with open(ass_file_ksnp3, 'w') as f:
            for ass in assembly_list:
                ass_name = '.'.join(os.path.basename(ass).split('.')[:-1])
                f.write('{}\t{}\n'.format(ass, ass_name))

        # Combine all fasta for
        os.chdir(output_folder)
        print('\tRunning MakeFasta')
        cmd_makefasta = ['MakeFasta',
                         ass_file_ksnp3,
                         'all.fasta']
        p = subprocess.Popen(cmd_makefasta, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
        stdout, sterr = p.communicate()
        if b'% of the median size of' in stdout:
            raise Exception('{} Cannot run kSNP3.'.format(stdout.split(b'\n')[1].decode('utf-8')))

        # Get optimal k value
        print('\tRunning Kchooser')
        cmd_kchooser = ['Kchooser',
                        'all.fasta']

        p = subprocess.Popen(cmd_kchooser, stdout=subprocess.PIPE)
        optimal_k = 19
        for line in io.TextIOWrapper(p.stdout, encoding="utf-8"):
            if 'The optimum value of K is' in line:
                optimal_k = line.split()[-1].replace('.', '')

        # Remove combined fasta file
        os.remove('all.fasta')

        print('\tRunning kSNP3')
        cmd_ksnp3 = ['kSNP3',
                     '-k', optimal_k,
                     '-in', ass_file_ksnp3,
                     '-outdir', output_folder,
                     '-core',
                     '-ML',
                     '-CPU', str(cpu),
                     '-NJ',
                     '-vcf']

        subprocess.run(cmd_ksnp3, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        # Remove file list
        os.remove(ass_file_ksnp3)

    @staticmethod
    def convert_xmfa_to_fastq(xmfa_input, fasta_output):
        cmd = ['harvesttools',
               '-x', xmfa_input,
               '-M', fasta_output]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

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
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def make_tree_fasttree(aligned_fasta, output_folder):
        cmd = ['FastTree',
               '-nt',
               '-gtr',
               aligned_fasta]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
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

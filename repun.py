#!/usr/bin/env python

# BSD 3-Clause License
#
# Copyright 2023 The University of Hong Kong, Department of Computer Science
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import os
import sys
import argparse
import shlex
import subprocess
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

from collections import defaultdict, namedtuple
from argparse import SUPPRESS
try:
    from packaging.version import parse as version_parse
except ModuleNotFoundError:
    from distutils.version import LooseVersion as version_parse

from time import time

import shared.param as param
from shared.interval_tree import bed_tree_from
from shared.utils import file_path_from, folder_path_from, subprocess_popen, str2bool, legal_range_from, log_error, log_warning

major_contigs = {"chr" + str(a) for a in list(range(1, 23)) + ["X", "Y"]}.union(
    {str(a) for a in list(range(1, 23)) + ["X", "Y"]})
major_contigs_order = ["chr" + str(a) for a in list(range(1, 23)) + ["X", "Y"]] + [str(a) for a in
                                                                                   list(range(1, 23)) + ["X", "Y"]]

file_directory = os.path.dirname(os.path.realpath(__file__))
main_entry = os.path.join(file_directory, "ru.py")
MAX_STEP = 10

OutputPath = namedtuple('OutputPath', [
    'log_path',
    'tmp_file_path',
    'split_bed_path',
    'candidates_path',
    'var_path',
    'phased_vcf_path',
    'phased_bam_path',
    'vcf_output_path',
])


class Tee(object):
    def __init__(self, name, mode):
        self.file = open(name, mode)
        self.stdout = sys.stdout
        sys.stdout = self

    def __del__(self):
        sys.stdout = self.stdout
        self.file.close()

    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)

    def flush(self):
        self.file.flush()

def logging(str):
    if args.tee is None:
        print(str)
    else:
        args.tee.stdin.write(bytes(str + '\n', encoding='utf8'))

def create_output_folder(args):
    # create temp file folder
    args.output_dir = folder_path_from(args.output_dir, create_not_found=True)
    log_path = folder_path_from(os.path.join(args.output_dir, 'logs'), create_not_found=True)
    tmp_file_path = folder_path_from(os.path.join(args.output_dir, 'tmp'), create_not_found=True)
    split_bed_path = folder_path_from(os.path.join(tmp_file_path, 'split_beds'), create_not_found=True)
    candidates_path = folder_path_from(os.path.join(tmp_file_path, 'candidates'), create_not_found=True)
    var_path = folder_path_from(os.path.join(tmp_file_path, 'var'), create_not_found=True)
    phased_vcf_path = folder_path_from(os.path.join(tmp_file_path, 'phased_vcf'), create_not_found=True)
    phased_bam_path = folder_path_from(os.path.join(tmp_file_path, 'phased_bam'), create_not_found=True)
    vcf_output_path = folder_path_from(os.path.join(tmp_file_path, 'vcf_output'), create_not_found=True)
    vcf_output_path = folder_path_from(os.path.join(tmp_file_path, 'af_distribution'), create_not_found=True)

    output_path = OutputPath(log_path=log_path,
                             tmp_file_path=tmp_file_path,
                             split_bed_path=split_bed_path,
                             candidates_path=candidates_path,
                             var_path=var_path,
                             phased_vcf_path=phased_vcf_path,
                             phased_bam_path=phased_bam_path,
                             vcf_output_path=vcf_output_path)
    return output_path

def check_version(tool, pos=None, is_pypy=False):
    try:
        if is_pypy:
            proc = subprocess.run("{} -c 'import sys; print (sys.version)'".format(tool), stdout=subprocess.PIPE,
                                  shell=True)
        else:
            proc = subprocess.run([tool, "--version"], stdout=subprocess.PIPE)
        if proc.returncode != 0:
            return None
        first_line = proc.stdout.decode().split("\n", 1)[0]
        version = first_line.split()[pos]
        version = version_parse(version)
    except Exception:
        return None

    return version

def check_skip_steps_legal(args):
    skip_steps = args.skip_steps
    skip_steps_list = skip_steps.rstrip().split(",")
    if len(skip_steps_list) == 0:
        sys.exit(log_error("[ERROR] --skip_steps option provided but no skip steps index found"))
    for step in skip_steps_list:
        if int(step) < 1 or int(step) > MAX_STEP:
            sys.exit(log_error("[ERROR] --skip_steps option provided but contains invalid skip steps index, should be 1-index"))

def check_python_path():
    python_path = subprocess.run("which python", stdout=subprocess.PIPE, shell=True).stdout.decode().rstrip()
    sys.exit(log_error("[ERROR] Current python execution path: {}".format(python_path)))

def check_python_version(python):
    python_path = subprocess.run("{} --version".format(python), stdout=subprocess.PIPE, shell=True).stdout.decode().rstrip()
    return python_path.split(' ')[1]


def check_tools_version(args):

    # Core required tools (always needed)
    required_tool_version = {
        'python': version_parse('3.9.0'),
        'samtools': version_parse('1.10'),
    }

    tool_version = {
        'python': version_parse(check_python_version(args.python)),
        'samtools': check_version(tool=args.samtools, pos=1),
    }

    # whatshap is needed for germline mode with phasing enabled
    if not getattr(args, 'somatic_mode', False) and not getattr(args, 'disable_phasing', False):
        required_tool_version['whatshap'] = version_parse('1.0')
        tool_version['whatshap'] = check_version(tool=args.whatshap, pos=1)

    # pypy3 is needed internally by somatic unify_representation_somatic
    if getattr(args, 'somatic_mode', False):
        required_tool_version['pypy'] = version_parse('3.6')
        tool_version['pypy'] = check_version(tool=args.pypy, pos=0, is_pypy=True)

    for tool, version in tool_version.items():
        required_version = required_tool_version[tool]
        if version is None:
            logging(log_error("[ERROR] {} not found, please check if you are in the conda virtual environment".format(tool)))
            check_python_path()
        elif version < required_version:
            logging(log_error("[ERROR] Tool version not match, please check if you are in the conda virtual environment"))
            logging(' '.join([str(item).ljust(10) for item in ["Tool", "Version", "Required"]]))
            error_info = ' '.join([str(item).ljust(10) for item in [tool, version, '>=' + str(required_version)]])
            logging(error_info)
            check_python_path()
    return


def check_contig_in_bam(bam_fn, sorted_contig_list, samtools, allow_none=False):
    if allow_none and bam_fn is None:
        return sorted_contig_list, True

    contig_with_read_support_set = set()
    is_cram = bam_fn.lower().endswith('.cram')

    if is_cram:
        # For CRAM files, samtools idxstats is slow as it reads the whole file.
        # We fallback to checking the header which is much faster.
        possible_contigs = set()
        process = subprocess_popen(shlex.split("{} view -H {}".format(samtools, bam_fn)))
        for row in process.stdout:
            if row.startswith('@SQ'):
                columns = row.rstrip().split('\t')
                for col in columns:
                    if col.startswith('SN:'):
                        contig_name = col[3:]
                        if contig_name in sorted_contig_list:
                            possible_contigs.add(contig_name)
                            contig_with_read_support_set.add(contig_name)
        process.stdout.close()
        process.wait()

        # Check for existence of reads using random access (fast for CRAM/BAM with index)
        # We only need to know if count > 0, so reading one read is enough.
        # contigs_to_check = [c for c in sorted_contig_list if c in possible_contigs]
        # for contig_name in contigs_to_check:
        #      # check if any read exists in the region
        #      cmd = "{} view - {} {} ".format(samtools, bam_fn, contig_name)
        #      proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        #      # Try to read one byte/line. If we get anything, the contig is not empty.
        #      has_reads = False
        #      try:
        #          if proc.stdout.read(1):
        #             has_reads = True
        #      finally:
        #          proc.kill()
        #          proc.wait() # avoid zombie processes
             
        #      if has_reads:
        #          contig_with_read_support_set.add(contig_name)

    else:
        bai_process = subprocess_popen(shlex.split("{} idxstats {}".format(samtools, bam_fn)))
        for row_id, row in enumerate(bai_process.stdout):
            row = row.split('\t')
            if len(row) != 4:
                continue
            contig_name, contig_length, mapped_reads, unmapped_reads = row
            if contig_name not in sorted_contig_list:
                continue
            if int(mapped_reads) > 0:
                contig_with_read_support_set.add(contig_name)

    for contig_name in sorted_contig_list:
        if contig_name not in contig_with_read_support_set:
            logging(log_warning(
                "[WARNING] Contig name {} provided but no mapped reads found in BAM/CRAM, skip!".format(contig_name)))
    
    filtered_sorted_contig_list = [item for item in sorted_contig_list if item in contig_with_read_support_set]

    found_contig = True
    if len(filtered_sorted_contig_list) == 0:
        found_contig = False
        logging(log_warning(
            "[WARNING] No mapped reads found in BAM/CRAM for provided contigs set {}".format(' '.join(sorted_contig_list))))

    return filtered_sorted_contig_list, found_contig


def check_threads(args):
    threads = args.threads
    #sched_getaffinity is not exist in pypy
    try:
        sched_getaffinity_list = list(os.sched_getaffinity(0))
        num_cpus = len(sched_getaffinity_list)
    except:
        num_cpus = int(subprocess.run(args.python + " -c \"import os; print(len(os.sched_getaffinity(0)))\"", \
                                      stdout=subprocess.PIPE, shell=True).stdout.decode().rstrip())

    if threads > num_cpus:
        logging(log_warning(
            '[WARNING] Threads setting {} is larger than the number of available threads {} in the system,'.format(
                threads, num_cpus)))
        logging(log_warning('Set --threads={} for better parallelism.'.format(num_cpus)))
        args.threads = num_cpus
    return args

def split_extend_vcf(truth_vcf_fn, output_fn):
    expand_region_size = param.no_of_positions
    output_ctg_dict = defaultdict(list)
    unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (truth_vcf_fn)))

    for row_id, row in enumerate(unzip_process.stdout):
        if row[0] == '#':
            continue
        columns = row.strip().split(maxsplit=3)
        ctg_name = columns[0]

        center_pos = int(columns[1])
        ctg_start, ctg_end = center_pos - 1, center_pos
        if ctg_start < 0:
            sys.exit(
                log_error("[ERROR] Invalid VCF input at the {}-th row {} {}".format(row_id + 1, ctg_name, center_pos)))
        if ctg_start - expand_region_size < 0:
            continue
        expand_ctg_start = ctg_start - expand_region_size
        expand_ctg_end = ctg_end + expand_region_size

        output_ctg_dict[ctg_name].append(
            ' '.join([ctg_name, str(expand_ctg_start), str(expand_ctg_end)]))

    for key, value in output_ctg_dict.items():
        ctg_output_fn = os.path.join(output_fn, key)
        with open(ctg_output_fn, 'w') as output_file:
            output_file.write('\n'.join(value))

    unzip_process.stdout.close()
    unzip_process.wait()

    know_vcf_contig_set = set(list(output_ctg_dict.keys()))

    return know_vcf_contig_set


def split_extend_bed(bed_fn, output_fn, contig_set=None):
    expand_region_size = param.no_of_positions
    output_ctg_dict = defaultdict(list)
    unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (bed_fn)))
    for row_id, row in enumerate(unzip_process.stdout):
        if row[0] == '#':
            continue
        columns = row.strip().split()
        ctg_name = columns[0]
        if contig_set and ctg_name not in contig_set:
            continue

        ctg_start, ctg_end = int(columns[1]), int(columns[2])

        if ctg_end < ctg_start or ctg_start < 0 or ctg_end < 0:
            sys.exit(log_error(
                "[ERROR] Invalid BED input at the {}-th row {} {} {}".format(row_id + 1, ctg_name, ctg_start, ctg_end)))
        expand_ctg_start = max(0, ctg_start - expand_region_size)
        expand_ctg_end = max(0, ctg_end + expand_region_size)
        output_ctg_dict[ctg_name].append(
            ' '.join([ctg_name, str(expand_ctg_start), str(expand_ctg_end)]))

    for key, value in output_ctg_dict.items():
        ctg_output_fn = os.path.join(output_fn, key)
        with open(ctg_output_fn, 'w') as output_file:
            output_file.write('\n'.join(value))

    unzip_process.stdout.close()
    unzip_process.wait()

def write_region_bed(region):

    try:
        ctg_name, start_end = region.split(':')
        ctg_start, ctg_end = int(start_end.split('-')[0]) - 1, int(start_end.split('-')[1]) - 1  # bed format
    except:
        sys.exit("[ERROR] Please use the correct format for --region: ctg_name:start-end, your input is {}".format(
            region))
    if ctg_end < ctg_start or ctg_start < 0 or ctg_end < 0:
        sys.exit("[ERROR] Invalid region input: {}".format(region))

    output_bed_path = os.path.join(args.output_dir, 'tmp', 'region.bed')
    with open(output_bed_path, 'w') as f:
        f.write('\t'.join([ctg_name, str(ctg_start), str(ctg_end)]) + '\n')
    return output_bed_path

def check_contigs_intersection(args, fai_fn):

    MIN_CHUNK_LENGTH = 200000
    MAX_CHUNK_LENGTH = 20000000
    is_include_all_contigs = args.include_all_ctgs
    is_bed_file_provided = args.bed_fn is not None or args.region is not None
    is_known_vcf_file_provided = args.truth_vcf_fn is not None
    is_ctg_name_list_provided = args.ctg_name is not None

    if args.region is not None:
        args.bed_fn = write_region_bed(args.region)

    split_bed_path = os.path.join(args.output_dir, 'tmp', 'split_beds')
    tree = bed_tree_from(bed_file_path=args.bed_fn)
    know_vcf_contig_set = split_extend_vcf(truth_vcf_fn=args.truth_vcf_fn, output_fn=split_bed_path) if is_known_vcf_file_provided else set()
    contig_set = set(args.ctg_name.split(',')) if is_ctg_name_list_provided else set()

    if not args.include_all_ctgs:
        logging("[INFO] --include_all_ctgs not enabled, use chr{1..22,X,Y} and {1..22,X,Y} by default")
    else:
        logging("[INFO] --include_all_ctgs enabled")

    if is_ctg_name_list_provided and is_bed_file_provided:
        logging(log_warning("[WARNING] both --ctg_name and --bed_fn provided, will only proceed with the contigs appeared in both"))

    if is_ctg_name_list_provided and is_known_vcf_file_provided:
        logging(log_warning("[WARNING] both --ctg_name and --truth_vcf_fn provided, will only proceed with the contigs appeared in both"))

    if is_ctg_name_list_provided:
        contig_set = contig_set.intersection(
            set(tree.keys())) if is_bed_file_provided else contig_set
        contig_set = contig_set.intersection(
            know_vcf_contig_set) if is_known_vcf_file_provided else contig_set
    else:
        contig_set = contig_set.union(
            set(tree.keys())) if is_bed_file_provided else contig_set

        contig_set = contig_set.union(
            know_vcf_contig_set) if is_known_vcf_file_provided else contig_set

    # if each split region is too small(long) for given default chunk num, will increase(decrease) the total chunk num
    default_chunk_num = 0
    DEFAULT_CHUNK_SIZE = args.chunk_size
    contig_length_list = []
    contig_chunk_num = {}

    with open(fai_fn, 'r') as fai_fp:
        for row in fai_fp:
            columns = row.strip().split("\t")
            contig_name, contig_length = columns[0], int(columns[1])
            if not is_include_all_contigs and (
            not (is_bed_file_provided or is_ctg_name_list_provided or is_known_vcf_file_provided)) and str(
                    contig_name) not in major_contigs:
                continue

            if is_bed_file_provided and contig_name not in tree:
                continue
            if is_ctg_name_list_provided and contig_name not in contig_set:
                continue
            if is_known_vcf_file_provided and contig_name not in contig_set:
                continue

            contig_set.add(contig_name)
            contig_length_list.append(contig_length)
            chunk_num = int(
                contig_length / float(DEFAULT_CHUNK_SIZE)) + 1 if contig_length % DEFAULT_CHUNK_SIZE else int(
                contig_length / float(DEFAULT_CHUNK_SIZE))
            contig_chunk_num[contig_name] = max(chunk_num, 1)

    if default_chunk_num > 0:
        min_chunk_length = min(contig_length_list) / float(default_chunk_num)
        max_chunk_length = max(contig_length_list) / float(default_chunk_num)

    contigs_order = major_contigs_order + list(contig_set)

    sorted_contig_list = sorted(list(contig_set), key=lambda x: contigs_order.index(x))

    found_contig = True
    if not len(contig_set):
        if is_bed_file_provided:
            all_contig_in_bed = ' '.join(list(tree.keys()))
            logging(log_warning("[WARNING] No contig in --bed_fn was found in the reference, contigs in BED {}: {}".format(args.bed_fn, all_contig_in_bed)))
        if is_known_vcf_file_provided:
            all_contig_in_vcf = ' '.join(list(know_vcf_contig_set))
            logging(log_warning("[WARNING] No contig in --truth_vcf_fn was found in the reference, contigs in VCF {}: {}".format(args.truth_vcf_fn, all_contig_in_vcf)))
        if is_ctg_name_list_provided:
            all_contig_in_ctg_name = ' '.join(args.ctg_name.split(','))
            logging(log_warning("[WARNING] No contig in --ctg_name was found in the reference, contigs in contigs list: {}".format(all_contig_in_ctg_name)))
        found_contig = False
    else:
        for c in sorted_contig_list:
            if c not in contig_chunk_num:
                logging(log_warning(("[WARNING] Contig {} given but not found in the reference".format(c))))
        sorted_contig_list = [c for c in sorted_contig_list if c in contig_chunk_num]

        # check contig in bam have support reads
        sorted_contig_list, found_contig = check_contig_in_bam(bam_fn=args.bam_fn, sorted_contig_list=sorted_contig_list,
                                                               samtools=args.samtools)

    if not found_contig:
        log_warning("[WARNING] Exit calling because no contig was found in BAM!")
        sys.exit(0)
    logging('[INFO] Execute RU in contigs: {}'.format(' '.join(sorted_contig_list)))
    logging('[INFO] Number of chunks for each contig: {}'.format(
        ' '.join([str(contig_chunk_num[c]) for c in sorted_contig_list])))

    if default_chunk_num > 0 and max_chunk_length > MAX_CHUNK_LENGTH:
        logging(log_warning(
            '[WARNING] The maximum chunk size set {} is larger than the suggested maximum chunk size {}, consider setting a larger --chunk_num= instead for better parallelism.'.format(
                min_chunk_length, MAX_CHUNK_LENGTH)))

    elif default_chunk_num > 0 and min_chunk_length < MIN_CHUNK_LENGTH:
        logging(log_warning(
            '[WARNING] The minimum chunk size set {} is smaller than the suggested  minimum chunk size {}, consider setting a smaller --chunk_num= instead.'.format(
                min_chunk_length, MIN_CHUNK_LENGTH)))

    if default_chunk_num == 0 and max(contig_length_list) < DEFAULT_CHUNK_SIZE / 5:
        logging(log_warning(
            '[WARNING] The length of the longest contig {} is more than five times smaller smaller than the default chunk size {}, consider setting a smaller --chunk_size= instead for better parallelism.'.format(
                max(contig_length_list), DEFAULT_CHUNK_SIZE)))

    contig_path = os.path.join(args.output_dir, 'tmp', 'CONTIGS')
    with open(contig_path, 'w') as output_file:
        output_file.write('\n'.join(sorted_contig_list))

    chunk_list = []
    chunk_list_path = os.path.join(args.output_dir, 'tmp', 'CHUNK_LIST')
    with open(chunk_list_path, 'w') as output_file:
        for contig_name in sorted_contig_list:
            chunk_num = contig_chunk_num[contig_name] if args.chunk_num is None else args.chunk_num
            for chunk_id in range(1, chunk_num + 1):
                output_file.write(contig_name + ' ' + str(chunk_id) + ' ' + str(chunk_num) + '\n')
                chunk_list.append((contig_name, chunk_id, chunk_num))
    args.chunk_list = chunk_list

    return args


def check_args(args):

    if args.conda_prefix is None:
        if 'CONDA_PREFIX' in os.environ:
            args.conda_prefix = os.environ['CONDA_PREFIX']
        else:
            try:
                python_path = subprocess.run('which python', stdout=subprocess.PIPE, shell=True).stdout.decode().rstrip()
                args.conda_prefix = os.path.dirname(os.path.dirname(python_path))
            except:
                sys.exit(log_error("[ERROR] Conda prefix not found, please activate a correct conda environment."))

    args.bam_fn = file_path_from(file_name=args.bam_fn, exit_on_not_found=True)
    bai_fn = file_path_from(file_name=args.bam_fn, suffix=".bai", exit_on_not_found=False, sep='.')
    crai_fn = file_path_from(file_name=args.bam_fn, suffix=".crai", exit_on_not_found=False, sep='.')
    csi_fn = file_path_from(file_name=args.bam_fn, suffix=".csi", exit_on_not_found=False, sep='.')

    args.ref_fn = file_path_from(file_name=args.ref_fn, exit_on_not_found=True)
    fai_fn = file_path_from(file_name=args.ref_fn, suffix=".fai", exit_on_not_found=True, sep='.')

    # args.bed_fn = file_path_from(file_name=args.bed_fn, exit_on_not_found=True)
    args.bed_fn = file_path_from(file_name=args.bed_fn, exit_on_not_found=False)

    args.truth_vcf_fn = file_path_from(file_name=args.truth_vcf_fn, exit_on_not_found=True)

    if bai_fn is None and crai_fn is None and csi_fn is None:
        sys.exit(log_error("[ERROR] BAM index file {} or {} not found. Please run `samtools index $BAM` first.".format(args.bam_fn + '.bai',
                                                                                      args.bam_fn + '.crai')))

    if args.min_af is None:
        args.min_af = param.min_af

    if args.min_coverage is None:
        args.min_coverage = 2

    if args.chunk_size is None:
        args.chunk_size = 5000000

    if args.max_af_for_somatic_unification is None:
        args.max_af_for_somatic_unification = 0.08
    if args.vaf_threshold_for_pass is None:
        args.vaf_threshold_for_pass = 0.08

    if args.snv_maximum_edit_distance is None:
        args.snv_maximum_edit_distance = 0

    if args.allow_candidate_haplotype_shorter_than_truth_haplotype is None:
        args.allow_candidate_haplotype_shorter_than_truth_haplotype = True

    if args.indel_maximum_edit_distance is None:
        args.indel_maximum_edit_distance = 4
    if args.platform not in {'ont', 'hifi', 'ilmn'}:
        logging(log_error('[ERROR] Invalid platform input, optional: {ont, hifi, ilmn}'))

    if args.skip_steps is not None:
        check_skip_steps_legal(args)

    legal_range_from(param_name="threads", x=args.threads, min_num=1, exit_out_of_range=True)
    legal_range_from(param_name="min_coverage", x=args.min_coverage, min_num=0, exit_out_of_range=True)
    legal_range_from(param_name="min_af", x=args.min_af, min_num=0, max_num=1, exit_out_of_range=True)
    legal_range_from(param_name="chunk_size", x=args.chunk_size, min_num=0, exit_out_of_range=True)

    args.output_path = create_output_folder(args)
    check_tools_version(args=args)
    args = check_contigs_intersection(args=args, fai_fn=fai_fn)

    return args


def print_args(args):
    version = '0.1.4'
    logging("")
    logging("[INFO] RU VERSION: {}".format(version))
    logging("[INFO] BAM FILE PATH: {}".format(args.bam_fn))
    logging("[INFO] REFERENCE FILE PATH: {}".format(args.ref_fn))
    logging("[INFO] TRUTH VCF FILE PATH: {}".format(args.truth_vcf_fn))
    logging("[INFO] PLATFORM: {}".format(args.platform))
    logging("[INFO] THREADS: {}".format(args.threads))
    logging("[INFO] OUTPUT FOLDER: {}".format(args.output_dir))
    logging("[INFO] OUTPUT VCF PATH: {}".format(os.path.join(args.output_dir, 'unified.vcf.gz')))
    logging("[INFO] BED FILE PATH: {}".format(args.bed_fn))
    logging("[INFO] REGION FOR EXECUTING: {}".format(args.region))
    logging("[INFO] CONTIGS FOR EXECUTING: {}".format(args.ctg_name))
    logging("[INFO] CONDA BINARY PREFIX: {}".format(args.conda_prefix))
    logging("[INFO] SAMTOOLS BINARY PATH: {}".format(args.samtools))
    logging("[INFO] WHATSHAP BINARY PATH: {}".format(args.whatshap))
    logging("[INFO] PYTHON BINARY PATH: {}".format(args.python))
    if args.somatic_mode:
        logging("[INFO] PYPY BINARY PATH: {}".format(args.pypy))
    logging("[INFO] CHUNK SIZE: {}".format(args.chunk_size))
    logging("[INFO] MINIMUM AF: {}".format(args.min_af))
    logging("[INFO] DISABLE PHASING: {}".format(args.disable_phasing))
    logging("[INFO] ENABLE DRY RUN: {}".format(args.dry_run))
    logging("[INFO] ENABLE INCLUDING ALL CTGS FOR EXECUTING: {}".format(args.include_all_ctgs))
    logging("[INFO] ENABLE REMOVING INTERMEDIATE FILES: {}".format(args.remove_intermediate_dir))
    logging("")

    if args.cmdline is not None and args.cmdline != "":
        with open(args.output_dir + '/tmp/CMD', 'w') as f:
            f.write(args.cmdline + '\n')

    return args


def print_command_line(args):

    try:

        cmdline = os.path.realpath(__file__)
        cmdline += ' --bam_fn {} '.format(args.bam_fn)
        cmdline += '--ref_fn {} '.format(args.ref_fn)
        cmdline += '--truth_vcf_fn {} '.format(args.truth_vcf_fn)
        cmdline += '--threads {} '.format(args.threads)
        cmdline += '--platform {} '.format(args.platform)
        cmdline += '--output_dir {} '.format(args.output_dir)
        cmdline += '--bed_fn {} '.format(args.bed_fn) if args.bed_fn is not None else ""
        cmdline += '--ctg_name {} '.format(args.ctg_name) if args.ctg_name is not None else ""
        cmdline += '--region {} '.format(args.region) if args.region is not None else ""
        cmdline += '--min_af {} '.format(args.snv_min_af) if args.min_af is not None else ""
        cmdline += '--min_coverage {} '.format(args.min_coverage) if args.min_coverage is not None else ""
        cmdline += '--chunk_size {} '.format(args.chunk_size) if args.chunk_size is not None else ""
        cmdline += '--sample_name {} '.format(args.sample_name) if args.sample_name != "SAMPLE" else ""
        cmdline += '--output_prefix {} '.format(args.output_prefix) if args.output_prefix != "output" else ""
        cmdline += '--remove_intermediate_dir ' if args.remove_intermediate_dir else ""
        cmdline += '--include_all_ctgs ' if args.include_all_ctgs else ""
        cmdline += '--python {} '.format(args.python) if args.python != "python3" else ""
        cmdline += '--pypy {} '.format(args.pypy) if args.pypy != "pypy3" else ""
        cmdline += '--samtools {} '.format(args.samtools) if args.samtools != "samtools" else ""
        cmdline += '--whatshap {} '.format(args.whatshap) if args.whatshap != "whatshap" else ""
        cmdline += '--disable_phasing ' if args.disable_phasing else ""
        cmdline += '--skip_steps {} '.format(args.skip_steps) if args.skip_steps is not None else ""
        cmdline += '--conda_prefix {} '.format(args.conda_prefix) if args.conda_prefix is not None else ""
        args.cmdline = cmdline
    except:
        return args
    logging("[COMMAND] " + cmdline + '\n')


# ---------------------------------------------------------------------------
# Top-level worker functions for concurrent.futures parallelization.
# These must be module-level (not nested) so they can be pickled by
# ProcessPoolExecutor.
# ---------------------------------------------------------------------------

def _run_split_extend_bed(bed_fn, output_fn, ctg_name):
    """Worker: run SplitExtendBed for a single contig."""
    import argparse as _ap
    from src.SplitExtendBed import split_extend_bed as _fn
    _fn(_ap.Namespace(
        bed_fn=bed_fn, output_fn=output_fn, ctgName=ctg_name,
        ctgStart=None, ctgEnd=None, expand_region_size=100,
    ))


def _run_get_truth(vcf_fn, var_fn, ctg_name):
    """Worker: run GetTruth for a single contig."""
    import argparse as _ap
    from src.GetTruth import OutputVariant as _fn
    _fn(_ap.Namespace(
        vcf_fn=vcf_fn, var_fn=var_fn, ctgName=ctg_name,
        ctgStart=None, ctgEnd=None, truth_vcf_fn=None,
    ))


def _run_extract_candidates(bam_fn, ref_fn, ctg_name, chunk_id, chunk_num,
                             samtools, min_af, extend_bed, unify_repre_fn,
                             platform, min_coverage):
    """Worker: run ExtractCandidates for a single chunk."""
    import argparse as _ap
    import shared.param as _param
    from src.ExtractCandidates import ExtractCandidates as _fn
    _fn(_ap.Namespace(
        platform=platform, bam_fn=bam_fn, ref_fn=ref_fn,
        tensor_can_fn='PIPE', vcf_fn=None, min_af=min_af,
        ctgName=ctg_name, ctgStart=None, ctgEnd=None,
        bed_fn=None, gvcf=False, sampleName='SAMPLE',
        samtools=samtools, minCoverage=min_coverage,
        minMQ=5, minBQ=0, max_depth=144,
        phasing_info_in_bam=True,
        phasing_window_size=_param.phasing_window_size,
        extend_bed=extend_bed, indel_fn=None,
        base_err=0.001, gq_bin_size=5, bp_resolution=False,
        zstd=_param.zstd, test_pos=False,
        chunk_num=chunk_num, chunk_id=chunk_id,
        phased_vcf_fn=None, add_no_phasing_data_training=False,
        unify_repre=True, unify_repre_fn=unify_repre_fn,
        full_aln_regions=None, need_phasing=False,
        need_realignment=False, min_bq_cut=0,
    ))


def _run_unify_representation(vcf_fn, candidates_fn, ref_fn, output_vcf_fn,
                               platform, ctg_name, bed_fn, chunk_id, chunk_num,
                               sample_name):
    """Worker: run UnifyRepresentation for a single chunk."""
    import argparse as _ap
    from src.UnifyRepresentation import UnifyRepresentation as _fn
    _fn(_ap.Namespace(
        platform=platform, vcf_fn=vcf_fn, ref_fn=ref_fn,
        candidates_fn=candidates_fn, output_vcf_fn=output_vcf_fn,
        sampleName=sample_name, ctgName=ctg_name,
        ctgStart=None, ctgEnd=None,
        max_candidates_distance=100, max_calculate_count=10000,
        partition_size=15, minimum_allele_gap=0.15,
        bed_fn=bed_fn, subsample_ratio=1000,
        test_pos=False, chunk_num=chunk_num, chunk_id=chunk_id,
    ))


def _run_cal_af_distribution(tumor_bam_fn, truth_vcf_fn, output_path,
                              ctg_name, chunk_id, chunk_num, samtools):
    """Worker: run cal_af_distribution for a single chunk."""
    import argparse as _ap
    from src.cal_af_distribution import cal_af as _fn
    _fn(_ap.Namespace(
        tumor_bam_fn=tumor_bam_fn, normal_bam_fn=None,
        input_vcf_fn=None, bed_fn=None,
        truth_vcf_fn=truth_vcf_fn,
        input_filter_tag=None, truth_filter_tag=None,
        ctg_name=ctg_name, output_path=output_path,
        threads=1, phase_output=False, samtools=samtools,
        min_mq=5, min_bq=0, min_bq_cut=0,
        chunk_num=chunk_num, chunk_id=chunk_id, debug=False,
    ))


def _run_get_truth_somatic(vcf_fn, var_fn, output_bed_fn, ctg_name, max_af):
    """Worker: run get_truth_somatic for a single contig."""
    import argparse as _ap
    from src.get_truth_somatic import OutputVariant as _fn
    _fn(_ap.Namespace(
        vcf_fn=vcf_fn, var_fn=var_fn,
        output_bed_fn=output_bed_fn, ctgName=ctg_name,
        ctgStart=None, ctgEnd=None,
        max_af=max_af, truth_vcf_fn=None,
    ))


def _run_unify_representation_somatic(bam_fn, ref_fn, var_fn, output_vcf_fn,
                                       platform, ctg_name, bed_fn, extend_bed,
                                       chunk_id, chunk_num, sample_name,
                                       snv_maximum_edit_distance,
                                       indel_maximum_edit_distance,
                                       allow_shorter, samtools, pypy):
    """Worker: run unify_representation_somatic for a single chunk."""
    import argparse as _ap
    from src.unify_representation_somatic import UnifyRepresentation as _fn
    _fn(_ap.Namespace(
        platform=platform, var_fn=var_fn, ref_fn=ref_fn,
        candidate_details_fn=None,
        output_vcf_fn=output_vcf_fn, output_info_fn=None,
        sampleName=sample_name, ctgName=ctg_name,
        ctgStart=None, ctgEnd=None,
        find_match_with_most_read_support=None,
        max_candidates_distance=100, max_calculate_count=10000,
        partition_size=15, minimum_allele_gap=0.15,
        bed_fn=bed_fn, vcf_fn=None, extend_bed=extend_bed,
        subsample_ratio=1000, test_pos=0,
        chunk_num=chunk_num, chunk_id=chunk_id,
        snv_maximum_edit_distance=snv_maximum_edit_distance,
        indel_maximum_edit_distance=indel_maximum_edit_distance,
        allow_candidate_haplotype_shorter_than_truth_haplotype=allow_shorter,
        minimum_supporting_reads=None, minimum_supporting_af=None,
        min_af=None, bam_fn=bam_fn, pypy=pypy, samtools=samtools,
    ))


def _run_subprocess_task(cmd):
    """Worker: run a single shell command, raising CalledProcessError on failure."""
    subprocess.run(cmd, shell=True, check=True)


def _run_whatshap_haplotag_and_index(ctg, output_dir, ref_fn, bam_fn,
                                      whatshap, samtools, sample_name):
    """Worker: run whatshap haplotag + samtools index for a single contig."""
    haplotag_cmd = (
        '{} haplotag --output {}/tmp/phased_bam/{}_{}.bam '
        '--reference {} --regions {} --ignore-read-groups '
        '{}/tmp/phased_vcf/phased_{}.vcf.gz {}'
    ).format(whatshap, output_dir, sample_name, ctg,
             ref_fn, ctg, output_dir, ctg, bam_fn)
    subprocess.run(haplotag_cmd, shell=True, check=True)

    index_cmd = '{} index {}/tmp/phased_bam/{}_{}.bam'.format(
        samtools, output_dir, sample_name, ctg)
    subprocess.run(index_cmd, shell=True, check=True)


def _get_sorted_contig_list(chunk_list):
    """Return unique contig names from chunk_list in their original order."""
    seen = set()
    contigs = []
    for ctg, _cid, _cnum in chunk_list:
        if ctg not in seen:
            seen.add(ctg)
            contigs.append(ctg)
    return contigs


def _run_parallel(tasks, worker_fn, max_workers, use_processes=True):
    """
    Execute *tasks* in parallel using ProcessPoolExecutor (default) or
    ThreadPoolExecutor.  Each element of *tasks* is a tuple of positional
    arguments forwarded to *worker_fn*.  Re-raises the first exception
    encountered in any worker.
    """
    if not tasks:
        return
    executor_cls = ProcessPoolExecutor if use_processes else ThreadPoolExecutor
    with executor_cls(max_workers=max_workers) as executor:
        futures = [executor.submit(worker_fn, *task) for task in tasks]
        errors = []
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as exc:
                errors.append(exc)
    if errors:
        raise errors[0]


# ---------------------------------------------------------------------------
# Workflow functions
# ---------------------------------------------------------------------------


def ru_somatic(args):
    """Somatic-mode workflow using Python-internal parallelization."""

    output_dir = args.output_dir
    sorted_contig_list = _get_sorted_contig_list(args.chunk_list)
    skip_steps = args.skip_steps.rstrip().split(',') if args.skip_steps else None
    step = 1

    # STEP 1: Split bed file regions (parallel over contigs)
    echo = "[INFO] STEP {}: Split bed file regions according to contig name".format(step)
    logging(echo)
    logging("[INFO] Running SplitExtendBed for {} contigs with {} workers".format(
        len(sorted_contig_list), args.threads))
    logging("")
    if not args.dry_run:
        if skip_steps is not None and str(step) in skip_steps:
            logging("[INFO] --skip_steps is enabled, skip running step {}.".format(step))
        else:
            seb_tasks = [
                (args.bed_fn,
                 os.path.join(output_dir, 'tmp', 'split_beds', ctg),
                 ctg)
                for ctg in sorted_contig_list
            ]
            try:
                _run_parallel(seb_tasks, _run_split_extend_bed, args.threads)
            except Exception as exc:
                sys.stderr.write("ERROR in STEP {}: {}\n".format(step, exc))
                sys.exit(1)
    logging("")
    step += 1

    # STEP 2: Calculate AF from BAM (parallel over chunks)
    echo = "[INFO] STEP {}: Calculate AF from truth".format(step)
    logging(echo)
    logging("[INFO] Running cal_af_distribution for {} chunks with {} workers".format(
        len(args.chunk_list), args.threads))
    logging("")
    if not args.dry_run:
        if skip_steps is not None and str(step) in skip_steps:
            logging("[INFO] --skip_steps is enabled, skip running step {}.".format(step))
        else:
            caf_tasks = [
                (args.bam_fn, args.truth_vcf_fn,
                 os.path.join(output_dir, 'tmp', 'af_distribution',
                              'truth_{}_{}'.format(ctg, chunk_id)),
                 ctg, chunk_id, chunk_num, args.samtools)
                for ctg, chunk_id, chunk_num in args.chunk_list
            ]
            try:
                _run_parallel(caf_tasks, _run_cal_af_distribution, args.threads)
                # Concatenate per-chunk AF files into a single truths file
                truths_path = os.path.join(output_dir, 'tmp', 'af_distribution', 'truths')
                import glob as _glob
                with open(truths_path, 'w') as _out:
                    for _f in sorted(_glob.glob(
                            os.path.join(output_dir, 'tmp', 'af_distribution', 'truth_*'))):
                        with open(_f) as _fh:
                            _out.write(_fh.read())
            except Exception as exc:
                sys.stderr.write("ERROR in STEP {}: {}\n".format(step, exc))
                sys.exit(1)
    logging("")
    step += 1

    # STEP 3: Get true variant label information (parallel over contigs)
    echo = "[INFO] STEP {}: Get true variant label information from VCF file".format(step)
    logging(echo)
    logging("[INFO] Running get_truth_somatic for {} contigs with {} workers".format(
        len(sorted_contig_list), args.threads))
    logging("")
    if not args.dry_run:
        if skip_steps is not None and str(step) in skip_steps:
            logging("[INFO] --skip_steps is enabled, skip running step {}.".format(step))
        else:
            truths_path = os.path.join(output_dir, 'tmp', 'af_distribution', 'truths')
            gtv_tasks = [
                (truths_path,
                 os.path.join(output_dir, 'tmp', 'var', 'var_{}'.format(ctg)),
                 os.path.join(output_dir, 'tmp', 'split_beds', 'bed_{}'.format(ctg)),
                 ctg, args.max_af_for_somatic_unification)
                for ctg in sorted_contig_list
            ]
            try:
                _run_parallel(gtv_tasks, _run_get_truth_somatic, args.threads)
            except Exception as exc:
                sys.stderr.write("ERROR in STEP {}: {}\n".format(step, exc))
                sys.exit(1)
    logging("")
    step += 1

    # STEP 4: Unify somatic variants (parallel over chunks)
    echo = "[INFO] STEP {}: Unify somatic variant".format(step)
    logging(echo)
    logging("[INFO] Running unify_representation_somatic for {} chunks with {} workers".format(
        len(args.chunk_list), args.threads))
    logging("")
    if not args.dry_run:
        if skip_steps is not None and str(step) in skip_steps:
            logging("[INFO] --skip_steps is enabled, skip running step {}.".format(step))
        else:
            lv_tasks = [
                (args.bam_fn, args.ref_fn,
                 os.path.join(output_dir, 'tmp', 'var', 'var_{}'.format(ctg)),
                 os.path.join(output_dir, 'tmp', 'vcf_output',
                              'vcf_{}_{}.vcf'.format(ctg, chunk_id)),
                 args.platform, ctg,
                 os.path.join(output_dir, 'tmp', 'split_beds', 'bed_{}'.format(ctg)),
                 os.path.join(output_dir, 'tmp', 'split_beds', 'bed_{}'.format(ctg)),
                 chunk_id, chunk_num, args.sample_name,
                 args.snv_maximum_edit_distance,
                 args.indel_maximum_edit_distance,
                 args.allow_candidate_haplotype_shorter_than_truth_haplotype,
                 args.samtools, args.pypy)
                for ctg, chunk_id, chunk_num in args.chunk_list
            ]
            try:
                _run_parallel(lv_tasks, _run_unify_representation_somatic, args.threads)
            except Exception as exc:
                sys.stderr.write("ERROR in STEP {}: {}\n".format(step, exc))
                sys.exit(1)
    logging("")
    step += 1

    # STEP 5: Merge and Sort VCFs
    echo = "[INFO] STEP {}: Merge and Sort VCFs".format(step)
    logging(echo)
    logging("[INFO] Running sort_vcf_somatic")
    logging("")
    if not args.dry_run:
        if skip_steps is not None and str(step) in skip_steps:
            logging("[INFO] --skip_steps is enabled, skip running step {}.".format(step))
        else:
            import argparse as _ap
            from src.sort_vcf_somatic import sort_vcf_from as _sort_somatic
            sort_args = _ap.Namespace(
                output_fn=os.path.join(output_dir, 'unified_truths.vcf'),
                output_candidate_fn=os.path.join(output_dir, 'unified_candidates.vcf'),
                truth_vcf_fn=os.path.join(
                    output_dir, 'tmp', 'af_distribution', 'truths'),
                input_dir=os.path.join(output_dir, 'tmp', 'vcf_output'),
                vcf_fn_prefix='vcf_',
                vcf_fn_suffix='.vcf',
                ref_fn=args.ref_fn,
                sampleName=args.sample_name,
                contigs_fn=os.path.join(output_dir, 'tmp', 'CONTIGS'),
                cmd_fn=os.path.join(output_dir, 'tmp', 'CMD'),
                vaf_threshold_for_pass=args.vaf_threshold_for_pass,
            )
            try:
                _sort_somatic(sort_args)
            except Exception as exc:
                sys.stderr.write("ERROR in STEP {}: {}\n".format(step, exc))
                sys.exit(1)
    logging("")

    if args.remove_intermediate_dir:
        logging("[INFO] Removing intermediate files in {}/tmp ...".format(output_dir))
        subprocess.run('rm -rf {}/tmp'.format(output_dir), shell=True)


def ru(args):
    """Germline-mode workflow using Python-internal parallelization."""

    output_dir = args.output_dir
    sorted_contig_list = _get_sorted_contig_list(args.chunk_list)
    skip_steps = args.skip_steps.rstrip().split(',') if args.skip_steps else None
    step = 1

    if not args.disable_phasing:
        # STEP 1: Whatshap UnPhasing (serial, single command)
        echo = "[INFO] STEP {}: Whatshap UnPhasing".format(step)
        logging(echo)
        wu_cmd = '{} unphase {} > {}'.format(
            args.whatshap, args.truth_vcf_fn,
            os.path.join(output_dir, 'tmp', 'INPUT.vcf.gz'))
        logging("[INFO] RUN: " + wu_cmd)
        logging("")
        if not args.dry_run:
            if skip_steps is not None and str(step) in skip_steps:
                logging("[INFO] --skip_steps is enabled, skip running step {}.".format(step))
            else:
                try:
                    subprocess.check_call(wu_cmd, shell=True)
                except subprocess.CalledProcessError as exc:
                    sys.stderr.write("ERROR in STEP {}: {}\n".format(step, exc))
                    sys.exit(1)
        logging("")
        step += 1

        # STEP 2: Whatshap Phasing (parallel over contigs)
        echo = "[INFO] STEP {}: Whatshap Phasing INPUT VCF file".format(step)
        logging(echo)
        logging("[INFO] Phasing {} contigs with {} workers".format(
            len(sorted_contig_list), args.threads))
        logging("")
        if not args.dry_run:
            if skip_steps is not None and str(step) in skip_steps:
                logging("[INFO] --skip_steps is enabled, skip running step {}.".format(step))
            else:
                phase_tasks = [
                    (('{} phase --output {}/tmp/phased_vcf/phased_{}.vcf.gz '
                      '--reference {} --chromosome {} --ignore-read-groups '
                      '--distrust-genotypes {}/tmp/INPUT.vcf.gz {}').format(
                          args.whatshap, output_dir, ctg, args.ref_fn, ctg,
                          output_dir, args.bam_fn),)
                    for ctg in sorted_contig_list
                ]
                try:
                    _run_parallel(phase_tasks,
                                  _run_subprocess_task, args.threads,
                                  use_processes=False)
                except Exception as exc:
                    sys.stderr.write("ERROR in STEP {}: {}\n".format(step, exc))
                    sys.exit(1)
        logging("")
        step += 1

        # STEP 3: Index phased VCF files (parallel over contigs)
        echo = "[INFO] STEP {}: Index VCF file".format(step)
        logging(echo)
        logging("[INFO] Indexing {} phased VCF files".format(len(sorted_contig_list)))
        logging("")
        if not args.dry_run:
            if skip_steps is not None and str(step) in skip_steps:
                logging("[INFO] --skip_steps is enabled, skip running step {}.".format(step))
            else:
                tabix_tasks = [
                    ('tabix -f -p vcf {}/tmp/phased_vcf/phased_{}.vcf.gz'.format(
                        output_dir, ctg),)
                    for ctg in sorted_contig_list
                ]
                try:
                    _run_parallel(tabix_tasks, _run_subprocess_task,
                                  args.threads, use_processes=False)
                except Exception as exc:
                    sys.stderr.write("ERROR in STEP {}: {}\n".format(step, exc))
                    sys.exit(1)
        logging("")
        step += 1

        # STEP 4: Whatshap Haplotagging + BAM indexing (parallel over contigs)
        echo = "[INFO] STEP {}: Whatshap Haplotagging BAM file".format(step)
        logging(echo)
        logging("[INFO] Haplotagging {} contigs with {} workers".format(
            len(sorted_contig_list), args.threads))
        logging("")
        if not args.dry_run:
            if skip_steps is not None and str(step) in skip_steps:
                logging("[INFO] --skip_steps is enabled, skip running step {}.".format(step))
            else:
                haplotag_tasks = [
                    (ctg, output_dir, args.ref_fn, args.bam_fn,
                     args.whatshap, args.samtools, args.sample_name)
                    for ctg in sorted_contig_list
                ]
                try:
                    _run_parallel(haplotag_tasks, _run_whatshap_haplotag_and_index,
                                  args.threads, use_processes=False)
                except Exception as exc:
                    sys.stderr.write("ERROR in STEP {}: {}\n".format(step, exc))
                    sys.exit(1)
        logging("")
        step += 1

    # STEP: Split bed file regions (parallel over contigs)
    echo = "[INFO] STEP {}: Split bed file regions according to contig name".format(step)
    logging(echo)
    logging("[INFO] Running SplitExtendBed for {} contigs with {} workers".format(
        len(sorted_contig_list), args.threads))
    logging("")
    if not args.dry_run:
        if skip_steps is not None and str(step) in skip_steps:
            logging("[INFO] --skip_steps is enabled, skip running step {}.".format(step))
        else:
            seb_tasks = [
                (args.bed_fn,
                 os.path.join(output_dir, 'tmp', 'split_beds', ctg),
                 ctg)
                for ctg in sorted_contig_list
            ]
            try:
                _run_parallel(seb_tasks, _run_split_extend_bed, args.threads)
            except Exception as exc:
                sys.stderr.write("ERROR in STEP {}: {}\n".format(step, exc))
                sys.exit(1)
    logging("")
    step += 1

    # STEP: Get true variant labels (parallel over contigs)
    echo = "[INFO] STEP {}: Get true variant label information from VCF file".format(step)
    logging(echo)
    logging("[INFO] Running GetTruth for {} contigs with {} workers".format(
        len(sorted_contig_list), args.threads))
    logging("")
    if not args.dry_run:
        if skip_steps is not None and str(step) in skip_steps:
            logging("[INFO] --skip_steps is enabled, skip running step {}.".format(step))
        else:
            gtv_tasks = [
                (args.truth_vcf_fn,
                 os.path.join(output_dir, 'tmp', 'var', 'var_{}'.format(ctg)),
                 ctg)
                for ctg in sorted_contig_list
            ]
            try:
                _run_parallel(gtv_tasks, _run_get_truth, args.threads)
            except Exception as exc:
                sys.stderr.write("ERROR in STEP {}: {}\n".format(step, exc))
                sys.exit(1)
    logging("")
    step += 1

    # STEP: Extract variant candidates from BAM (parallel over chunks)
    echo = "[INFO] STEP {}: Extract variant candidates from BAM file".format(step)
    logging(echo)
    logging("[INFO] Running ExtractCandidates for {} chunks with {} workers".format(
        len(args.chunk_list), args.threads))
    logging("")
    if not args.dry_run:
        if skip_steps is not None and str(step) in skip_steps:
            logging("[INFO] --skip_steps is enabled, skip running step {}.".format(step))
        else:
            ec_tasks = [
                (os.path.join(output_dir, 'tmp', 'phased_bam',
                              '{}_{}.bam'.format(args.sample_name, ctg))
                 if not args.disable_phasing else args.bam_fn,
                 args.ref_fn, ctg, chunk_id, chunk_num,
                 args.samtools, args.min_af,
                 os.path.join(output_dir, 'tmp', 'split_beds', ctg),
                 os.path.join(output_dir, 'tmp', 'candidates',
                              '{}_{}_{}'.format(args.sample_name, chunk_id, ctg)),
                 args.platform, args.min_coverage)
                for ctg, chunk_id, chunk_num in args.chunk_list
            ]
            try:
                _run_parallel(ec_tasks, _run_extract_candidates, args.threads)
            except Exception as exc:
                sys.stderr.write("ERROR in STEP {}: {}\n".format(step, exc))
                sys.exit(1)
    logging("")
    step += 1

    # STEP: Unify variant representation (parallel over chunks)
    echo = "[INFO] STEP {}: Label variant".format(step)
    logging(echo)
    logging("[INFO] Running UnifyRepresentation for {} chunks with {} workers".format(
        len(args.chunk_list), args.threads))
    logging("")
    if not args.dry_run:
        if skip_steps is not None and str(step) in skip_steps:
            logging("[INFO] --skip_steps is enabled, skip running step {}.".format(step))
        else:
            lv_tasks = [
                (os.path.join(output_dir, 'tmp', 'var', 'var_{}'.format(ctg)),
                 os.path.join(output_dir, 'tmp', 'candidates',
                              '{}_{}_{}'.format(args.sample_name, chunk_id, ctg)),
                 args.ref_fn,
                 os.path.join(output_dir, 'tmp', 'vcf_output',
                              'vcf_{}_{}'.format(ctg, chunk_id)),
                 args.platform, ctg,
                 os.path.join(output_dir, 'tmp', 'split_beds', ctg),
                 chunk_id, chunk_num, args.sample_name)
                for ctg, chunk_id, chunk_num in args.chunk_list
            ]
            try:
                _run_parallel(lv_tasks, _run_unify_representation, args.threads)
            except Exception as exc:
                sys.stderr.write("ERROR in STEP {}: {}\n".format(step, exc))
                sys.exit(1)
    logging("")
    step += 1

    # STEP: Sort, compress, and index VCF
    echo = "[INFO] STEP {}: Sort VCF file".format(step)
    logging(echo)
    logging("[INFO] Running SortVcf and compressing output")
    logging("")
    if not args.dry_run:
        if skip_steps is not None and str(step) in skip_steps:
            logging("[INFO] --skip_steps is enabled, skip running step {}.".format(step))
        else:
            import argparse as _ap
            from src.SortVcf import sort_vcf_from as _sort_vcf
            sort_args = _ap.Namespace(
                output_fn=os.path.join(output_dir, 'unified.vcf'),
                input_dir=os.path.join(output_dir, 'tmp', 'vcf_output'),
                vcf_fn_prefix='vcf_',
                vcf_fn_suffix='',
                ref_fn=args.ref_fn,
                sampleName=args.sample_name,
                contigs_fn=os.path.join(output_dir, 'tmp', 'CONTIGS'),
                cmd_fn=os.path.join(output_dir, 'tmp', 'CMD'),
            )
            try:
                _sort_vcf(sort_args)
            except Exception as exc:
                sys.stderr.write("ERROR in STEP {}: {}\n".format(step, exc))
                sys.exit(1)
    logging("")

    if args.remove_intermediate_dir:
        logging("[INFO] Removing intermediate files in {}/tmp ...".format(output_dir))
        subprocess.run('rm -rf {}/tmp'.format(output_dir), shell=True)




def ru_parser():

    parser = argparse.ArgumentParser(
        description="Run Clair3 Representation Unification. Example run: repun -b BAM_FN -r REF_FN --truth_vcf_fn TRUTH_VCF_FN -o OUTPUT_DIR -t THREADS -p PLATFORM")

    # print version
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {}'.format('0.0.1'))

    required_params = parser.add_argument_group('Required parameters')
    required_params.add_argument(
        '-b',
        "--bam_fn",
        type=str,
        required=True,
        default=None,
        help="BAM file input. The input file must be samtools indexed."
    )

    required_params.add_argument(
        "-r",
        "--ref_fn",
        type=str,
        required=True,
        default=None,
        help="FASTA reference file input. The input file must be samtools indexed."
    )

    required_params.add_argument(
        "--truth_vcf_fn",
        type=str,
        required=True,
        default=None,
        help="Truth VCF file input."
    )

    required_params.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=True,
        default=None,
        help="Output directory."
    )

    required_params.add_argument(
        "-t",
        "--threads",
        required=True,
        type=int,
        default=None,
        help="Max #threads to be used."
    )

    required_params.add_argument(
        "-p",
        "--platform",
        required=True,
        type=str,
        default=None,
        help="Select the sequencing platform of the input. Possible options: {ont, hifi, ilmn}."
    )

    optional_params = parser.add_argument_group('Optional parameters')
    optional_params.add_argument(
        "-c",
        "--ctg_name",
        type=str,
        default=None,
        help="The name of the contigs to be processed. Split by ',' for multiple contigs. Default: all contigs will be processed."
    )

    optional_params.add_argument(
        "--bed_fn",
        type=str,
        default=None,
        help="Path to a BED file. Execute RU only in the provided BED regions."
    )

    optional_params.add_argument(
        "--region",
        type=str,
        default=None,
        help="A region to be processed. Format: `ctg_name:start-end` (start is 1-based)."
    )

    optional_params.add_argument(
        "--min_af",
        type=float,
        default=None,
        help="Minimal AF required for a variant to be called. Default: 0.08."
    )

    optional_params.add_argument(
        "--min_coverage",
        type=int,
        default=None,
        help="Minimal coverage required for a variant to be called. Default: 4."
    )

    optional_params.add_argument(
        "-s",
        "--sample_name",
        type=str,
        default="SAMPLE",
        help="Define the sample name to be shown in the VCF file. Default: SAMPLE."
    )

    optional_params.add_argument(
        "--output_prefix",
        type=str,
        default="output",
        help="Prefix for output VCF filename. Default: output."
    )

    optional_params.add_argument(
        "--remove_intermediate_dir",
        action='store_true',
        help="Remove intermediate directory before finishing to save disk space."
    )

    optional_params.add_argument(
        "--include_all_ctgs",
        action='store_true',
        help="Execute RU on all contigs, otherwise call in chr{1..22,X,Y} and {1..22,X,Y}."
    )

    optional_params.add_argument(
        '-d',
        "--dry_run",
        action='store_true',
        help="Print the commands that will be ran."
    )

    optional_params.add_argument(
        "--python",
        type=str,
        default="python3",
        help="Absolute path of python, python3 >= 3.9 is required."
    )

    optional_params.add_argument(
        "--pypy",
        type=str,
        default="pypy3",
        help="Absolute path of pypy3 (>= 3.6). Used internally for somatic mode."
    )

    optional_params.add_argument(
        "--samtools",
        type=str,
        default="samtools",
        help="Absolute path of samtools, samtools version >= 1.10 is required."
    )

    optional_params.add_argument(
        "--whatshap",
        type=str,
        default="whatshap",
        help="Absolute path of whatshap, whatshap >= 1.0 is required."
    )

    optional_params.add_argument(
        "--parallel",
        type=str,
        default="parallel",
        help=SUPPRESS,
    )

    optional_params.add_argument(
        "--disable_phasing",
        action='store_true',
        help="Disable phasing with whatshap. Work with the haplotagged BAM to speed up unification"
    )

    optional_params.add_argument(
        "--chunk_size",
        type=int,
        default=None,
        help=SUPPRESS
    )

    optional_params.add_argument(
        "--chunk_num",
        type=int,
        default=None,
        help=SUPPRESS
    )

    optional_params.add_argument(
        "--output_path",
        type=str,
        default=None,
        help=SUPPRESS
    )

    optional_params.add_argument(
        "--skip_steps",
        type=str,
        default=None,
        help=SUPPRESS
    )

    optional_params.add_argument(
        "--tee",
        type=str,
        default=None,
        help=SUPPRESS
    )

    optional_params.add_argument(
        "--conda_prefix",
        type=str,
        default=None,
        help=SUPPRESS
    )

    optional_params.add_argument(
        "--cmdline",
        type=str,
        default=None,
        help=SUPPRESS
    )

    optional_params.add_argument(
        "--somatic_mode",
        action='store_true',
        help="Enable somatic mode. Default: False."
    )

    somatic_params = parser.add_argument_group('Somatic parameters')

    somatic_params.add_argument(
        "--max_af_for_somatic_unification",
        type=float,
        default=None,
        help='Maximum allelic fraction for a somatic variant to be unified. Default: 0.08'
    )

    somatic_params.add_argument(
        "--vaf_threshold_for_pass",
        type=float,
        default=None,
        help="If set, variants with >VAF will be marked as PASS, or LowVAF otherwise. Default: 0.08"
    )

    somatic_params.add_argument(
        '--snv_maximum_edit_distance',
        type=int,
        default=None,
        help="Maximum SNV edit distance that allow to be unified, default: 0")

    somatic_params.add_argument(
        '--indel_maximum_edit_distance',
        type=int,
        default=None,
        help="Maximum Indel edit distance that allow to be unified, default: 4")

    somatic_params.add_argument(
        '--allow_candidate_haplotype_shorter_than_truth_haplotype',
        type=str2bool,
        default=None,
        help="Allow the candidate haplotype shorter than truth haploytpe, default: True")
    return parser

def main():
    """
    Main interface for Clair3 RU.
    """

    global args

    call_start_time = time()

    parser = ru_parser()
    args = parser.parse_args()

    args.output_dir = folder_path_from(args.output_dir, create_not_found=True)
    tee_logger = os.path.join(args.output_dir, 'ru.log' if not args.dry_run else "ru_dry_run.log")
    if os.path.exists(tee_logger):
        subprocess.run("mv {} {}".format(tee_logger, tee_logger + '.bak'), shell=True)
    try:
        args.tee = subprocess.Popen(['tee', tee_logger], stdin=subprocess.PIPE, bufsize=0)
    except:
        logging(log_warning("[WARNING] `tee` not found, disable `tee` logging!"))
        args.tee = None

    logging("")

    args = print_command_line(args)
    args = check_args(args)
    args = print_args(args)
    if args.somatic_mode:
        ru_somatic(args)
    else:
        ru(args)

    runtime = time() - call_start_time
    logging("[INFO] Total time elapsed: %im%.2fs\n" % (int(runtime/60), int(runtime % 60)))
    if args.somatic_mode:
        logging("[INFO] Finish Representation Unification, output vcf file using truth coordinate: {}/unified_truths.vcf.gz\n".format(args.output_dir))
        logging("[INFO] Finish Representation Unification, output vcf file using candidate coordinate: {}/unified_candidates.vcf.gz\n".format(args.output_dir))
    else:
        logging("[INFO] Finish Representation Unification, output vcf file: {}/unified.vcf.gz\n".format(args.output_dir))
    if args.tee is not None:
        args.tee.stdin.close()

if __name__ == '__main__':
    main()

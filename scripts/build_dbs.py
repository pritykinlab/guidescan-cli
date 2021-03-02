'''
Builds a Guidescan database from scratch with
default parameters.

Works by recursively calling the script from LSF 
job completion hooks with varying state parameters 
to specify the part of the state machine the process 
is in. This script essentially initializes a 
distributed state machine to perform its work.

Author: 

  Henri Schmidt
  henri.schmidt@tufts.edu
'''

import os
import sys
import argparse
import random
import subprocess as sp
import fcntl 
import time
import re

from pathlib import Path

def get_bsub_args(state, bsub_files, name, threads=1):
    mem = STATE_PARAMS[state]['mem']
    time = STATE_PARAMS[state]['time']
    return list(map(str, [
        '-W', time, 
        '-R', f'rusage[mem={mem}]',
        '-R', 'span[hosts=1]',
        '-n', threads,
        '-J', name,
        '-o', f'{bsub_files}/%J.stdout',
        '-eo', f'{bsub_files}/%J.stderr'
    ]))

'''
Unparses command line arguments to list.
'''
def unparse_to_list(args):
    args = vars(args)
    xs = [args['organism']]
    for flag, item in args.items():
        if flag == 'organism':
            continue
        xs += [f'--{flag}', item]
    return xs

def make_cmd(bsub_args, sp_args):
    args = ['bsub']
    args += bsub_args
    args.append(f'python {sys.argv[0]} {" ".join(map(str, sp_args))}')
    return list(map(str, args))

def log_state(state_file, message):
    with open(state_file, 'a') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.write(message)

def start_bsub(old_state, state_file, bsub_args, sp_args):
    cmd_args = make_cmd(bsub_args, sp_args)

    try:
        sp.run(cmd_args, check=True)
    except sp.CalledProcessError as e:
        log_state(state_file, f'FAILED\t{old_state}\n')
        sys.exit(1)

'''
Sets up the distributed system by spawning the 
watcher as a standalone job.
'''
def initial_state(args):
    if args.job_id == None:
        args.job_id = random.randint(0, 1_000_000)

    if args.state_file == None:
        args.state_file = f"state_file_{args.job_id}.txt"
    
    log_state(args.state_file, f'STARTED\tinitial\n')
    args.state = 'watcher'

    sp_args = unparse_to_list(args)
    bsub_args = get_bsub_args(args.state, args.bsub_files, f'watcher-{args.job_id}')

    Path(args.bsub_files).mkdir(parents=True, exist_ok=True)
    Path(args.results).mkdir(parents=True, exist_ok=True)

    start_bsub('initial', args.state_file, bsub_args, sp_args)
    log_state(args.state_file, f'STARTED\twatcher\n')
    log_state(args.state_file, f'COMPLETED\tinitial\n')
    sys.exit(0)
    
'''
Watches and manages the state of the distributed 
system. Responsible for spawning new jobs once the 
required jobs have been completed. 
'''
def watcher_state(args):
    def get_status():
        with open(args.state_file, 'r') as f:
            states = list(map(lambda l: l.split(), f.readlines()))
            started_states = map(lambda s: s[1], filter(lambda s: s[0] == 'STARTED', states))
            failed_states = map(lambda s: s[1], filter(lambda s: s[0] == 'FAILED', states))
            completed_states = map(lambda s: s[1], filter(lambda s: s[0] == 'COMPLETED', states))
            return list(started_states), list(failed_states), list(completed_states)
        
    while True:
        started, failed, completed = get_status()

        if failed:
            log_state(args.state_file, f'FAILED\twatcher\n')
            sys.exit(1)
            
        if 'initial' in completed and 'gen-kmers' not in started:
            args.state = 'gen-kmers'
            sp_args = unparse_to_list(args)
            bsub_args = get_bsub_args(args.state, args.bsub_files, 
                                      f'gen-kmers-{args.job_id}' )
            start_bsub('watcher', args.state_file, bsub_args, sp_args)
            log_state(args.state_file, f'STARTED\tgen-kmers\n')
           
        if 'initial' in completed and 'gen-idx' not in started:
            args.state = 'gen-idx'
            sp_args = unparse_to_list(args)
            bsub_args = get_bsub_args(args.state, args.bsub_files,
                                      f'gen-idx-{args.job_id}' )
            start_bsub('watcher', args.state_file, bsub_args, sp_args)
            log_state(args.state_file, f'STARTED\tgen-idx\n')

        if 'gen-kmers' in completed and 'build-dbs' not in started:
            args.state = 'build-dbs'
            sp_args = unparse_to_list(args)
            bsub_args = get_bsub_args(args.state, args.bsub_files,
                                      f'build-dbs-{args.job_id}' )
            start_bsub('watcher', args.state_file, bsub_args, sp_args)
            log_state(args.state_file, f'STARTED\tbuild-dbs\n')

        if 'build-dbs' in completed and 'merge-dbs' not in started:
            args.state = 'merge-dbs'
            sp_args = unparse_to_list(args)
            bsub_args = get_bsub_args(args.state, args.bsub_files,
                                      f'merge-dbs-{args.job_id}' )
            start_bsub('watcher', args.state_file, bsub_args, sp_args)
            log_state(args.state_file, f'STARTED\tmerge-dbs\n')

        if 'merge-dbs' in completed:
            break

        time.sleep(0.25)

    log_state(args.state_file, f'COMPLETED\twatcher\n')
    sys.exit(0)

def generate_kmers_state(args):
    kmer_dir = f'{args.results}/kmers'
    kmer_file = f'{kmer_dir}/all_kmers.txt'
    kmer_file_shuffled = f'{kmer_dir}/all_kmers_shuffled.txt'

    kmer_prefix = f'{kmer_dir}/kmer_'
    kmer_suffix = '.txt'
    
    Path(kmer_dir).mkdir(parents=True, exist_ok=True)

    gs_args = [
        'guidescan', 'kmers', 
        '-o', kmer_file,
        args.organism
    ]
    
    try:
        gs = sp.run(gs_args, check=True)
        with open(kmer_file_shuffled, 'w') as f:
            sp.run(['shuf', kmer_file], check=True, stdout=f) 
        sp.run(['split', '-n', f'l/{NUM_PARTITIONS}', '-d', '-a', 
                str(6),  '--additional-suffix', kmer_suffix,
                kmer_file_shuffled, kmer_prefix], check=True)
    except sp.CalledProcessError as e:
        log_state(state_file, f'FAILED\tgen-kmers\n')
        sys.exit(1)

    log_state(args.state_file, f'COMPLETED\tgen-kmers\n')

def generate_index_state(args):
    tmp_kmers = f'tmp_kmers_{args.job_id}.txt'
    tmp_out = f'tmp_kmers_{args.job_id}.txt'

    tmp_fd = open(tmp_kmers, 'w')
    tmp_fd.close()

    gs_args = [
        'guidescan', 'build', 
        '-f', tmp_kmers, 
        '-n', '1', 
        '-o', tmp_out,
        args.organism
    ]
    
    try:
        gs = sp.run(gs_args, check=True)
    except sp.CalledProcessError as e:
        log_state(state_file, f'FAILED\tgen-idx\n')
        sys.exit(1)
    
    os.remove(tmp_out)
    log_state(args.state_file, f'COMPLETED\tgen-idx\n')

def build_dbs_state(args):
    if args.data == 'built-dbs':
        log_state(args.state_file, f'COMPLETED\tbuild-dbs\n')
        return

    kmer_dir = f'{args.results}/kmers'
    db_dir = f'{args.results}/split_dbs/'
    
    Path(db_dir).mkdir(parents=True, exist_ok=True)
    
    jobs = []
    for path in os.listdir(kmer_dir):
        res = re.search('kmer_(?P<id>\\d+)\\.txt', path)
        if res is None:
            continue
    
        kmer_split_id = int(res.group('id'))
        kmer_job_name = f'build-split-db-{kmer_split_id}'
        bsub_args = get_bsub_args('build-split-db', args.bsub_files, kmer_job_name)
        sp_args = ['guidescan', 'build', '-n', '1', '-f', f'{kmer_dir}/{path}', 
                   '-o', f'{db_dir}/{kmer_split_id}.sam', args.organism]

        try:
            cmd_args = ['bsub'] + bsub_args + [' '.join(sp_args)]
            sp.run(cmd_args, check=True)
        except sp.CalledProcessError as e:
            log_state(state_file, f'FAILED\tbuild-dbs\n')
            sys.exit(1)

        jobs.append(kmer_job_name)

    bsub_args = get_bsub_args('build-dbs', args.bsub_files, f'built-dbs-{args.job_id}')
    bsub_args += ['-w', ' && '.join([f'done({job})' for job in jobs])]

    args.state = 'build-dbs'
    args.data = 'built-dbs' 
    sp_args = unparse_to_list(args)
    start_bsub('build-dbs', args.state_file, bsub_args, sp_args)

def merge_dbs_state(args):
    guide_db = open(f'{args.results}/guide_db.sam', 'w')

    db_dir = f'{args.results}/split_dbs/'

    sam_files = os.scandir(db_dir)
    sam_db = next(sam_files)

    with open(sam_db, 'r') as sf:
        for line in sf:
            guide_db.write(line)
    
    for sam_db in sam_files:
        with open(sam_db, 'r') as sf:
            for line in sf:
                if line.startswith('@'): continue
                guide_db.write(line)

    guide_db.close()
    log_state(args.state_file, f'COMPLETED\tmerge-dbs\n')

STATE_MACHINE = {
    'initial': initial_state,
    'watcher': watcher_state,
    'gen-kmers': generate_kmers_state,
    'gen-idx': generate_index_state,
    'build-dbs': build_dbs_state,
    'merge-dbs': merge_dbs_state,
}

STATE_PARAMS = {
    'initial':        {'mem': '512MB', 'time': '1:00'},
    'watcher':        {'mem': '512MB', 'time': '8:00'},
    'gen-kmers':      {'mem': '2GB', 'time': '8:00'},
    'gen-idx':        {'mem': '12GB', 'time': '8:00'},
    'build-dbs':      {'mem': '512MB', 'time': '8:00'},
    'build-split-db': {'mem': '12GB', 'time': '8:00'},
    'merge-dbs':      {'mem': '4GB', 'time': '8:00'},
}

NUM_PARTITIONS=2000

def parse_arguments():
    parser = argparse.ArgumentParser(description='Build a Guidescan database from scratch')
    parser.add_argument('--state', help='state of database construction',
                        choices=['gen-idx', 'gen-kmers', 'shuf-kmers',
                                 'split-kmers', 'build-dbs', 'merge-dbs',
                                 'initial', 'watcher'],
                        default='initial')
    parser.add_argument('organism', help='location of FASTA file for organism')
    parser.add_argument('--state_file', help='stores the current state of the system',
                        default=None)
    parser.add_argument('--job_id', help='the identity of the system',
                        default=None)
    parser.add_argument('--bsub_files', help='location of LSF bsub file directory',
                        default='bsub_files')
    parser.add_argument('--results', help='location of result file directory',
                        default='results')
    parser.add_argument('--data', help='Arbitrary data to be passed between states.')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_arguments()
    STATE_MACHINE[args.state](args)

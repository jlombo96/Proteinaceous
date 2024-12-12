import os
import subprocess as sp
import multiprocessing as mp
from glob import glob as glob

def run_makeblastdb(db=None):
    """
    Arguments:
    """
    if db is None:
        raise ValueError("ERROR: No database file provided")
    if not os.path.exists(db):
        raise FileNotFoundError(f"ERROR: Database file {db} not found")
    cmd = f"makeblastdb -in {db} -dbtype prot"
    sp.run(cmd,shell=True)
    return cmd

def run_blastp(db=None,query=None,output_path=None,params=None,evalue=1E-4,threads=None):
    """
    """
    if db is None:
        raise ValueError("ERROR: No database file provided")
    if query is None:
        raise ValueError("ERROR: No query file provided")
    if output_path is None:
        raise ValueError("ERROR: No output file provided")
    if not os.path.exists(db):
        raise FileNotFoundError(f"ERROR: Database file {db} not found")
    if not os.path.exists(query):
        raise FileNotFoundError(f"ERROR: Query file {query} not found")
    if threads is None:
        threads = 1
    if params is None:
        params = ""
    cmd = f"blastp -db {db} -query {query} -evalue {evalue} -out {output_path}.blastp -num_threads {threads} {params}"
    print(cmd)
    sp.run(cmd,shell=True)
    return cmd

def run_blastp_wrapper(job):
    """
    """
    return run_blastp(**job)

def run_hmmpress(db=None):
    """
    Arguments:
    """
    if db is None:
        raise ValueError("ERROR: No HMM file provided")
    if not os.path.exists(db):
        raise FileNotFoundError(f"ERROR: HMM file {db} notfound")
    cmd = f"hmmpress {db}"
    extensions = ['.h3f','.h3i','.h3m','.h3p']
    index_files = [db + ext for ext in extensions]
    for file in index_files:
        if os.path.exists(file):
            print(f"{file} index already exists, removing.")
            os.remove(file)
    sp.run(cmd,shell=True)
    return cmd


def run_hmmsearch(db=None,query=None,output_path=None,params=None,threads=None,evalue=1e-4):
    """
    # hmmsearch :: search profile(s) against a sequence database
    # HMMER 3.4 (Aug 2023); http://hmmer.org/
    # Copyright (C) 2023 Howard Hughes Medical Institute.
    # Freely distributed under the BSD open source license.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Usage: hmmsearch [options] <hmmfile> <seqdb>

    Basic options:
    -h : show brief help on version and usage

    Options directing output:
    -o <f>           : direct output to file <f>, not stdout
    -A <f>           : save multiple alignment of all hits to file <f>
    --tblout <f>     : save parseable table of per-queryuence hits to file <f>
    --domtblout <f>  : save parseable table of per-domain hits to file <f>
    --pfamtblout <f> : save table of hits and domains to file, in Pfam format <f>
    --acc            : prefer accessions over names in output
    --noali          : don't output alignments, so output is smaller
    --notextw        : unlimit ASCII text output line width
    --textw <n>      : set max width of ASCII text output lines  [120]  (n>=120)

    Options controlling reporting thresholds:
    -E <x>     : report sequences <= this E-value threshold in output  [10.0]  (x>0)
    -T <x>     : report sequences >= this score threshold in output
    --domE <x> : report domains <= this E-value threshold in output  [10.0]  (x>0)
    --domT <x> : report domains >= this score cutoff in output

    Options controlling inclusion (significance) thresholds:
    --incE <x>    : consider sequences <= this E-value threshold as significant
    --incT <x>    : consider sequences >= this score threshold as significant
    --incdomE <x> : consider domains <= this E-value threshold as significant
    --incdomT <x> : consider domains >= this score threshold as significant

    Options controlling model-specific thresholding:
    --cut_ga : use profile's GA gathering cutoffs to set all thresholding
    --cut_nc : use profile's NC noise cutoffs to set all thresholding
    --cut_tc : use profile's TC trusted cutoffs to set all thresholding

    Options controlling acceleration heuristics:
    --max    : Turn all heuristic filters off (less speed, more power)
    --F1 <x> : Stage 1 (MSV) threshold: promote hits w/ P <= F1  [0.02]
    --F2 <x> : Stage 2 (Vit) threshold: promote hits w/ P <= F2  [1e-3]
    --F3 <x> : Stage 3 (Fwd) threshold: promote hits w/ P <= F3  [1e-5]
    --nobias : turn off composition bias filter

    Other expert options:
    --nonull2     : turn off biased composition score corrections
    -Z <x>        : set # of comparisons done, for E-value calculation
    --domZ <x>    : set # of significant seqs, for domain E-value calculation
    --seed <n>    : set RNG seed to <n> (if 0: one-time arbitrary seed)  [42]
    --tformat <s> : assert target <seqfile> is in format <s>: no autodetection
    --cpu <n>     : number of parallel CPU workers to use for multithreads  [2]
    """
    print(db)
    print(query)

    reserved_params = ' -o,--tblout','--domtblout','--pfamtblout', '--cpu'
    if params is None:
        params = ""
    if query is None:
        raise ValueError("No sequence file provided")
    if db is None:
        raise ValueError("No HMM database provided")
    if threads is None:
        threads = 1
    for reserved_param in reserved_params:
        if reserved_param in params:
            raise ValueError(f"Cannot use reserved parameter {reserved_param} in params")
    if evalue is not None:
        params += f" -E {evalue}"
    if output_path is None:
        raise ValueError("No output path provided")
    output_file = output_path + ".out"
    tblout = output_path + ".tblout"
    domtblout = output_path + ".domtblout"
    pfamtblout = output_path + ".pfamtblout"
    cmd = f"hmmsearch --tblout {tblout} --domtblout {domtblout} --pfamtblout {pfamtblout} -o {output_file} {params} {db} {query}"
    #print(cmd)
    hmmsearch_process = sp.run(cmd, shell=True,stderr=None,stdout=sp.PIPE)
    with open(output_path + ".log",'w') as logstream:
        logstream.write(hmmsearch_process.stdout.decode())
    logstream.close()
    return cmd

def run_hmmsearch_wrapper(job):
    """
    """
    return run_hmmsearch(**job)


def glob_protein_files(indir):
    """
    """
    extensions = ['.fa','.fasta','.faa']
    extensions = extensions + [ext + ".gz" for ext in extensions]
    filepaths = []
    for ext in extensions:
        filepaths += glob(f"{indir}/**/*{ext}",recursive=True)
    return filepaths

def glob_hmm_files(dir,ext=".tblout"):
    """
    """
    hmm_files = glob(f"{dir}/*.tblout",recursive=True)
    return hmm_files

def schedule_searches(file_paths,db_path,outprefix=None,cpus=1,params=None,method="hmmsearch",overwrite=False):
    """
    """
    # print the arguments of this function
    job_dictionary_list = []
    
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"ERROR: HMM database {db_path} not found")
    if params is None:
        params = ""
    if outprefix is None:
        output_prefix = "proteinsearch"
    for file in file_paths:
        if not os.path.exists(file):
            print(f"ERROR: File {file} not found")
        output_filename = f"{output_prefix}_{os.path.dirname(file).split('/')[-1]}_{os.path.basename(file).split('.')[0]}"
        output_path = f"{os.path.dirname(file)}/{output_filename}"
        if method == "hmmsearch":
            required_files = [f"{output_path}.tblout",f"{output_path}.domtblout",f"{output_path}.pfamtblout"]
        elif method == "blastp":
            required_files = [f"{output_path}.blastp"]
        if all([os.path.exists(f) for f in required_files]):
            if not overwrite:
                print(f"Search results already exist for {file}, skipping")
                continue
        job_dictionary = {"query":file,"db":db_path,"output_path":output_path,"params":params}
        job_dictionary_list.append(job_dictionary)
    return job_dictionary_list

def run_searches(indir=None,outdir=None,db_path=None,cpus=1,dry_run=False,search_method="hmmsearch",params=None,overwrite=False):
    """
    """
    if indir is None:
        raise ValueError("ERROR: No input directory provided")
    if db_path is None:
        raise ValueError("ERROR: No HMM database provided")
    if outdir is None:
        outdir = "hmmsearch_results"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if search_method == "hmmsearch":
        build_func = run_hmmpress
        search_func = run_hmmsearch_wrapper
    elif search_method == "blastp":
        build_func = run_makeblastdb
        search_func = run_blastp_wrapper
    else:
        raise ValueError("ERROR: Invalid search method")
    protein_files = glob_protein_files(indir)
    job_list = schedule_searches(protein_files,db_path,cpus=cpus,params=params,method=search_method,overwrite=overwrite)
    if not dry_run:
        #run_hmmpress(hmm=db_path)
        build_func(db=db_path)
        pool = mp.Pool(cpus)
        results = pool.map(search_func,job_list)
        pool.close()
        pool.join()
        #result_files = glob_hmm_files(outdir)
    else:
        print("Dry run")
        result_files = []
        for job in job_list:
            print(job)
    return job_list

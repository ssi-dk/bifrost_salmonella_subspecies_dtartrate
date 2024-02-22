

import subprocess
import getposfromsam
import traceback
from collections import Counter

def run_cmd(command, log):
    with open(log.out_file, "a+") as out, open(log.err_file, "a+") as err:
        command_log_out, command_log_err = subprocess.Popen(command, shell=True).communicate()
        if command_log_err == None:
            command_log_err = ""
        if command_log_out == None:
            command_log_out = ""
        out.write(command_log_out)
        err.write(command_log_err)

def rule__dtartrate(input: object, output: object, dtartrate_db: str, log: object) -> None:
    try:
        bwa = subprocess.Popen(
            ["bwa", "mem", dtartrate_db, input[0], input[1]], 
            stdout = subprocess.PIPE, 
            stderr = subprocess.PIPE).communicate()
        elprep = subprocess.Popen(
            ["elprep","filter", "/dev/stdin", "/dev/stdout", "--filter-unmapped-reads", "--sorting-order", "coordinate", "--nr-of-threads", "1"],
            stdin=subprocess.PIPE, 
            stdout = subprocess.PIPE).communicate(input=bwa[0])
        d_tartrate_bases = Counter(getposfromsam.get_pos_from_sam(10,11,elprep[0].decode()))
        with open(output._file, 'w') as fh:
            fh.write(str(d_tartrate_bases))
    except Exception:
        with open(log.err_file, "w+") as fh:
            fh.write(traceback.format_exc())

rule__dtartrate(
    snakemake.input,
    snakemake.output,
    snakemake.params.dtartrate_db,
    snakemake.log)



import subprocess
import NearestSubspecies
import traceback

def run_cmd(command, log):
    with open(log.out_file, "a+") as out, open(log.err_file, "a+") as err:
        command_log_out, command_log_err = subprocess.Popen(command, shell=True).communicate()
        if command_log_err == None:
            command_log_err = ""
        if command_log_out == None:
            command_log_out = ""
        out.write(command_log_out)
        err.write(command_log_err)

def rule__subspecies(input: object, output: object, params: object, log: object) -> None:
    try:
        ST = params.mlsttype
        mlstDir = params.mlst_db
        subspecies_ref = params.subspecies_reference
        try:
            int(ST)
            query = NearestSubspecies.query_from_ST(ST, mlstDir)
            result = NearestSubspecies.subspecies_from_query(query, subspecies_ref)
            subspecies = result[0][1] # Result is a sorted table of (score, subspecies) entries
        except ValueError:
            """ No known st """
            subspecies = 'NA'
        with open(output._file, 'w') as fh:
            fh.write(str(subspecies))
    except Exception:
        with open(log.err_file, "w+") as fh:
            fh.write(traceback.format_exc())

rule__subspecies(
    snakemake.input,
    snakemake.output,
    snakemake.params,
    snakemake.log)

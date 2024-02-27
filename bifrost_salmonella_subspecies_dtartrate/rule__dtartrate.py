
import subprocess
import sys
import traceback
from collections import Counter

from bifrost_salmonella_subspecies_dtartrate import getposfromsam


def rule__dtartrate(input: object, output: object, params: object, log: object) -> None:
    dtartrate_db = params.dtartrate_db
    try:
        bwa = subprocess.Popen(
            ["bwa", "mem", dtartrate_db, input[0], input[1]], 
            stdout = subprocess.PIPE, 
            stderr = subprocess.PIPE).communicate()
        elprep = subprocess.Popen(
            ["elprep","filter", "/dev/stdin", "/dev/stdout", "--filter-unmapped-reads", "--sorting-order", "coordinate", "--nr-of-threads", "1"],
            stdin=subprocess.PIPE, 
            stdout = subprocess.PIPE).communicate(input=bwa[0])
        elprep_output = elprep[0].decode().split('\n')
        d_tartrate_bases = Counter(getposfromsam.get_pos_from_sam(10,11,elprep_output))
        print(dict(d_tartrate_bases), file=sys.stderr)
        with open(output._file, 'w') as fh:
            fh.write(str(dict(d_tartrate_bases)))
    except Exception:
        with open(log.err_file, "w+") as fh:
            fh.write(traceback.format_exc())
        raise

if __name__ == "__main__":
    rule__dtartrate(
        snakemake.input,
        snakemake.output,
        snakemake.params,
        snakemake.log)

import subprocess
import logging
from multiprocessing import Pool
import argparse
import json
from typing import List

# logger
logger = logging.getLogger(__name__)
file_handler = logging.FileHandler('fetch-with-sratoolkit.log', 'a')
logger.addHandler(file_handler)
stream_handler = logging.StreamHandler()
logger.addHandler(stream_handler)

def read_sra_files(filename: str) -> List[str]:    
    try:
        with open(args['jsonfile']) as fobj:
            json_obj = json.load(fobj)
            return json_obj["sraIds"]
    except Exception as e:
        logger.error(str(e))
        return []
    
def execute_cmd(cmd: str):
    try:
        subprocess.call(cmd, shell=True)
    except Exception as e:
        logger.error("Error occurred in cmd: " + cmd + ". Error: " + str(e))

def download_sra_files(ids: List[str], output_dir: str, temp_dir: str):
    # generate commands
    cmds = []
    for id in ids:
        fasterq_dump = f"fasterq-dump -p -t {temp_dir} -O {output_dir} " + id
        cmds.append(fasterq_dump)

    with Pool(8) as p:
        result = p.map(execute_cmd, cmds)
        return result

if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('jsonfile', help="Json file specifying sra file ids to download")
    parser.add_argument(
        "-o", "--outputDir", 
        help="Output directory (optional), the default output directory is the current directory",
        default="./"
    )
    parser.add_argument(
        "-t", "--tempDir",
        help="Temporary directory (optional), the default temporary directory is the current directory",
        default="./"
    )
    args = vars(parser.parse_args())
    fn = args["jsonfile"]

    # read sra ids
    ids = read_sra_files(fn)

    # download sra files
    result = download_sra_files(ids, args["outputDir"], args["tempDir"])
    



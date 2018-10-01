''' Generic methods for methylSeq analysis.

    This module contains general methods used in sequence
    analysis pipelines
'''
import os
import sys
import re
import errno
import glob
import logging
import subprocess

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def slash_terminate(s):
    '''Add a trailing '/' to the string if not present
    '''
    return(s if s.endswith('/') else s + '/')

def runAndCheck(cmd, msg):
    '''Run the given command and check the return value.

        If the return value of the command is not 0 (success),
        print the given error message and return
        exit code 1.

        --cmd - the command to be run
        --msg - the message to print if the command returns an error
    '''
    logger = logging.getLogger("analysis.functions")
    logger.debug("Invoking: " + cmd)
    try:
        retcode = subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT)
        if retcode != 0:
            logger.debug("System call returned "+str(retcode))
            logger.error(msg+". Quitting.")
        return(retcode)
    except OSError as e:
        logger.exception("Error invoking command")
        print >>sys.stderr, "Execution exception:", e
        return(1)

def srun(cmd, ncores="1", mem="", sdout="", sderr="", exclusive=False, job_name=""):
    '''Run the given command through srun.

    If srun is available, the command will be run
    with the given cores and memory. If srun is not
    installed, the command is invoked directly.

    --cmd       - the command to be run (required)
    --ncores    - the number of cores (optional)
    --mem       - the memory in Gb (optional)
    --sdout     - the desired output file for sdout redirection (optional)
    --sderr     - the desired output file for sderr redirection (optional)
    --exclusive - if true, the exclusive flag will be set; only this job will use the node (optional)
    --job_name  - the name for the job to display in sview / squeue (optional)
    '''
    logger = logging.getLogger("analysis.functions")
    if(srun_is_installed()):
        s = "srun -c %s " % (str(ncores))
        s = s + " --mem="+str(mem)+"G " if mem else s
        s = s + (" -o %s " % sdout) if sdout else s
        s = s + (" -e %s " % sderr) if sderr else s
        s = s + " --exclusive " if exclusive else s
        s = s + (" -J %s " % job_name) if job_name else s
        cmd = s+cmd
    else:
        logger.debug("srun not found")

    return runAndCheck(cmd, "Error in srun command") 

def sbatch(cmd_list, sdout="", sderr=""):
    # logger = logging.getLogger("analysis.functions")
    sdout = sdout if sdout else "/dev/null"
    sderr = sderr if sderr else "/dev/null"
    run = "sbatch -o %s -e %s <<EOF\n#!/bin/sh\n%s\nEOF" % (sdout, sderr, "\n".join(cmd_list))
    # logger.debug("Invoking sbatch: "+run)
    return runAndCheck(run, "Error in sbatch command") 


def srunRscript(cmd, ncores=1, mem="", sdout="", sderr=""):
    ''' Run the given R script via srun.

    This is a wrapper to functions::srun, which simply adds
    the Rscript executable to the front of the command.

    --cmd - the command to be run
    --ncores - the number of cores
    --mem - the memory in Gb
    --log_file - the desired output file for sdout and sderr redirection
    '''
    logger = logging.getLogger("analysis.functions")
    cmd = '/usr/bin/Rscript '+ cmd
    logger.debug("Invoking Rscript: %s " % cmd)
    srun(cmd, ncores, mem, sdout, sderr)


def testLogging(logger):
    '''Testing for the logger
    ''' 
    logger.info("Testing logging info level")
    logger.debug("Testing logging debug level")
    logger.error("Testing logging error level")
    try:
        logger.info("Testing stack trace")
        raise RuntimeError
    except Exception, err:
        logger.exception("Expected exception")

def srun_is_installed():
    '''Tests if srun is installed. 

        Uses 'command -v' for POSIX compliance; works in sh and bash.
        When srun is not found, the result will be 1, else null.
    '''
    cmd = 'command -v srun >/dev/null 2>&1 || { echo "1" >&2; }'
    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    output = p.stdout.read()
    return(output != "1\n")
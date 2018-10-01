''' Generic methods for bioinformatic pipeline running

    These are designed to be used with a slurm cluster
'''
import os
import errno
import logging
import subprocess

__version__ = 'v01'

class SlurmRunner:

    _logger = ""

    def __init__(self, logger_hierarchy=""):
        '''If there is a logger hierarchy, append the name of the module
            and fetch the logger. Otherwise, create a logger instance for
            the module
        '''
        concat = "." if logger_hierarchy else ""
        _logger = logging.getLogger(logger_hierarchy+concat+__name__)

    def runAndCheck(self, cmd, msg):
        '''Run the given command and check the return value.

        If the return value of the command is not 0 (success),
        print the given error message and return
        exit code 1.

        --cmd - the command to be run
        --msg - the message to print if the command returns an error
        '''
        _logger.debug("Invoking: " + cmd)
        try:
            retcode = subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT)
            if retcode != 0:
                _logger.debug("System call returned "+str(retcode))
                _logger.error(msg+". Quitting.")
                return(1)
        except OSError as e:
            _logger.exception("Error invoking command")
            print >>sys.stderr, "Execution exception:", e
            return(1)

    def srun(self, cmd, ncores="1", mem="", sdout="", sderr=""):
        '''Run the given command through srun.

        If srun is available, the command will be run
        with the given cores and memory. If srun is not
        installed, the command is invoked directly.

        --cmd - the command to be run
        --ncores - the number of cores
        --mem - the memory in Gb
        --sdout - the desired output file for sdout redirection (optional)
        --sderr - the desired output file for sderr redirection (optional)
        '''
        if(srun_is_installed()):
            s = "srun -c %s " % (str(ncores))
            s = s + " --mem="+str(mem)+"G " if mem else s
            s = s + (" -o %s " % sdout) if sdout else s
            s = s + (" -e %s " % sderr) if sderr else s
            cmd = s+cmd
        else:
            _logger.debug("srun not found")

        runAndCheck(cmd, "Error in command") 

    def srunRscript(self, cmd, ncores=1, mem="", sdout="", sderr=""):
        ''' Run the given R script via srun.

        This is a wrapper to functions::srun, which simply adds
        the Rscript executable to the front of the command.

        --cmd - the command to be run
        --ncores - the number of cores
        --mem - the memory in Gb
        --log_file - the desired output file for sdout and sderr redirection
        '''
        cmd = '/usr/bin/Rscript '+ cmd
        _logger.debug("Invoking Rscript: %s " % cmd)
        srun(cmd, ncores, mem, sdout, sderr)

    def srun_is_installed(self):
        '''Tests if srun is installed. 

            Uses 'command -v' for POSIX compliance; works in sh and bash.
            When srun is not found, the result will be 1, else null.
        '''
        cmd = 'command -v srun >/dev/null 2>&1 || { echo "1" >&2; }'
        p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
        output = p.stdout.read()
        return(output != "1\n")
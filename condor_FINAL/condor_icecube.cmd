###############################################
# Submit description file for condor-example.py
###############################################
# The "executable", "input", and the list in transfer_input_files are
# always copied by condor. It will look in the directory "initialdir"
# for these files. If there's no initialdir, it will look in the
# directory from which you submit the job.

# The program we're going to execute. As with most
# condor tasks, we're going to submit a shell script
# that will organize the environment for the 'real'
# program, condor-example.py
executable     = /nevis/amsterdam/share/dcianci/IcecubePack/condor_FINAL/condor_icecube.sh

# The program that will be executed by the shell script.
transfer_input_files = icecubescan_FINAL_DIAGNOSTIC,/nevis/amsterdam/share/dcianci/IcecubePack/data.tar.gz

# The arguments to the executable. If you look at the
# script condor-example.sh, you'll see it expects
# one argument. See below for the value of $(Process).
arguments      = $(Process)

# Leave these unchanged. 
universe       = vanilla
should_transfer_files = YES
when_to_transfer_output = ON_EXIT

# Requirements are a restriction on what kinds of machines can run
# your job. Here's an example; it means to make sure the job runs on a
# 64-bit computer.
requirements = ( Arch == "X86_64" )

# You always want to specify where the program's
# text output will go, where its error messages
# will go, and where condor's log messages will go.
TAG=timestamper_dec22
output         = /a/data/westside/dcianci/condor-icecube-${TAG}_$(Process).out
error          = /a/data/westside/dcianci/condor-icecube_${TAG}_$(Process).err
log            = /a/data/westside/dcianci/condor-icecube_${TAG}_$(Process).log

# This keeps you from getting 10,000 emails if you submit 10,000
# jobs... unless each of those 10,000 jobs crash!
notification   = Error
     
# The number following "queue" is how many times this process will be
# submitted to condor. The value of $(Process) will be 0 for the first
# process, 1 for the second, etc.  
queue 1

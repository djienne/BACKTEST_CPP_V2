#!/usr/bin/env python

import platform
from subprocess import call
from functools import partial
from multiprocessing.dummy import Pool

EXECUTABLE = "./F_3EMA_SRSI_ATR.exe"

# number of CPU to use (for local run only)
NB_CPU_LOCAL = 2

print(f'Number of threads {NB_CPU_LOCAL}')

# defining a list of commands to be run in parallel

commands = []

# loops over required initial parameters list

for _ in range(NB_CPU_LOCAL):
    commands.append("./" + EXECUTABLE)

# print(commands)
command_number = len(commands)
print(f'Number of commands required {command_number}')

####################################
computer_name = platform.node()

nb_thread = NB_CPU_LOCAL  # number of threads (cpu) to run

# Making an array where each element is the list of command for a given thread

pool = Pool(nb_thread)  # to be always set to 1 for this MPI case
for i, returncode in enumerate(pool.imap(partial(call, shell=True), commands)):
    if returncode != 0:
        print("%d command failed: %d" % (i, returncode))

#!/bin/bash

HEADER="\033[95m"
OKBLUE="\033[94m"
OKGREEN="\033[92m"
WARNING="\033[93m"
FAIL="\033[91m"
ENDC="\033[0m"
BOLD="\033[1m"
UNDERLINE="\033[4m"

echo "# FINGERPRINTIPY DIRECTORY" >> ~/.bashrc
echo 'export PATH=$PATH:'"$(pwd)" >> ~/.bashrc
echo "" >> ~/.bashrc

echo -e '['${OKGREEN}'Success'${ENDC}'] Finished linking the script.'


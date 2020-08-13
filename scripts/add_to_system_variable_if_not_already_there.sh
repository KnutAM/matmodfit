#!/bin/bash
# Script to add folder to variable if not already included
# Input 1 is system variable
# Input 2 is folder to add

system_variable=$1
folder=$2
if [[ :$"$system_variable": != *:"$folder":* ]] ; then
    # Modify the $HOME/.[bash_]profile file to contain the necessary environmental variables
    ADD_STRING="export "$system_variable"=$"$system_variable:"$folder"

    # We need to check if .bash_profile exists. If it does, it will override content of .profile
    # We should therefore not create it if it doesn't exist.
    # Using the bashrc should also work fine, and is therefore favoured over .profile if .bashrc exists
    
    if [ -f $HOME/.bash_profile ]; then
        path_file="$HOME/.bash_profile"
    elif [ -f $HOME/.bashrc ]; then
        path_file="$HOME/.bashrc"
    else
        path_file="$HOME/.profile"
    fi;
    
    echo "$ADD_STRING" >>$path_file
fi;
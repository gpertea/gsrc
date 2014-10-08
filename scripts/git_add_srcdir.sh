#!/bin/sh
if [[ -z "$1" ]]; then
 echo "Usage: git_add_srcdir.sh <subdir>"
 exit 1
fi
cat << EOF > $1/.gitignore
# Ignore everything in this directory
*
# Except these files
!*.c
!*.cpp
!*.h
!*.hh
!*.pl
!*.py
!*.sh
!Makefile
EOF

git add "$1"
mv $1/.gitignore $1/add.gitignore

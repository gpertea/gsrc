#!/bin/sh
# Manages privoxy daemon, allows switching between different configs if
# setup right.  Either use the -k option to kill any running privoxy, or
# supply the name of a privoxy config file name under PRIVOXYRC to start
# privoxy with.

# where actions, templates, etc. are stored
# need to chdir into here so paths work properly
PRIVOXYCONF=$HOME/ann/misc/cfg/privoxy

# where config files are stored.  these need to
# reference PRIVOXYCONF directory properly
PRIVOXYRC=$PRIVOXYCONF
DEFAULTCFG=config
# where to find privoxy application
PRIVOXY=$HOME/ann/misc/privoxy

RETURN=0

# sanity checks
if [ ! -x $PRIVOXY ]; then
  echo "error: no executable program: $PRIVOXY"
  exit 1
fi

if [ ! -e $PRIVOXYCONF ]; then
  echo "error: no such directory: $PRIVOXYCONF"
  exit 1
fi
if [ ! -d $PRIVOXYCONF ]; then
  echo "error: not a directory: $PRIVOXYCONF"
  exit 1
fi

#if [ -z "$1" ]; then
for opt in "-h" "-help" "--h" "--help"; do
  if [ "x"$opt = "x"$1 ]; then
    echo "usage: `basename $0` [-k] configname"
    exit 1
  fi
done
# routine to find, kill privoxy
killit () {
  ps xwwo pid,command | while read pid command; do
    if echo $command | grep -- "$PRIVOXY" >/dev/null; then
      for signal in "TERM" "INT" "HUP" "KILL"; do
        kill -$signal $pid
        RETURN=$?
        if [ $RETURN -eq 0 ]; then
          break
        fi
        echo "failed: kill $signal $pid" >&2
        sleep 1
      done
    fi
  done
}

OPT=
while getopts "k" OPT; do
  case $OPT in
    k)
      killit
      exit $RETURN
    ;;
  esac
done
shift $(($OPTIND - 1))

MODE=$1
shift

if [ -z "$MODE" ]; then
 MODE=$DEFAULTCFG
fi

if [ ! -s $PRIVOXYRC/$MODE ]; then
 echo "Error: config file not found ($PRIVOXYRC/$MODE)"
 exit 1
fi
killit
cd $PRIVOXYRC
#exec $PRIVOXY --no-daemon $MODE
$PRIVOXY --no-daemon $MODE

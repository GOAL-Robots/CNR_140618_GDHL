#!/bin/bash

set -e 

rm -fr *Rout *log  END matlab.out

FIRST=$1
COLLECT=$2
NEW=$3
ANALYSIS=$4
SAMPLES=$5
STORING=$6
ALL_SAMPLES=$7; ALL_SAMPLES=${ALL_SAMPLES:=0}


if [ "$#" != 6 ] && [ "$#" != 7 ] ; then
  echo
  echo "  usage: $0 <first> <collect> <new> <analysis> <samples> <storing> [all]"
  echo
  exit;
fi

if [ "$#" != 0 ] && [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
  echo
  echo "  usage: $0 <first> <collect> <new> <analysis> <samples> <storing> [all]"
  echo
  exit;
fi



MATLAB="$(which matlab) -nodesktop -nosplash -nodisplay" 

PROCS=6
TRIALS=10
NTRIALS=20
GA_SRC=ga_analysis
GA_TST=ga_test
GA_TST_ign=ga_test_ignore_components
EXP_DATA=exp_data
SOURCES=../code
CURRENT=$(pwd)

let P[1]=8;    
let P[2]=28;    
let P[3]=56;    
let P[4]=70;    
let P[5]=56;    
let P[6]=28;    
let P[7]=8;    
let P[8]=1;    

exps=" BiPoo2001 CassenaerLaurent2007 FroemkeDan2002  \
      # WittenbergWang20062 WoodinGangulyPoo2003  \
      # NishiyamaHongKato2000 TzounopoulosKimTrussell20042  \
      # ZhangTaoPoo1998 ZhouAckerWhite2005reg  \
      # ZhouAckerWhite2005broad HaasNowotnyAbarbanel2006"


# number of processes
count=0;
if [ "$FIRST" == 1 ]; then
    for exp_name in  $exps; do
        for t in $(seq 1 8); do
            for i in $(seq 1 ${P[t]}); do
                for k in $(seq 1 $TRIALS); do
                    let count+=1
                done 
            done
        done
    done
    echo "number of processes $count";
fi



if [ "$FIRST" == 1 ]; then

    echo "genetic algorithms"
    rm -fr ${SOURCES}/data ${SOURCES}/data-* 

    for exp_name in  $exps; do
        for t in $(seq 1 8); do
            for i in $(seq 1 ${P[t]}); do
                for k in $(seq 1 $TRIALS); do

                    cd $SOURCES/matlab
                    sleep .3;
                    echo "$MATLAB -r \" ${GA_SRC}($t,$i,'$exp_name'); exit;\""
                    $MATLAB -r " ${GA_SRC}($t,$i,'$exp_name'); exit;" >> ${CURRENT}/logs/matlab_log 2>&1 &
                    cd ${CURRENT}

                    while [ "$(ps x| grep MATLAB|wc -l)" -gt "$PROCS" ]; 
                    do 
                        printf "\rwaiting..."; 
                        sleep 1; 
                    done
                    printf "\n";

                done 
            done
        done
    done

fi

while [  "$(ps axf | grep -v grep|grep MATLAB)" != "" ]; do
    echo >/dev/null;
done

if [ "$COLLECT" == 1 ]; then


    rm -fr *.png >/dev/null 2>&1
    rm -fr data_all >/dev/null 2>&1
    rm -fr sample-* >/dev/null 2>&1


    echo "collect data"

    if [ ! -d "${SOURCES}/old" ]; then
        mkdir ${SOURCES}/old
    fi
    if [ -d "${SOURCES}/data" ]; then
        mv ${SOURCES}/data old
    fi  

    mkdir  ${SOURCES}/data
    for f in ${SOURCES}/matlab/data-*; 
    do
        mv $f/* ${SOURCES}/data
        rm -fr $f
    done

    for f in $(find ${SOURCES}/data); 
    do 
        if [ -f "$f" ];  then 
            echo $(echo $(basename $f) | sed -e"s/-/ /g" | awk '{print $2,$3,$4}') $(cat $f|head -n 1)>>data_all; 
        fi; 
    done

fi

if [ "$ANALYSIS" == 1 ]; then

    echo "data analysis"

    ${SOURCES}/R/data_analysis.r >${CURRENT}/logs/r_log 2>&1

fi

if [ "$NEW" == 1 ]; then

    echo "new data" 

    for f  in sample-*; 
    do 
        exp_name=$(echo -n $f|sed -e"s/sample-//")  

        while read s 
        do

            t=$(echo $s| awk '{print $2}')
            i=$(echo $s| awk '{print $3}')

            for k in $(seq 1 $NTRIALS); do
                cd $SOURCES/matlab
                sleep .1;
                $MATLAB -r " ${GA_SRC}($t,$i,'$exp_name'); exit;" >> ${CURRENT}/logs/matlab_log 2>&1 &
                cd ${CURRENT}

                while [ "$(ps x| grep MATLAB|wc -l)" -gt "$PROCS" ]; 
                do 
                    printf "\rwaiting..."; sleep 1; 
                done
                printf "\n";
            done

        done <<< "$(cat $f| grep "^3 ")"

    done

fi

if [ "$SAMPLES" == 1 ]; then

    echo "best examples"

    for f  in sample-*; 
    do 
        nm=$(echo -n $f|sed -e"s/sample-//")  

        s="$(cat $f| grep "^1 \|^2 ")"

        while read r
        do
            if ! [ -z "$(echo $r| grep NA)" ]; then continue; fi


            kin=$(echo $r| awk '{print $1}')
            vnr=$(echo $r| awk '{print $2}')
            ind=$(echo $r| awk '{print $3}')
            see=$(echo $r| awk '{print $4}')

            expr="${GA_TST}($vnr,$ind,'$nm',$see); exit;"

            if [ "$(echo $r| grep mask )" != "" ]; then
                gpp=$(echo $r| awk '{print $7}')
                gnp=$(echo $r| awk '{print $8}')
                gpn=$(echo $r| awk '{print $9}')
                gnn=$(echo $r| awk '{print $10}')
                eps=$(echo $r| awk '{print $11}')
                esp=$(echo $r| awk '{print $12}')
                ens=$(echo $r| awk '{print $13}')
                esn=$(echo $r| awk '{print $14}')
                expr="${GA_TST_ign}($vnr,$ind,'$nm',$gpp,$gnp,$gpn,$gnn,$eps,$esp,$ens,$esn,$see); exit;"
            fi

            if [ "$kin" == 1 ]; then
                cd ${SOURCES}/matlab
                $MATLAB -r "$expr" > ${CURRENT}/logs/matlab.log 2>&1 
                cd ${CURRENT}
                mv  -f ${SOURCES}/matlab/curve-${nm}-* . 
            fi

            if [ "$ALL_SAMPLES" == 1 ]; then
                if [ "$kin" == 2 ]; then

                    echo "begin $vnr,$ind,$nm,$see"

                    cd ${SOURCES}
                    $MATLAB -r "$expr" > ${CURRENT}/logs/matlab.log 2>&1 
                    cd ${CURRENT}
                    mv  -f ${SOURCES}/matlab/curve-${nm}-* . 

                    echo "end"

                fi
            fi

        done <<< "$(echo -e "$s")";

    done

fi

if [ "$STORING" == 1 ]; then

    echo "store data analysis"

    analysisdir="data_analysis-$(date +%Y%m%d)" 
    if ! [ -d  "$analysisdir" ]; then
        mkdir $analysisdir
    fi

    mv data_all curve-* ga-*.png sample-* *_params *bestscores $analysisdir

fi

touch END

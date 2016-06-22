#!/bin/bash


if [ ! -n $1 ] || [ -z $1 ] || [ $1 == "?" ] || [ $1 == "h" ]; then  
    echo "Usage: lcheck.sh <path> [assem|launch|splash]"
    echo "       log check for assemble.py, launch.py and splashdown.py logs." 
    exit 1
fi

type="unknown"
if [ $# -ge 1 ]; then
    path=$1
    if [ $# -eq 2 ]; then
        type=$2
    else
        if [[ $path == *"assem"* ]]; then 
            type="assem"
        elif [[ $path == *"launch"* ]]; then 
            type="launch"
        elif [[ $path == *"splash"* ]]; then 
            type="splash"
        elif [[ $path == *"scrub"* ]]; then 
            type="scrub"
        fi
    fi
fi
    
if [ "$type" != "assem" ]  && [ "$type" != "a" ] \
&& [ "$type" != "launch" ] && [ "$type" != "l" ] \
&& [ "$type" != "splash" ] && [ "$type" != "s" ] \
&& [ "$type" != "scrub" ]; then
    echo "Usage: lcheck.sh <path> [assem|launch|splash|scrub]"
    echo "       log check for assemble.py, launch.py and splashdown.py logs." 
    exit 1
fi
if [[ $path == *"ENCSR"* ]] && [[ $path == *".log"* ]]; then
    cat $path
fi
path=${path%/*}

exps=`ls -1 $path/*.log | wc -l`
echo "Exps:      " $exps
echo "ERRORS:    " `grep -i ERROR $path/*.log | grep -v Gateway | grep -v "500 Internal" | wc -l` "---" `grep -i ERROR $path/*.log | grep -v Gateway` | grep -v "500 Internal"
echo "WARNINGS:  " `grep -i WARN $path/*.log | grep -v retry | grep -v "Waiting 60" | wc -l` `grep -i WARN $path/*.log | grep -v retry` | grep -v "Waiting 60" 

if [ "$type" == "assem" ] || [ "$type" == "a" ]; then

    echo "Replicates: 1_1:"`grep rep1_1 $path/*.log | grep -v fastq | wc -l` "1_2:"`grep rep1_2 $path/*.log | grep -v fastq | wc -l` \
                     "2_1:"`grep rep2_1 $path/*.log | grep -v fastq | wc -l` "2_2:"`grep rep2_2 $path/*.log | grep -v fastq | wc -l`
    echo "Finished:  " `grep finished $path/*.log | wc -l`
    grep Processed $path/*.log
    grep "No files need" $path/*.log

elif [ "$type" == "launch" ] || [ "$type" == "l" ]; then  

    echo "TEST ONLY: " `grep "TEST ONLY" $path/*.log | wc -l`
    echo "Lift-off:  " `grep "We have liftoff" $path/*.log | wc -l`
    steps=`grep "will be run" $path/*.log | wc -l`
    flows=`ls -1 $path/*.log | wc -l`
    # integer math:
    steps_per_10x=`expr $steps \* 10 / $flows`
    steps_per=`expr $steps / $flows`
    steps_per_frac=`expr $steps_per_10x - $steps_per \* 10`
    #echo "Steps: ${steps}"
    #echo "Flows: ${flows}"
    #echo "steps_per_10x: ${steps_per_10x}"
    #echo "steps_per_frac: ${steps_per_frac}"
    echo "Steps Per:  ${steps_per}.${steps_per_frac}"
    #grep Processed $path/*.log

elif [ "$type" == "splash" ] || [ "$type" == "s" ]; then

    echo "RETRIES:   " `grep -i WARN $path/*.log | grep -E 'retry|Waiting' | wc -l`
    echo "HALTS:     " `grep -i HALTING $path/*.log | wc -l` `grep -i HALTING $path/*.log`
    files=`grep Flagged $path/*.log | wc -l`
    files_per_10x=`expr $files \* 10 / $exps`
    files_per=`expr $files / $exps`
    files_per_frac=`expr $files_per_10x - $files_per \* 10`
    echo "Files:     " $files "per exp:  ${files_per}.${files_per_frac}"
    qcs=`grep "Posted qc_metric" $path/*.log | wc -l`
    qcs_per_10x=`expr $qcs \* 10 / $exps`
    qcs_per=`expr $qcs / $exps`
    qcs_per_frac=`expr $qcs_per_10x - $qcs_per \* 10`
    echo "QC metrics:" $qcs "per exp:  ${qcs_per}.${qcs_per_frac}"
    echo "Finished:  " `grep finished $path/*.log | wc -l`
    grep cost $path/*.log | grep -v analysis | sed s/^.*\\/// | sed s/\.log:cost://  ## If ENCSR000ZZZ.log:cost:
    grep Processed $path/*.log | sed s/^.*\\/// | sed s/\.log//
    #grep 'pipeline completed' $path/*.log | sed s/^.*\\/// | sed s/\.log//
    echo "COMPLETED: " `grep 'pipeline completed' $path/*.log | uniq | wc -l`

elif [ "$type" == "scrub" ]; then

    expected=`grep scrubbing $path/*.log | grep -v b\.log | awk '{ sum += $3 } END { print sum }'`
    exps=`grep processed $path/*.log | wc -l`
    files=`grep processed $path/*.log | awk '{ sum += $8 } END { print sum }'` 
    files_per_10x=`expr $files \* 10 / $exps`
    files_per=`expr $files / $exps`
    files_per_frac=`expr $files_per_10x - $files_per \* 10`
    echo "Experiments:" $exps "of" $expected
    echo "Files:      " $files "per exp:  ${files_per}.${files_per_frac}"

fi  


#!/bin/bash

if [ ! -n $1 ] || [ -z $1 ] || [ $1 == "?" ] || [ $1 == "h" ]; then  
    echo "Usage: lcheck.sh <logpath> [assem|launch|splash|scrub|rollover]"
    echo "       log check for assemble.py, launch.py and splashdown.py logs." 
    exit 1
fi

type="unknown"
if [ $# -ge 1 ]; then
    logpath=$1
    if [ $# -eq 2 ]; then
        type=$2
    else
        if [[ $logpath == *"assem"* ]]; then 
            type="assem"
        elif [[ $logpath == *"launch"* ]]; then 
            type="launch"
        elif [[ $logpath == *"splash"* ]]; then 
            type="splash"
        elif [[ $logpath == *"scrub"* ]]; then 
            type="scrub"
        elif [[ $logpath == *"rollover"* ]]; then 
            type="rollover"
        fi
    fi
fi
    
if [ "$type" != "assem" ]  && [ "$type" != "a" ] \
&& [ "$type" != "launch" ] && [ "$type" != "l" ] \
&& [ "$type" != "splash" ] && [ "$type" != "s" ] \
&& [ "$type" != "scrub" ] && [ "$type" != "rollover" ] ; then
    echo "Usage: lcheck.sh <logpath> [assem|launch|splash|scrub|rollover]"
    echo "       log check for assemble.py, launch.py and splashdown.py logs." 
    exit 1
fi
if [[ $logpath == *"ENCSR"* ]] && [[ $logpath == *".log"* ]]; then
    cat $logpath
fi
logpath=${logpath%/*}

exps=`ls -1 $logpath/*.log | wc -l`
if [ "$type" != "scrub" ] &&  [ "$type" != "rollover" ]; then
    echo "Exps:      " $exps
    echo "ERRORS:    " `grep -i ERROR $logpath/*.log | grep -v Gateway | grep -v "500 Internal" | wc -l` "---" `grep -i ERROR $logpath/*.log | grep -v Gateway` | grep -v "500 Internal"
    echo "WARNINGS:  " `grep -i WARN $logpath/*.log | grep -v retry | grep -v "Waiting 60" | wc -l` `grep -i WARN $logpath/*.log | grep -v retry` | grep -v "Waiting 60" 
fi

if [ "$type" == "assem" ] || [ "$type" == "a" ]; then

    echo "Replicates: 1_1:"`grep rep1_1 $logpath/*.log | grep -v fastq | wc -l` "1_2:"`grep rep1_2 $logpath/*.log | grep -v fastq | wc -l` \
                     "2_1:"`grep rep2_1 $logpath/*.log | grep -v fastq | wc -l` "2_2:"`grep rep2_2 $logpath/*.log | grep -v fastq | wc -l`
    echo "Finished:  " `grep finished $logpath/*.log | wc -l`
    grep Processed $logpath/*.log
    grep "No files need" $logpath/*.log

elif [ "$type" == "launch" ] || [ "$type" == "l" ]; then  

    echo "TEST ONLY: " `grep "TEST ONLY" $logpath/*.log | wc -l`
    echo "Lift-off:  " `grep "We have liftoff" $logpath/*.log | wc -l`
    steps=`grep "will be run" $logpath/*.log | wc -l`
    flows=`ls -1 $logpath/*.log | wc -l`
    # integer math:
    steps_per_10x=`expr $steps \* 10 / $flows`
    steps_per=`expr $steps / $flows`
    steps_per_frac=`expr $steps_per_10x - $steps_per \* 10`
    #echo "Steps: ${steps}"
    #echo "Flows: ${flows}"
    #echo "steps_per_10x: ${steps_per_10x}"
    #echo "steps_per_frac: ${steps_per_frac}"
    echo "Steps Per:  ${steps_per}.${steps_per_frac}"
    #grep Processed $logpath/*.log

elif [ "$type" == "splash" ] || [ "$type" == "s" ]; then

    echo "RETRIES:   " `grep -i WARN $logpath/*.log | grep -E 'retry|Waiting' | wc -l`
    echo "HALTS:     " `grep -i HALTING $logpath/*.log | wc -l` `grep -i HALTING $logpath/*.log`
    files=`grep Flagged $logpath/*.log | wc -l`
    files_per_10x=`expr $files \* 10 / $exps`
    files_per=`expr $files / $exps`
    files_per_frac=`expr $files_per_10x - $files_per \* 10`
    echo "Files:     " $files "per exp:  ${files_per}.${files_per_frac}"
    qcs=`grep "Posted qc_metric" $logpath/*.log | wc -l`
    qcs_per_10x=`expr $qcs \* 10 / $exps`
    qcs_per=`expr $qcs / $exps`
    qcs_per_frac=`expr $qcs_per_10x - $qcs_per \* 10`
    echo "QC metrics:" $qcs "per exp:  ${qcs_per}.${qcs_per_frac}"
    echo "Finished:  " `grep finished $logpath/*.log | wc -l`
    grep cost $logpath/*.log | grep -v analysis | sed s/^.*\\/// | sed s/\.log:cost://  ## If ENCSR000ZZZ.log:cost:
    grep Processed $logpath/*.log | sed s/^.*\\/// | sed s/\.log//
    #grep 'pipeline completed' $logpath/*.log | sed s/^.*\\/// | sed s/\.log//
    echo "COMPLETED: " `grep 'pipeline completed' $logpath/*.log | uniq | wc -l`

elif [ "$type" == "scrub" ]; then

    testing=`grep "Would remove" $logpath/*.log | wc -l`
    if [ $testing -gt 0 ]; then
        testing="Testing "
    else
        testing=""
    fi
    echo "${testing}Scrub:"
    expected=`grep scrubbing $logpath/*.log | grep -v b\.log | grep -v c\.log | awk '{ sum += $3 } END { print sum }'`
    examined=`grep "Examining ENCODE" $logpath/*.log | wc -l`
    exps=`grep processed $logpath/*.log | wc -l`
    if [ $exps -eq 0 ]; then
        exps=$examined # --test
    fi
    #files=`grep processed $logpath/*.log | awk '{ sum += $8 } END { print sum }'` 
    files=`grep "Removing file" $logpath/*.log | wc -l`
    if [ $files -eq 0 ]; then
        files=`grep "Would remove file" $logpath/*.log | wc -l`
        if [ $files -eq 0 ] && [ $exps -gt 0 ]; then
            files=-999 # --remove_all was requested
        fi
    fi
    echo "${testing}Exps:       $exps complete.  $examined of $expected examined."
    if [ $files -ne -999 ]; then
        files_per=0
        files_per_frac=0
        if [ $exps -gt 0 ]; then
            files_per_10x=`expr $files \* 10 / $exps`
            files_per=`expr $files / $exps`
            files_per_frac=`expr $files_per_10x - $files_per \* 10`
        fi
        echo "${testing}Files:      $files scrubbed.  ${files_per}.${files_per_frac} per experiment."
    else
        echo "${testing}Files:      --remove_all scrubs experiment directories"
    fi
    grep "(finished)" $logpath/*.log

elif [ "$type" == "rollover" ]; then

    testing=`grep "Would move" $logpath/*.log | wc -l`
    if [ $testing -gt 0 ]; then
        testing="Testing "
    else
        testing=""
    fi
    echo "${testing}Rollover:"
    #logpath="logs/scrub/lrna/2016-06-21/"
    expected=`grep Requested $logpath/*.log | grep -v b\.log | grep -v c\.log | awk '{print $3}'`
    examined=`grep Working $logpath/*.log | wc -l`
    exps=`grep Found $logpath/*.log | grep -v "Found 0" | wc -l`
    files=`grep Moving $logpath/*.log | awk '{ sum += $4 } END { print sum }'`
    if [ "$files" == "" ]; then
        files=`grep "Would move" $logpath/*.log | awk '{ sum += $5 } END { print sum }'`
    fi
    files_per_10x=`expr $files \* 10 / $exps`
    files_per=`expr $files / $exps`
    files_per_frac=`expr $files_per_10x - $files_per \* 10`
    echo "${testing}Exps:       $exps with fastqs.  $examined of $expected examined."
    echo "${testing}Files:      $files moving.  ${files_per}.${files_per_frac} per experiment."
    grep "(finished)" $logpath/*.log

fi  


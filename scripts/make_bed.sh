#!/bin/bash
(( num = $1 ))
(( max = $2 ))
while (( max >= num )); do
    (( prev = num  ))
    (( num = prev + ${3}-1 ))
    if (( num >= max )); then
            echo $4 $prev $max
    else
            echo $4 $prev $num
    fi
done

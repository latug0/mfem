#!/bin/bash
for i in $(seq 6); do 
   time ./lamine -o 3 -t $i -m cube_2mat_per.msh 1>out${i} 2>err${i}
   echo -n "  $i "; tail -n 3 err${i}
done

#!/bin/bash

#./runfit ../Out/toyMC_txt/200_100/0.txt ../Out/template_txt/background.txt ../Out/template_txt/signal.txt 2 500_1000
#"
#for bkg in 10 100 500 1000 5000 10000
#do
#    for sign in 1 5 10 20 40 60 100 200 400 600 1000
#    do
#	for i in {0..100}
#	do
#	    ./runfit ../Out_V2/toyMC_txt/${sign}_${bkg}/${i}.txt ../Out_V2/template_txt/background.txt ../Out_V2/template_txt/signal.txt ${i} ${sign}_${bkg}
#	done
#    done
#done

for bkg in 10 #100 500 1000 5000 10000 100000 1000000
do
    for sign in 400 600 1000
    do
	for i in {0..5} 
	do
	    #echo ./runfit ../Out_V2/toyMC_txt/${sign}_${bkg}/${i}.txt ../Out_V2/template_txt/background.txt ../Out_V2/template_txt/signal.txt ${i} ${sign}_${bkg}
	    ./runfit ../toyMC_txt/${sign}_${bkg}/${i}.txt ../background.txt ../signal.txt ${i} ${sign}_${bkg} NID
	done
    done
done




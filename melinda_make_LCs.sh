#!/bin/bash
INLIST=lc_cads.dat
LC=/Users/msoaresfurtado/Dropbox/Research_NGC6791/Notebooks/b0i0d10/
INPATH=/Users/msoaresfurtado/Dropbox/Research_NGC6791/Notebooks/subtraction_photometry_output/photometry_output_
for base in $(cat $INLIST);
do
    JDs=${base: -5}
    #JDs=$(grep $BASE G559_jd_recond.txt | cut -d ' ' -f 2,3)
    #cat $frames | grtrans -i - --col-xy 2,3 --col-out 16,17,18 --input-transformation $TRANSFILE -o - |  awk -v JDs="$JDs" '{print JDs,$0}' ; 
    #cat $frames | grtrans -i - --col-xy 2,3 --col-out 16,17,18 --input-transformation $TRANSFILE -o - ; 
    #/home/hatuser/venv/bin/python2.7 photbin2txt.py $frame | grtrans -i - --col-xy 2,3 --col-out 19,20,21 --input-transformation $TRANSFILE -o - ; 
    #cat $frames | grtrans -i - --col-xy 2,3 --col-out 19,20,21 --input-transformation $TRANSFILE -o - ; 
    grep 'G' $INPATH$base.out | awk -v JDs="$JDs" '{print JDs,$0}'
done | grcollect - --col-base 2 --prefix $LC --extension rlc

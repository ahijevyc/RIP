#!/bin/csh
#This should make 2 ps files for each location listed between the parentheses.

foreach f (330_260 350_260 370_260 390_260 410_260 430_260 330_270 350_270 370_270 390_270 410_270 430_270 330_280 350_280 370_280 390_280 410_280 430_280 330_290 350_290 370_290 390_290 410_290 430_290 330_300 350_300 370_300 390_300 410_300 430_300 330_310 350_310 370_310 390_310 410_310 430_310 )
psmerge ts_tke_${f}_tot.ps ts_bvf_${f}_tot.ps ts_ri_${f}_tot.ps ts_shr_${f}_tot.ps | psnup -4 > $f.1.ps
psmerge ts_thp_${f}_tot.ps ts_www_${f}_tot.ps ts_thp_${f}_upr.ps ts_wsp_${f}_upr.ps | psnup -4 > $f.2.ps
end 

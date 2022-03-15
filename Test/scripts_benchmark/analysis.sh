#!/bin/bash

nsims=100
forcelist="1 5 10 15 20"
fbase="/media/Storage/MADna/v0.0/AMBER_protocol"
fscripts="/media/Storage/MADna/v0.0/AMBER_protocol/Scripts_check"
imin=1000

seqlist="AA AC AG AT CG GG A4GGA4 A4TA4 A8GG A8T DDD DUE G4AAG4 G4CG4 TATA TFBS"
#seqlist="AA"
handle=4
for sequenza in $seqlist; do 
	lseq=$(awk '{print length($1)}' $fbase/Simulations/$sequenza/initialization/sequence.dat)
	i0=$(echo $handle | awk '{print $1+1}')
	i1=$(echo $handle $lseq | awk '{print $2-$1}')

	#### extract parameters
	for sim in `seq 1 $nsims`; do 
		echo "$sequenza $sim" >&2
		iid=1
		for force in $forcelist; do
			fold="$fbase/Simulations/$sequenza/f$force/sim$sim"
			${fscripts}/Compute_ext_twist_crook $fold/dump.lammpstrj $RANDOM $i0 $i1 > $fold/ext_twist_crook & 
			pids[$iid]=$!
			iid=$((iid+1))
		done
		fold="$fbase/Simulations/$sequenza/f1/sim$sim"
		${fscripts}/ComputeHelicalParameters $fold/dump.lammpstrj $RANDOM $fold/hel &
		pids[$iid]=$!
		iid=$((iid+1))
		fold="$fbase/Simulations/$sequenza/f1/sim$sim"
		${fscripts}/ComputeSBBS $fold/dump.lammpstrj > $fold/SBBS &
		pids[$iid]=$!
		iid=$((iid+1))
		for pid in ${pids[*]}; do
			wait $pid
		done
	done

	#### compute averages
	mkdir -p $fbase/Simulations/$sequenza/results
	echo "Compute extension for $sequenza" >&2
	for force in $forcelist; do
		echo $force $(for sim in `seq 1 $nsims`; do awk -v imin=$imin 'NR>imin{print $1}' $fbase/Simulations/$sequenza/f$force/sim$sim/ext_twist_crook | media; done | media)
	done > $fbase/Simulations/$sequenza/results/extension
	
	echo "Compute twist for $sequenza" >&2
	for force in $forcelist; do
		echo $force $(for sim in `seq 1 $nsims`; do awk -v imin=$imin 'NR>imin{print $2}' $fbase/Simulations/$sequenza/f$force/sim$sim/ext_twist_crook | media; done | media)
	done > $fbase/Simulations/$sequenza/results/cumulate_twist

	echo "Compute crookedness for $sequenza" >&2
	for force in $forcelist; do
		echo $force $(for sim in `seq 1 $nsims`; do awk -v imin=$imin 'NR>imin{print $3}' $fbase/Simulations/$sequenza/f$force/sim$sim/ext_twist_crook | media; done | media)
	done > $fbase/Simulations/$sequenza/results/crookedness
	
	echo "Compute diameter for $sequenza and f=1pN" >&2
	for sim in `seq 1 $nsims`; do awk -v imin=$imin 'NR>imin{print $2}' $fbase/Simulations/$sequenza/f1/sim$sim/hel_diameter | media; done | media > $fbase/Simulations/$sequenza/results/diameter
	
	for xxx in hrise htwist; do
		echo "Compute $xxx for $sequenza and f=1pN" >&2
		echo "#sequence step i1 i2 $xxx error" > $fbase/Simulations/$sequenza/results/$xxx
		for i1step in `seq $((handle+1)) $((lseq-handle-1))`; do
			var=$(for sim in `seq 1 $nsims`; do awk -v imin=$imin -v i=$i1step 'NR>imin{print $(i+1)}' $fbase/Simulations/$sequenza/f1/sim$sim/hel_$xxx | media; done | media)
			step=$(awk -v i=$i1step '{print substr($1,i,2)}' $fbase/Simulations/$sequenza/initialization/sequence.dat)
			echo "$sequenza $step $i1step $((i1step+1)) $var"
		done >> $fbase/Simulations/$sequenza/results/$xxx
	done
	
	for xxx in major_groove_width minor_groove_width major_groove_depth minor_groove_depth; do
		echo "Compute $xxx for $sequenza and f=1pN" >&2
		echo "#sequence nucleotide i $xxx error" > $fbase/Simulations/$sequenza/results/$xxx
		for i1 in `seq $((handle+1)) $((lseq-handle))`; do
			var=$(for sim in `seq 1 $nsims`; do awk -v imin=$imin -v i=$i1 'NR>imin&&$(i+1)>-1{print $(i+1)}' $fbase/Simulations/$sequenza/f1/sim$sim/hel_$xxx | media; done | media)
			bp=$(awk -v i=$i1 '{print substr($1,i,1)}' $fbase/Simulations/$sequenza/initialization/sequence.dat)
			echo "$sequenza $bp $i1 $var"
		done >> $fbase/Simulations/$sequenza/results/$xxx
	done

	for xxx in SBBS; do
                echo "Compute $xxx for $sequenza and f=1pN" >&2
                echo "#sequence nucleotide i $xxx error" > $fbase/Simulations/$sequenza/results/$xxx
                for i1 in `seq $((handle+1)) $((lseq-handle))`; do
                        var=$(for sim in `seq 1 $nsims`; do awk -v imin=$imin -v i=$i1 'NR>imin{print $(i+1)}' $fbase/Simulations/$sequenza/f1/sim$sim/$xxx | media; done | media)
                        bp=$(awk -v i=$i1 '{print substr($1,i,1)}' $fbase/Simulations/$sequenza/initialization/sequence.dat)
                        echo "$sequenza $bp $i1 $var"
                done >> $fbase/Simulations/$sequenza/results/$xxx
        done

	#### compute elastic constants
        echo "Compute A1 for $sequenza" >&2
	res1=$(${fscripts}/FitWWithMyError_General.py $fbase/Simulations/$sequenza/results/extension)
	A1=$(echo $res1 | awk '{print $1}')
	dA1=$(echo $res1 | awk '{print $2}')
	L0=$(echo $res1 | awk '{print $3}')
	dL0=$(echo $res1 | awk '{print $4}')

        echo "Compute A2 for $sequenza" >&2
	res2=$(${fscripts}/FitWWithMyError_General.py $fbase/Simulations/$sequenza/results/cumulate_twist)
	A2=$(echo $res2 | awk '{print $1}')
	dA2=$(echo $res2 | awk '{print $2}')

	echo "Compute theta_var for $sequenza" >&2
        thetamedio=$(head -1 $fbase/Simulations/$sequenza/results/cumulate_twist | awk '{print $2}')
        for sim in `seq 1 $nsims`; do
                awk -v imin=$imin -v t0=$thetamedio 'NR>imin{print ($2-t0)^2}' $fbase/Simulations/$sequenza/f1/sim$sim/ext_twist_crook | media
        done | awk '{print $1}' > $fbase/Simulations/$sequenza/results/thetavar_list
        $fscripts/bootstrap $fbase/Simulations/$sequenza/results/thetavar_list $RANDOM 10000 | tail -1 | awk '{print $2,$3}' > $fbase/Simulations/$sequenza/results/thetavar

        echo "Compute kbeta for $sequenza" >&2
	res3=$(${fscripts}/FitWWithMyError_Crook_New.py $fbase/Simulations/$sequenza/results/crookedness)
	kbeta=$(echo $res3 | awk '{print $1}')
	dkbeta=$(echo $res3 | awk '{print $2}')

        echo "Compute elastic constants for $sequenza" >&2
	echo $L0 $dL0 $A1 $dA1 $A2 $dA2 $kbeta $dkbeta $(cat $fbase/Simulations/$sequenza/results/thetavar) | awk '\
		{L0=$1; dL0=$2; A1=$3; dA1=$4; A2=$5; dA2=$6; kb=$7; dkb=$8; tvar=$9; dtvar=$10;\
		s = L0/A1; ds=dL0/A1+L0/A1^2*dA1;\
		c = 4.14*L0/tvar; dc=4.14*dL0/tvar + 4.14*L0/tvar^2*dtvar;\
		g = -A2*c*s/L0; dg=dA2*c*s/L0 + A2*dc*s/L0 + A2*c*ds/L0 + A2*c*s/L0^2*dL0;
		print "Stilde(pN) ",s,ds;\
                print "C(pN nm^2) ",c,dc;\
                print "g(pN nm) ",g,dg;
                print "kbeta(pN) ",kb,dkb}' > $fbase/Simulations/$sequenza/results/constants
done

seqlist="DNAall"
handle=3
for sequenza in $seqlist; do 
	lseq=$(awk '{print length($1)}' $fbase/Simulations/$sequenza/initialization/sequence.dat)
	i0=$(echo $handle | awk '{print $1+1}')
	i1=$(echo $handle $lseq | awk '{print $2-$1}')

	#### extract parameters
	for sim in `seq 1 $nsims`; do 
		iid=1
		for force in $forcelist; do
			fold="$fbase/Simulations/$sequenza/f$force/sim$sim"
			echo "$sequenza ${force}pN $sim extension, twist and crookedness" >&2
			${fscripts}/Compute_ext_twist_crook $fold/dump.lammpstrj $RANDOM $i0 $i1 > $fold/ext_twist_crook & 
			pids[$iid]=$!
			iid=$((iid+1))
		done
		echo "$sequenza 1pN $sim helical parameters" >&2
		fold="$fbase/Simulations/$sequenza/f1/sim$sim"
		${fscripts}/ComputeHelicalParameters $fold/dump.lammpstrj $RANDOM $fold/hel &
		pids[$iid]=$!
		iid=$((iid+1))
		echo "$sequenza 1pN $sim SBBS" >&2
		fold="$fbase/Simulations/$sequenza/f1/sim$sim"
		${fscripts}/ComputeSBBS $fold/dump.lammpstrj > $fold/SBBS &
		pids[$iid]=$!
		iid=$((iid+1))
		for pid in ${pids[*]}; do
			wait $pid
		done
	done

	#### compute averages
	mkdir -p $fbase/Simulations/$sequenza/results
	echo "Compute extension for $sequenza" >&2
	for force in $forcelist; do
		echo $force $(for sim in `seq 1 $nsims`; do awk -v imin=$imin 'NR>imin{print $1}' $fbase/Simulations/$sequenza/f$force/sim$sim/ext_twist_crook | media; done | media)
	done > $fbase/Simulations/$sequenza/results/extension
	
	echo "Compute twist for $sequenza" >&2
	for force in $forcelist; do
		echo $force $(for sim in `seq 1 $nsims`; do awk -v imin=$imin 'NR>imin{print $2}' $fbase/Simulations/$sequenza/f$force/sim$sim/ext_twist_crook | media; done | media)
	done > $fbase/Simulations/$sequenza/results/cumulate_twist

	echo "Compute crookedness for $sequenza" >&2
	for force in $forcelist; do
		echo $force $(for sim in `seq 1 $nsims`; do awk -v imin=$imin 'NR>imin{print $3}' $fbase/Simulations/$sequenza/f$force/sim$sim/ext_twist_crook | media; done | media)
	done > $fbase/Simulations/$sequenza/results/crookedness
	
	echo "Compute diameter for $sequenza and f=1pN" >&2
	for sim in `seq 1 $nsims`; do awk -v imin=$imin 'NR>imin{print $2}' $fbase/Simulations/$sequenza/f1/sim$sim/hel_diameter | media; done | media > $fbase/Simulations/$sequenza/results/diameter
	
	for xxx in hrise htwist; do
		echo "Compute $xxx for $sequenza and f=1pN" >&2
		echo "#sequence step i1 i2 $xxx error" > $fbase/Simulations/$sequenza/results/$xxx
		for i1step in `seq $((handle+1)) $((lseq-handle-1))`; do
			var=$(for sim in `seq 1 $nsims`; do awk -v imin=$imin -v i=$i1step 'NR>imin{print $(i+1)}' $fbase/Simulations/$sequenza/f1/sim$sim/hel_$xxx | media; done | media)
			step=$(awk -v i=$i1step '{print substr($1,i,2)}' $fbase/Simulations/$sequenza/initialization/sequence.dat)
			echo "$sequenza $step $i1step $((i1step+1)) $var"
		done >> $fbase/Simulations/$sequenza/results/$xxx
	done
	
	for xxx in major_groove_width minor_groove_width major_groove_depth minor_groove_depth; do
		echo "Compute $xxx for $sequenza and f=1pN" >&2
		echo "#sequence nucleotide i $xxx error" > $fbase/Simulations/$sequenza/results/$xxx
		for i1 in `seq $((handle+1)) $((lseq-handle))`; do
			var=$(for sim in `seq 1 $nsims`; do awk -v imin=$imin -v i=$i1 'NR>imin&&$(i+1)>-1{print $(i+1)}' $fbase/Simulations/$sequenza/f1/sim$sim/hel_$xxx | media; done | media)
			bp=$(awk -v i=$i1 '{print substr($1,i,1)}' $fbase/Simulations/$sequenza/initialization/sequence.dat)
			echo "$sequenza $bp $i1 $var"
		done >> $fbase/Simulations/$sequenza/results/$xxx
	done

	for xxx in SBBS; do
                echo "Compute $xxx for $sequenza and f=1pN" >&2
                echo "#sequence nucleotide i $xxx error" > $fbase/Simulations/$sequenza/results/$xxx
                for i1 in `seq $((handle+1)) $((lseq-handle))`; do
                        var=$(for sim in `seq 1 $nsims`; do awk -v imin=$imin -v i=$i1 'NR>imin{print $(i+1)}' $fbase/Simulations/$sequenza/f1/sim$sim/$xxx | media; done | media)
                        bp=$(awk -v i=$i1 '{print substr($1,i,1)}' $fbase/Simulations/$sequenza/initialization/sequence.dat)
                        echo "$sequenza $bp $i1 $var"
                done >> $fbase/Simulations/$sequenza/results/$xxx
        done

	#### compute elastic constants
        echo "Compute A1 for $sequenza" >&2
	res1=$(${fscripts}/FitWWithMyError_General.py $fbase/Simulations/$sequenza/results/extension)
	A1=$(echo $res1 | awk '{print $1}')
	dA1=$(echo $res1 | awk '{print $2}')
	L0=$(echo $res1 | awk '{print $3}')
	dL0=$(echo $res1 | awk '{print $4}')

        echo "Compute A2 for $sequenza" >&2
	res2=$(${fscripts}/FitWWithMyError_General.py $fbase/Simulations/$sequenza/results/cumulate_twist)
	A2=$(echo $res2 | awk '{print $1}')
	dA2=$(echo $res2 | awk '{print $2}')

        echo "Compute theta_var for $sequenza" >&2
	thetamedio=$(head -1 $fbase/Simulations/$sequenza/results/cumulate_twist | awk '{print $2}')
	for sim in `seq 1 $nsims`; do
		awk -v imin=$imin -v t0=$thetamedio 'NR>imin{print ($2-t0)^2}' $fbase/Simulations/$sequenza/f1/sim$sim/ext_twist_crook | media
	done | awk '{print $1}' > $fbase/Simulations/$sequenza/results/thetavar_list
	$fscripts/bootstrap $fbase/Simulations/$sequenza/results/thetavar_list $RANDOM 10000 | tail -1 | awk '{print $2,$3}' > $fbase/Simulations/$sequenza/results/thetavar

        echo "Compute kbeta for $sequenza" >&2
	res3=$(${fscripts}/FitWWithMyError_Crook_New.py $fbase/Simulations/$sequenza/results/crookedness)
	kbeta=$(echo $res3 | awk '{print $1}')
	dkbeta=$(echo $res3 | awk '{print $2}')

        echo "Compute elastic constants for $sequenza" >&2
	echo $L0 $dL0 $A1 $dA1 $A2 $dA2 $kbeta $dkbeta $(cat $fbase/Simulations/$sequenza/results/thetavar) | awk '\
		{L0=$1; dL0=$2; A1=$3; dA1=$4; A2=$5; dA2=$6; kb=$7; dkb=$8; tvar=$9; dtvar=$10;\
		s = L0/A1; ds=dL0/A1+L0/A1^2*dA1;\
		c = 4.14*L0/tvar; dc=4.14*dL0/tvar + 4.14*L0/tvar^2*dtvar;\
		g = -A2*c*s/L0; dg=dA2*c*s/L0 + A2*dc*s/L0 + A2*c*ds/L0 + A2*c*s/L0^2*dL0;
		print "Stilde(pN) ",s,ds;\
		print "C(pN nm^2) ",c,dc;\
		print "g(pN nm) ",g,dg;
		print "kbeta(pN) ",kb,dkb}' > $fbase/Simulations/$sequenza/results/constants
done



####### collect results for comparison
for xxx in major_groove_depth major_groove_width minor_groove_depth minor_groove_width SBBS; do 
	echo "Collect results for $xxx" >&2
	echo "#sequence nucleotide i simulation error benchmark error" > Simulations/Comparison/$xxx
	for x in AA AC AG AT CG GG A4GGA4 A4TA4 A8GG A8T DDD DUE G4AAG4 G4CG4 TATA TFBS DNAall; do
		paste Simulations/$x/results/$xxx Benchmark/$x/$xxx | grep -v "#" | awk '{print $1,$2,$3,$4,$5,$9,$10}'
	done >> Simulations/Comparison/$xxx
done

for xxx in hrise htwist; do 
	echo "Collect results for $xxx" >&2
	echo "#sequence nucleotide i1 i2 simulation error benchmark error" > Simulations/Comparison/$xxx
	for x in AA AC AG AT CG GG A4GGA4 A4TA4 A8GG A8T DDD DUE G4AAG4 G4CG4 TATA TFBS DNAall; do
		paste Simulations/$x/results/$xxx Benchmark/$x/$xxx | grep -v "#" | awk '{print $1,$2,$3,$4,$5,$6,$11,$12}'
	done >> Simulations/Comparison/$xxx
done

for xxx in diameter; do 
	echo "Collect results for $xxx" >&2
	echo "#sequence simulation error benchmark error" > Simulations/Comparison/$xxx
	for x in AA AC AG AT CG GG A4GGA4 A4TA4 A8GG A8T DDD DUE G4AAG4 G4CG4 TATA TFBS DNAall; do
		paste Simulations/$x/results/$xxx Benchmark/$x/$xxx | awk -v x=$x '{print x,$0}'
	done >> Simulations/Comparison/$xxx
done

for xxx in Stilde C g kbeta; do 
	echo "Collect results for $xxx" >&2
	echo "#sequence simulation error benchmark error" > Simulations/Comparison/$xxx
	for x in AA AC AG AT CG GG A4GGA4 A4TA4 A8GG A8T DDD DUE G4AAG4 G4CG4 TATA TFBS DNAall; do
		paste <(awk '{print $1, $(NF-1),$NF}' Simulations/$x/results/constants) <(awk '{print $(NF-1),$NF}' Benchmark/$x/constants) | grep $xxx | awk -v x=$x '{print x,$2,$3,$4,$5}'
	done >> Simulations/Comparison/$xxx
done


####### create graphs
for xxx in major_groove_depth major_groove_width minor_groove_depth minor_groove_width; do
echo "Creating graphs for $xxx" >&2
title=$(echo $xxx | awk -F "_" '{print $1,$2,$3}')
gnuplot << daqui
set xlabel "Simulations"
set ylabel "Benchmark"
set title "$title (nm)"
set term png
set output "Simulations/Comparison/scatter_plot_${xxx}.png"
p "Simulations/Comparison/$xxx" u 4:6:5:7 t '' w xye pt 7 ps 2 lw 2, x t '' lw 2
daqui
done

xxx="hrise"
echo "Creating graphs for $xxx" >&2
gnuplot << daqui
set xlabel "Simulations"
set ylabel "Benchmark"
set title "$xxx (nm)"
set term png
set output "Simulations/Comparison/scatter_plot_${xxx}.png"
p "Simulations/Comparison/$xxx" u 5:7:6:8 t '' w xye pt 7 ps 2 lw 2, x t '' lw 2
daqui

for xxx in htwist; do
echo "Creating graphs for $xxx" >&2
gnuplot << daqui
set xlabel "Simulations"
set ylabel "Benchmark"
set title "$xxx (deg)"
set term png
set output "Simulations/Comparison/scatter_plot_${xxx}.png"
p "Simulations/Comparison/$xxx" u 5:7:6:8 t '' w xye pt 7 ps 2 lw 2, x t '' lw 2
daqui
done

for xxx in SBBS; do
echo "Creating graphs for $xxx" >&2
gnuplot << daqui
set xlabel "Simulations"
set ylabel "Benchmark"
set title "$xxx (deg)"
set term png
set output "Simulations/Comparison/scatter_plot_${xxx}.png"
p "Simulations/Comparison/$xxx" u 4:6:5:7 t '' w xye pt 7 ps 2 lw 2, x t '' lw 2
daqui
done

for xxx in Stilde kbeta; do
echo "Creating graphs for $xxx" >&2
gnuplot << daqui
set xlabel "Simulations"
set ylabel "Benchmark"
set title "$xxx (pN)"
set term png
set output "Simulations/Comparison/scatter_plot_${xxx}.png"
p "Simulations/Comparison/$xxx" u 2:4:3:5 t '' w xye pt 7 ps 2 lw 2, x t '' lw 2
daqui
done

for xxx in C; do
echo "Creating graphs for $xxx" >&2
gnuplot << daqui
set xlabel "Simulations"
set ylabel "Benchmark"
set title "$xxx (pN nm^2)"
set term png
set output "Simulations/Comparison/scatter_plot_${xxx}.png"
p "Simulations/Comparison/$xxx" u 2:4:3:5 t '' w xye pt 7 ps 2 lw 2, x t '' lw 2
daqui
done

for xxx in g; do
echo "Creating graphs for $xxx" >&2
gnuplot << daqui
set xlabel "Simulations"
set ylabel "Benchmark"
set title "$xxx (pN nm)"
set term png
set output "Simulations/Comparison/scatter_plot_${xxx}.png"
p "Simulations/Comparison/$xxx" u 2:4:3:5 t '' w xye pt 7 ps 2 lw 2, x t '' lw 2
daqui
done

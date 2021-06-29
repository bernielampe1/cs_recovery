for num in 2 50 100 125 150
do
    for j in 0 1 
    do
	rm -f loqo_res_${num}_${j}_el1.out
	printfiles_el1 $num $j
	time ampl loqo_el1.mod > loqo_res_${num}_${j}_el1.out
    done
done
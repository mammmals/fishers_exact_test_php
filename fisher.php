<?php
//Copyright (c) 2017 Devin K Schweppe, Gygi Lab
//Written by Devin K Schweppe
//Original binomial calculation by Ed Huttlin
//DIstributed under MIT License
	function binomial($N, $x){
		$product = 1;
		for($ii = 0; $ii < $x; $ii++){
			$product = $product * ($N - $ii) / ($ii + 1);
		}
		return $product;
	}
	function log_binomial($N, $x){
		$product = log10(1);
		for($ii = 0; $ii < $x; $ii++){
			$iterated_product = ($N - $ii) / ($ii + 1);
			$product = $product + log10($iterated_product);
		}
		return $product;
	}

	function single_value_fisher($net_go,$all_net,$bp_go,$all_bp) {
		$cell_a=$net_go;
		$cell_b=$bp_go;
		$cell_c=$all_net-$net_go;
		$cell_d=$all_bp-$bp_go;
		
		$a_p_b = $cell_a + $cell_b;
		$c_p_d = $cell_c + $cell_d;
		$a_p_c = $cell_a + $cell_c;
		
		$n_total=$cell_a+$cell_b+$cell_c+$cell_d;
		
		$bn_1=binomial($a_p_b,$cell_a);
		$bn_2=binomial($c_p_d,$cell_c);
		$bn_3=binomial($n_total,$a_p_c);
		
		$fisher_p = ($bn_1 * $bn_2)/$bn_3;
		
		return($fisher_p);
	}
		
	function fisher($net_go,$all_net,$bp_go,$all_bp) {
		$fisher_p=0;		
		foreach(range($net_go,$all_net) as $fisher_var){

			$cell_a=$fisher_var;
			$cell_b=$bp_go;
			$cell_c=$all_net-$fisher_var;
			$cell_d=$all_bp-$bp_go;
			
			$a_p_b = $cell_a + $cell_b;
			$c_p_d = $cell_c + $cell_d;
			$a_p_c = $cell_a + $cell_c;
			
			$n_total=$cell_a+$cell_b+$cell_c+$cell_d;
			
			$bn_1=binomial($a_p_b,$cell_a);
			$bn_2=binomial($c_p_d,$cell_c);
			$bn_3=binomial($n_total,$a_p_c);

			$fisher_p += ($bn_1 * $bn_2)/$bn_3;
		}
		return($fisher_p);
	}
	echo(fisher(1,50,100,10000)."<br />");
	
	function log_fisher($net_go,$all_net,$bp_go,$all_bp) {
		$fisher_p=0;		
		foreach(range($net_go,$all_net) as $fisher_var){

			$cell_a=$fisher_var;
			$cell_b=$bp_go;
			$cell_c=$all_net-$fisher_var;
			$cell_d=$all_bp-$bp_go;
			
			$a_p_b = $cell_a + $cell_b;
			$c_p_d = $cell_c + $cell_d;
			$a_p_c = $cell_a + $cell_c;
			
			$n_total=$cell_a+$cell_b+$cell_c+$cell_d;
			
			$bn_1=log_binomial($a_p_b,$cell_a);
			$bn_2=log_binomial($c_p_d,$cell_c);
			$bn_3=log_binomial($n_total,$a_p_c);

			$fisher_p_adder = $bn_1 + $bn_2 - $bn_3;
			$fisher_p = $fisher_p + pow(10,$fisher_p_adder);
		}
		return($fisher_p);
	}
	
	echo(log_fisher(1,50,100,10000));
?>
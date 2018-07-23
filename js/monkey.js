function InputTextBoxes()
{
	var N;
	var E;
    var mbuffer;
    var T;
    var mfilter_per_entry;
    var P;
    var leveltier;  // tiered is 0, leveled is 1, lazy_level is 2, fluid_LSM-Tree is 3
    var isLeveled;
		var isOptimalFPR; // 0 -> false; 1 -> true
		var fluidK;
		var fluidZ;
		var Mu;
		var s;
		var Y;
		var w;
		var r;
		var v;
		var qL;
		var qS;
}


function parseInputTextBoxes()
{
	var parsedBoxes = new InputTextBoxes();
    parsedBoxes.N = parseInt(document.getElementById("N").value.replace(/\D/g,''),10);
    parsedBoxes.E = parseInt(document.getElementById("E").value.replace(/\D/g,''),10);
    parsedBoxes.mbuffer = parseFloat(document.getElementById("mbuffer").value.replace(/\D/g,''))*1048576;
    parsedBoxes.T = parseInt(document.getElementById("T").value.replace(/\D/g,''), 10);
    parsedBoxes.mfilter_per_entry = parseFloat(document.getElementById("mfilter_per_entry").value);
    parsedBoxes.P = parseInt(document.getElementById("P").value.replace(/\D/g,''), 10);
		parsedBoxes.L = parseInt(document.getElementById("L").value.replace(/\D/g,''), 10);
		parsedBoxes.s = parseInt(document.getElementById("s").value.replace(/\D/g,''), 10);
    parsedBoxes.isLeveled = isRadioLeveled("ltradio");  // tiered is 0, leveled is 1
    parsedBoxes.leveltier = getRadioValueByName("ltradio");
		parsedBoxes.isOptimalFPR = getRadioValueByName("fpr_radio");
		parsedBoxes.fluidK = parseInt(document.getElementById("Fluid LSM-Tree K").value.replace(/\D/g,''), 10);
		parsedBoxes.fluidZ = parseInt(document.getElementById("Fluid LSM-Tree Z").value.replace(/\D/g,''), 10);
		parsedBoxes.Mu = parseFloat(document.getElementById("Mu").value)
		parsedBoxes.w = parseFloat(document.getElementById("w").value);
		parsedBoxes.r = parseFloat(document.getElementById("r").value);
		parsedBoxes.v = parseFloat(document.getElementById("v").value);
		parsedBoxes.qL = parseFloat(document.getElementById("qL").value);
		parsedBoxes.qS = parseFloat(document.getElementById("qS").value);
    return parsedBoxes;
}

function Filter() {
    var nokeys;
    var fp;
    var mem;
}

function getTotalNonExistingPointLookupCost(i, L, filter_array, N, T, B, Y, K, Z, s, Mu, isOptimalFPR){
	var EULER = 2.71828182845904523536;
	var result = 0;
	for (var j = 0; j < L - Y ; j++)
	{
			var val = calc_R(filter_array[j]);
			if (j < L - Y - 1) {  // tiering
					result += val * K;
			} else {
					result += val * Z;
			}
	}
	return result+Y*Z;
}

function getTotalExistingPointLookupCost(i, L, filter_array, N, T, B, Y, K, Z, s, Mu, isOptimalFPR){
	var EULER = 2.71828182845904523536;
	var result = 0;
	for(j=0;j<L;j++){
		var tmp = calc_R(filter_array[j]);
		if(j < L-Y-1){
			tmp *= K;
		}else if(j < L-1){
			tmp *= Z;
		}else{
			tmp *= Math.floor(Z/2);
			tmp += 1;
		}
		result += tmp;
	}
	return result;
}

function getTotalShortRangeLookupCost(i, L, filter_array, N, T, B, Y, K, Z, s, Mu, isOptimalFPR){
	return (Y+1)*Z + K*(L - Y - 1);
}

function getTotalLongRangeLookupCost(i, L, filter_array, N, T, B, Y, K, Z, s, Mu, isOptimalFPR){
	var total = 0;
	for(j = 1; j < L; j++){
		total += s/B/Math.pow(T, L - j)/Mu;
	}
	total += Z*s/B/Mu
	return Math.ceil(total);
}

function getTotalUpdateCost(i, L, filter_array, N, T, B, Y, K, Z, s, Mu, isOptimalFPR){
	return ((T - 1)*((L - Y - 1)/(K + 1) + (Y + 1.0)/(Z + 1)))/Mu/B;
}

function getTotalSpaceAmp(i, L, filter_array, N, T, B, Y, K, Z, s, Mu, isOptimalFPR){
	var result = 0.0;
	for(j = 1; j < L; j++){
		result += Math.pow(T, j - L)
	}
	result += Z - 1;
	return result;
}

function getLeveledLevel(i, L, filter_array, N, T, B,Y,  K, Z, s, Mu, isOptimalFPR){
	return i;
}

function getLeveledNonExistingPointLookupCost(i, L, filter_array, N, T, B,Y, K, Z, s, Mu, isOptimalFPR){
	var result;
	result = calc_R(filter_array[i - 1]);
	if(i < L-Y){
		result *= K;
	}else if(i == L - Y){
		result *= Z;
	}else{
		result = Z;
	}
	return result;
}


function getLeveledExistingPointLookupCost(i, L, filter_array, N, T, B, Y, K, Z, s, Mu, isOptimalFPR){
	var EULER = 2.71828182845904523536;
	var result;
	var result;
	result = calc_R(filter_array[i - 1]);
	if(i < L-Y){
		result *= K;
	}else if(i < L){
		result *= Z;
	}else{
		result *= Math.floor(Z/2);
		result += 1;
	}
	return result;
}

function getLeveledShortRangeLookupCost(i, L, filter_array, N, T, B, Y, K, Z, s, Mu, isOptimalFPR){
	if(i >= L-Y){
		return Z;
	}else{
		return K;
	}
}

function getLeveledLongRangeLookupCost(i, L, filter_array, N, T, B, Y, K, Z, s, Mu, isOptimalFPR){
	if(i == L){
		return Math.floor(Z*s/B/Mu);
	}else{
		return Math.floor(s/B/Math.pow(T, L - i)/Mu);
	}
}

function getLeveledUpdateCost(i, L, filter_array, N, T, B, Y, K, Z, s, Mu, isOptimalFPR){
	if(i >= L-Y){
		return (T - 1)/(B*(Z + 1)*Mu);
	}else{
		return (T - 1)/(B*(K + 1)*Mu);
	}
}

function getLeveledSpaceAmp(i, L, filter_array, N, T, B, Y, K, Z, s, Mu, isOptimalFPR){
	if(i == L){
		return Z-1;
	}else{
		return Math.pow(T, i - L);
	}
}

function myCeil(n, i){
	var n_tmp = n;
	var j = 0;
	while(n_tmp < 1){
		n_tmp *= 10;
		j += 1;
	}
	n_tmp *= Math.pow(10, i);
	n_tmp = Math.ceil(n_tmp);
	n_tmp /= Math.pow(10, i+j);
	return n_tmp;
}

function myFloor(n, i){
	var n_tmp = n;
	var j = 0;
	while(n_tmp < 1){
		n_tmp *= 10;
		j += 1;
	}
	n_tmp *= Math.pow(10, i);
	n_tmp = Math.floor(n_tmp);
	n_tmp /= Math.pow(10, i+j);
	return n_tmp;
}

var log10 = Math.log(10);
var lastTreeType = 0;
function getSignificantDigitCount(n) {
    n = Math.abs(String(n).replace(".", "")); //remove decimal and make positive
    if (n == 0) return 0;
    while (n != 0 && n % 10 == 0) n /= 10; //kill the 0s at the end of n

    return Math.floor(Math.log(n) / log10) + 1; //get number of digits
}

function getMostSignificantDigit(d)
{
	var i=0;
	if (d>1 || d<=0)
		return -1;
	else
	{
		while (d<1)
		{
			d=d*10;
			i++;
		}
		return i;
	}
}

function formatBytes(bytes,decimals) {
   if(bytes == 0) return '0 Byte';
   var k = 1024; // or 1024 for binary
   var dm = decimals + 1 || 3;
   var sizes = ['Bytes', 'KB', 'MB', 'GB', 'TB', 'PB', 'EB', 'ZB', 'YB'];
   var i = Math.floor(Math.log(bytes) / Math.log(k));
   return parseFloat((bytes / Math.pow(k, i)).toFixed(dm)) + ' ' + sizes[i];
}

function numberWithCommas(x) {
    return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
}

function pad(n, width, z) {
  z = z || '0';
  n = n + '';
  return n.length >= width ? n : new Array(width - n.length + 1).join(z) + n;
}

function isRadioLeveled(radioName)
{
    var leveltier=getRadioValueByName(radioName);
    return (leveltier==1);
}

function getRadioValueByName(radioName)
{
    var radios = document.getElementsByName(radioName);
    var val;
    for(var i = 0; i < radios.length; i++){
        if(radios[i].checked){
            val = radios[i].value;
        }
    }
    return val;
}

function myLog(x, base)
{
    return Math.log(x) / Math.log(base);
}


function calcBits(filter) {
    var denom=Math.pow(Math.log(2),2);
    return filter.fp == 0 ? 0 : ( - filter.nokeys * Math.log(filter.fp) ) / denom;
}


function calcTotalMemBits(filtersArray) {
    var total = 0;
    for (var i = 0; i < filtersArray.length ; i++) {
        // printf("%d   %d  %f \n", i,  rates[i].size, rates[i].false_positive_rate);
        var val = calcBits(filtersArray[i]);

        total += val;
    }
    return total;
}

function getRmax(lt,N,E,B,T)
{
    var R;
    if (lt==0)
        R=Math.ceil(myLog(N*E/B,T));
    else if (lt==1)
        R=(T-1)*Math.ceil(myLog(N*E/B,T));
    return R;
}




function calc_R(f)
{
    var denom = Math.pow(Math.log(2), 2);
		if(f.mem == undefined || isNaN(f.mem)){
			return 1;
		}
    var value = Math.exp(-(f.mem / f.nokeys) * denom);
    return value;
}

function calc_Y(b, K, Z, T, L){
	var flag = true;
	bb_tmp = Math.log(K)-Math.log(Z)+Math.log(T)/(T-1)
	bb_A = b*Math.pow(Math.log(2),2)*(Math.pow(T, L+1)-1) + bb_tmp*T
	bb_C = T*Math.log(T)
	bb_B = bb_tmp+Math.log(T)
	var Y = 0;
	if(flag){

		delta = Number.MAX_VALUE;
		for(var i=0;i<=L;i++){
			bb_X = L - i;
			tmp_delta = bb_A - (bb_B*Math.pow(T, bb_X) - bb_C*bb_X)
			if(tmp_delta > 0 && tmp_delta < delta){
				delta = tmp_delta;
				Y = i;
			}
		}
	}

	return Math.max(Y, 0);
}


// TODO for tiering, we now assume there are (T-1) runs in the last level are all equal to each other in size,
// but if the last level is not full to capacity, then there may in reality be less than (T-1) runs, and they
// may have different sizes. To fix this problem, we can insert multiple runs per filter to the filters_array
// until the number of remaining keys is 0.
function initFilters(N,E,mbuffer,T,K,Z,Y, mfilter_bits,P,leveltier, isOptimalFPR) {

    var filter_array = [];
    var remainingKeys=N-mbuffer/E;
    var level=0;
    //Calculate the number of keys per level in a almost-full in each level LSM tree
    while (remainingKeys>0)
    {
        level++;
        var levelKeys=Math.ceil(Math.min(Math.pow(T,level)*mbuffer/E,N));
        var newFilter = new Filter();
        newFilter.nokeys=levelKeys;
        newFilter.fp=0.0;
        // console.log("New remaining keys: "+(remainingKeys-levelKeys))
        if (remainingKeys-levelKeys<0)
            newFilter.nokeys=remainingKeys;
        //console.log(newFilter.nokeys)
        filter_array.push(newFilter);
        remainingKeys=remainingKeys-levelKeys;
        // console.log(levelKeys)
    }
		var mem_level = filter_array.length - Y;
		for (var i=0;i<mem_level;i++)
		{
				filter_array[i].mem=mfilter_bits/mem_level;
		}
		for (var i=mem_level;i<filter_array.length;i++){
			filter_array[i].mem= 0;
		}
/*
    //Initialize the memory per level to be equal
		if(!isOptimalFPR){
			for (var i=0;i<filter_array.length;i++)
	    {
	        filter_array[i].mem=mfilter_bits/filter_array.length;
	    }
	    return filter_array;
		}else{
			var EULER = 2.71828182845904523536;
			var C = Math.pow(Math.log(2), 2);
			var L = filter_array.length;
			var X = 1.0/C*(Math.log(T/(T - 1)) + Math.log(K/Z)/T)
			var Y = Math.ceil(Math.log(N*X/mfilter_bits)/Math.log(T))
			if(Y > 0){
				L -= Y;
			}
			var tmp = 0;
			for(var i=1;i<=L;i++){
				if(i == L){
					tmp += Z*filter_array[i-1].nokeys*Math.log(Z);
				}else{
					tmp += K*filter_array[i-1].nokeys*Math.log(K*Math.pow(T, L - i));
				}
			}
			var t = Math.pow(EULER, (tmp - mfilter_bits*C)/N);
			var equal_flag = false;
			var eqaul_mem = 0;
			var left_levels = 1;
			for(var i=1;i<=L;i++){
				var fpr_i;
				var runs;
				if(equal_flag){
					filter_array[i-1].mem = eqaul_mem;
				}else{
					if(i==L){
						fpr_i = Math.max(mfilter_bits, 0);
					}else{
						fpr_i = t/K/Math.pow(T, L-i);
					}
					if(fpr_i > 1){
						equal_flag = true;
						left_levels = L - i + 1;
						eqaul_mem = mfilter_bits/left_levels;
						filter_array[i-1].mem = eqaul_mem;
					}else{
						filter_array[i-1].mem=-filter_array[i-1].nokeys*Math.log(fpr_i)/C;
						mfilter_bits -= filter_array[i-1].mem;
					}
				}

			}
		}*/

		return filter_array;
}


function getBaselineFPassigment(N,E,mbuffer,T,K, Z, Y, mfilter_bits,P,leveltier)
{
	  var THRESHOLD = 1e-15
    var filter_array = initFilters(N,E,mbuffer,T,K, Z, Y, mfilter_bits,P,leveltier, false);

    // console.log(filter_array);
		var mem_level = filter_array.length - Y;
    var limit_on_M=mfilter_bits/mem_level; //amount of memory for filters in bits
    var diff = mfilter_bits;
    var change = true;
    var iteration = 0;

    while (diff > 1) {
        change = false;
        for (var i = 0; i < mem_level - 1; i++) {
            for (var j = i + 1; j < mem_level ; j++) {
							var flag = 0;
                var f1_orig = calc_R(filter_array[i]);
                var f2_orig = calc_R(filter_array[j]);
                // console.log(f1_orig+', '+f2_orig)
                var diff_orig = Math.abs(f1_orig - f2_orig);

                filter_array[i].mem += diff;
                filter_array[j].mem -= diff;

                var f1_new = calc_R(filter_array[i]);
                var f2_new = calc_R(filter_array[j]);
                var diff_new = Math.abs(f1_new - f2_new);
                // console.log(f1_new+', '+f2_new)


                if (diff_orig - diff_new >  THRESHOLD && filter_array[j].mem > 0 && filter_array[i].mem > 0) {
                    change = true;
										continue;
                }
                filter_array[i].mem -= diff * 2;
                filter_array[j].mem += diff * 2;

                f1_new = calc_R(filter_array[i]);
                f2_new = calc_R(filter_array[j]);
                diff_new = Math.abs(f1_new - f2_new);

                if (diff_orig - diff_new>  THRESHOLD	 && filter_array[j].mem > 0 && filter_array[i].mem > 0) {
                    change = true;
                    continue;
                }
								filter_array[i].mem += diff;
                filter_array[j].mem -= diff;

            }
        }
        if (!change) {
            diff /= 2;
        }
        iteration++;
    }
    for (var i = 0; i < filter_array.length; i++) {
        filter_array[i].fp = calc_R(filter_array[i]);
        // console.log(filter_array[i].mem+', '+filter_array[i].fp)
    }
    return filter_array;
}


function getMonkeyFPassigment(N, E, mbuffer, T, K, Z, Y, mfilter_bits, P, leveltier)
{
		var THRESHOLD = 1e-15
    var filter_array = initFilters(N,E,mbuffer,T,K, Z, Y, mfilter_bits,P,leveltier, true);
		var mem_level = filter_array.length - Y;
    // console.log(filter_array);
    var limit_on_M=mfilter_bits/filter_array.length; //amount of memory for filters in bits
    var diff = limit_on_M;
    var change = true;
    var iteration = 0;
    var current_R = eval_R(filter_array, leveltier, T, Y, K, Z);
    var original = current_R;
    var value = 0;
    while (diff > 1) {
        change = false;
        for (var i = 0; i < mem_level - 1; i++) {
            for (var j = i + 1; j < mem_level ; j++) {
							  var flag = 0;
                filter_array[i].mem += diff;
                filter_array[j].mem -= diff;
                value = eval_R(filter_array, leveltier, T, Y, K, Z);
                if (current_R - value > THRESHOLD && value > 0 && filter_array[j].mem > 0 ) {
                    current_R = value;
                    change = true;
										continue;
                }
                filter_array[i].mem -= diff * 2;
                filter_array[j].mem += diff * 2;

                value = eval_R(filter_array, leveltier, T, Y, K, Z);

                if (current_R - value > THRESHOLD && value > 0 && filter_array[i].mem > 0 ) {
                    current_R = value;
                    change = true;
										continue;
                }

								filter_array[i].mem += diff;
                filter_array[j].mem -= diff;

            }
        }
        if (!change) {
            diff /= 2;
        }
        iteration++;
    }

    for (var i = 0; i < filter_array.length; i++) {
        filter_array[i].fp = calc_R(filter_array[i]);
    }

    return filter_array;
}



function eval_R(filters, leveltier, T, Y, K = 9, Z = 1)
{
    var total = 0;
    var n = filters.length;
		if(leveltier == 0){
			K = T - 1;
			Z = T - 1;
		}else if(leveltier == 1){
			K = 1;
			Z = 1;
		}else if(leveltier == 2){
			K = T - 1;
			Z = 1;
		}
    n -= Y + 1;
    for (var i = 0; i < n ; i++)
    {
        var val = calc_R(filters[i]);
				total += val * K;
    }
    if(Y < filters.length){
			var val = calc_R(filters[n]);
	    total += val * Z;
		}


    return total+Y*Z;
}

function reset_button_colors()
{
	var color='#777';
    document.getElementById("scenario1").style.background=color;
    document.getElementById("scenario2").style.background=color;
    document.getElementById("scenario3").style.background=color;
}

function scenario1()
{
    document.getElementById("N").value=numberWithCommas(68719476736); //(10M values)
    document.getElementById("E").value=16;
    document.getElementById("mbuffer").value=2; //in MB
    document.getElementById("T").value=10;
		document.getElementById("L").value=6;
    document.getElementById("mfilter_per_entry").value=10; // bits per element
    document.getElementById("P").value=4096; //in B
    document.getElementsByName("ltradio")[0].checked=true;
    document.getElementsByName("ltradio")[1].checked=false;
		document.getElementById("Fluid LSM-Tree K").value = 1;
		document.getElementById("Fluid LSM-Tree Z").value = 1;
		document.getElementById("Mu").value = 1;
		document.getElementById("s").value = 8192;
		document.getElementById("w").value = 1;
		document.getElementById("r").value = 1;
		document.getElementById("v").value = 1;
		document.getElementById("qL").value = 1;
		document.getElementById("qS").value = 1;
		lastTreeType = 1;

    reset_button_colors()
    document.getElementById("scenario1").style.background='#000000';

    clickbloomTuningButton(true)
}

function scenario2()
{
    document.getElementById("N").value=numberWithCommas(68719476736); //(2^36 values)
    document.getElementById("E").value=16;
    document.getElementById("mbuffer").value=2; //in MB
    document.getElementById("T").value=4;
		document.getElementById("L").value=10;
    document.getElementById("mfilter_per_entry").value= 10 // bits per element
    document.getElementById("P").value=4096; //in B
    document.getElementsByName("ltradio")[0].checked=false;
    document.getElementsByName("ltradio")[1].checked=true;
		document.getElementById("Fluid LSM-Tree K").value = 3;
		document.getElementById("Fluid LSM-Tree Z").value = 3;
		document.getElementById("Mu").value = 1;
		document.getElementById("s").value = 8192;
		document.getElementById("w").value = 1;
		document.getElementById("r").value = 1;
		document.getElementById("v").value = 1;
		document.getElementById("qL").value = 1;
		document.getElementById("qS").value = 1;
		lastTreeType = 0;

    reset_button_colors()
    document.getElementById("scenario2").style.background='#000000';

    clickbloomTuningButton(true)
}

function scenario3()
{
    document.getElementById("N").value=numberWithCommas(68719476736); //(2^36 values)
    document.getElementById("E").value=16;
    document.getElementById("mbuffer").value=2; //in MB
    document.getElementById("T").value=2;
		document.getElementById("L").value=19;
    document.getElementById("mfilter_per_entry").value=0; //0 bits per element
    document.getElementById("P").value=4096; //in B
    document.getElementsByName("ltradio")[0].checked=true;
    document.getElementsByName("ltradio")[1].checked=false;
		document.getElementById("Fluid LSM-Tree K").value = 1;
		document.getElementById("Fluid LSM-Tree Z").value = 1;
		document.getElementById("Mu").value = 1;
		document.getElementById("s").value = 8192;
		document.getElementById("w").value = 1;
		document.getElementById("r").value = 1;
		document.getElementById("v").value = 1;
		document.getElementById("qL").value = 1;
		document.getElementById("qS").value = 1;
		lastTreeType = 1;

    reset_button_colors()
    document.getElementById("scenario3").style.background='#000000';

    clickbloomTuningButton(true)
}


function clickbloomTuningButton(move_to_anchor) {

	var inputParameters = parseInputTextBoxes();

    var N=inputParameters.N;
    var E=inputParameters.E;
    var mbuffer=inputParameters.mbuffer;
    var T=inputParameters.T;
    var mfilter_per_entry=inputParameters.mfilter_per_entry;
    var P=inputParameters.P;
    var leveltier=inputParameters.leveltier;
		var K = inputParameters.fluidK;
		var Z = inputParameters.fluidZ;
		var s = inputParameters.s;
		var isOptimalFPR = inputParameters.isOptimalFPR;
		var Mu = inputParameters.Mu;
		var B = P/E;
		var w=inputParameters.w;
		var v=inputParameters.v;
		var r=inputParameters.r;
		var qL=inputParameters.qL;
		var qS=inputParameters.qS;



    if (document.getElementById("N").value=="" || document.getElementById("E").value=="" || document.getElementById("T").value=="" || document.getElementById("P").value==""
        || document.getElementById("mbuffer").value=="" || document.getElementById("mfilter_per_entry").value=="" || isNaN(N)
        || isNaN(E) || isNaN(mbuffer) || isNaN(T) || isNaN(mfilter_per_entry) || isNaN(P) || isNaN(leveltier) || isNaN(K) || isNaN(Z) || isNaN(Mu) || isNaN(s) || isNaN(w) || isNaN(v) || isNaN(r) || isNaN(qL) || isNaN(qS))
	{
		alert("Some of the input boxes are empty! Either select values for all parameters or run one of the scenarios.");
		return;
	}

	if (leveltier==1) {
		K = 1;
		Z = 1;
		} else if (leveltier==0) {
			K = T - 1;
			Z = T - 1;
		} else if (leveltier==2) {
			K = T - 1;
			Z = 1;
		}
		lastTreeType = leveltier;
		var L = Math.ceil(Math.log(N*E*(T - 1)/mbuffer/T+ 1/T)/Math.log(T));
		var Y = calc_Y(mfilter_per_entry, K, Z, T, L);

    //get BF allocation
		var mfilter_bits = mfilter_per_entry*N;
		var filters;
		if(isOptimalFPR == "0"){
			filters = getBaselineFPassigment(N, E, mbuffer, T, K, Z, Y, mfilter_bits, P, leveltier);
		}else{
			filters = getMonkeyFPassigment(N, E, mbuffer, T, K, Z, Y, mfilter_bits, P, leveltier);
		}

    google.charts.setOnLoadCallback(drawChart);


	//present it nicely in the html!
    var result_div=document.getElementById("result_tuning")
	while (result_div.firstChild) {
	    result_div.removeChild(result_div.firstChild);
	}

	var THRESHOLD=1e-9;

    // var br=document.createElement("br");
    // result_div.appendChild(br);



    var first_col_alignment="center";
	//titles
    var div_new_row=document.createElement("div");
    div_new_row.setAttribute("class","row");
		var div_level = document.createElement("div");
		div_level.setAttribute("class","col-sm-1");
		div_level.setAttribute("style","width:5.333%");
		var p=document.createElement("p");
		p.setAttribute("style","text-align: center;font-size:15px;font-weight:bold")
		p.textContent=("Level")
		div_level.appendChild(p);
		var div_schema = document.createElement("div");
		div_schema.setAttribute("class","col-sm-6");
		div_schema.setAttribute("style","width:51%");
		var div_schema_row = document.createElement("div");
		div_schema_row.setAttribute("class","row");
		var div_title_col2=document.createElement("div");
		div_title_col2.setAttribute("class","col-sm-5")
		var lsm_header2=document.createElement("h5");
		lsm_header2.textContent="Entries per Level";
		lsm_header2.setAttribute("style","text-align: center;")
	div_title_col2.appendChild(lsm_header2);
		div_new_row.appendChild(div_level);
		div_new_row.appendChild(div_schema);
		div_schema.appendChild(div_schema_row);
		div_new_row.appendChild(div_title_col2);


		var div_col_tmp=document.createElement("div");
		div_col_tmp.setAttribute("class","col-sm-4");
	  var p=document.createElement("p");
		p.setAttribute("style","text-align: center;font-size:15px;font-weight:bold")
		p.textContent=("Point Lookup")
		div_col_tmp.appendChild(p);
		div_schema_row.appendChild(div_col_tmp);

		var div_col_tmp=document.createElement("div");
		div_col_tmp.setAttribute("class","col-sm-4")
		var p=document.createElement("p");
		p.setAttribute("style","text-align: center;font-size:15px;font-weight:bold")
		p.textContent=("Range Lookup")
		div_col_tmp.appendChild(p);
		div_schema_row.appendChild(div_col_tmp);

		var div_col_tmp=document.createElement("div");
		div_col_tmp.setAttribute("class","col-sm-2")
		div_col_tmp.setAttribute("style","width:13.33%")
		var p=document.createElement("p");
		var span = document.createElement("span");
		var msg = "In the worst case, an entry participates in O(T/K) merge operations within an active run across each of Levels 1 to L-1, " +
		"and in O(T/Z) merge operations within the active run at Level L"
		span.setAttribute("data-tooltip",msg);
		span.setAttribute("data-tooltip-position","bottom")
		p.setAttribute("style","text-align: center;font-size:15px;font-weight:bold")
		p.textContent=("Update")
		span.appendChild(p);
		div_col_tmp.appendChild(span);
		div_schema_row.appendChild(div_col_tmp);

		var div_col_tmp=document.createElement("div");
		div_col_tmp.setAttribute("class","col-sm-2")
		div_col_tmp.setAttribute("style","width:20%")
		var span = document.createElement("span");
		var msg = "In the worst case, every entry at Levels 1 to L-1 is an update to an existing entry at Level L. " +
		"Since the fraction of entries at Levels 1 to L-1 is 1/T of the overall number of entries, and at Level L, at most Z-1 of the runs " +
		"may be completely filled with obsolete entries. The space-amplification can be modeled as Z - 1 + 1/T"
		span.setAttribute("data-tooltip",msg);
		span.setAttribute("data-tooltip-position","bottom")
		var p=document.createElement("p");
		p.setAttribute("style","text-align: center;font-size:15px;font-weight:bold")
		p.textContent=("Space-Amp")
		span.appendChild(p);
		div_col_tmp.appendChild(span);
		div_schema_row.appendChild(div_col_tmp);


	    var second_text_array = [
				"Non-result",
				"Existing",
				"Short",
				"Long",
				"",
				""
			]
			var second_msg_array = [
				"The cost of non-result point lookups can be modeled as the sum of false positive rates across every run's Bloom filters",
				"The worst-case point lookup cost to an existing entry occurs when the target key is to the oldest run at the largest level. The expected I/O cost is one I/O to this target run plus the sum of FPRs across all other runs",
				"A short range lookup issues at most K I/Os per level to the smaller L-1 Levels and at most Z I/Os to the largest level for a total of Z+K*(L-1) random I/Os",
				"A long range lookup continues with a sequential scan to the relevant key range at each run issuing at least s/B sequential I/Os. The number of sequential I/Os is amplified by (1+1/T) for updated entries at Levels 1 to L-1 and Z for updated entries at Level L.",
				"",
				""
			]

			var class_array = [
				"col-sm-2",
				"col-sm-2",
				"col-sm-2",
				"col-sm-2",
				"col-sm-2",
				"col-sm-2"
			]
			var style_array= [
				"width:20%",
				"width:13.333%",
				"",
				"",
				"width:13.33%",
				"width:20%",
			]
			var text_array = [
				"Zero-result Point Lookup",
				"Existing Point Lookup",
				"Short Range Lookup",
				"Long Range Lookup",
				"Update",
				"Space Amplification"
			]
			var sum=r+qL+qS+v+w;
			var coefficient_array = [
				r/sum,
				v/sum,
				qS/sum,
				qL/sum,
				w/sum,
				0
			]
			var leveled_function_array = [
				getLeveledNonExistingPointLookupCost,
				getLeveledExistingPointLookupCost,
				getLeveledShortRangeLookupCost,
				getLeveledLongRangeLookupCost,
				getLeveledUpdateCost,
				getLeveledSpaceAmp
			]
			var total_function_array = [
				getTotalNonExistingPointLookupCost,
				getTotalExistingPointLookupCost,
				getTotalShortRangeLookupCost,
				getTotalLongRangeLookupCost,
				getTotalUpdateCost,
				getTotalSpaceAmp
			]

	result_div.appendChild(div_new_row);

	var div_new_row=document.createElement("div");
	div_new_row.setAttribute("class","row")
	var div_level = document.createElement("div");
	div_level.setAttribute("class","col-sm-1");
	div_level.setAttribute("style","width:5.33%");
	div_new_row.append(div_level);

	var div_col2 = document.createElement("div");
	div_col2.setAttribute("class","col-sm-6");
	div_col2.setAttribute("style","width:51%");
	var div_col2_row = document.createElement("div");
	div_col2_row.setAttribute("class","row");
	div_col2.append(div_col2_row);
	div_new_row.append(div_col2);

	for(i=0;i <= 5;i++){
		var div_col_tmp=document.createElement("div");
		div_col_tmp.setAttribute("class",class_array[i]);
		div_col_tmp.setAttribute("style",style_array[i]);
		var span = document.createElement("span");
		if(second_msg_array[i] != ""){
			span.setAttribute("data-tooltip",second_msg_array[i]);
			span.setAttribute("data-tooltip-position","bottom")
		}

		var p = document.createElement("p");
		p.setAttribute("style","text-align: center;font-size:15px")
		p.textContent=(second_text_array[i])
		span.appendChild(p)
		div_col_tmp.appendChild(span);
		div_col2_row.appendChild(div_col_tmp);
	}
 result_div.appendChild(div_new_row);



//adding the buffer row
    var div_new_row=document.createElement("div");
    div_new_row.setAttribute("class","row")


				var div_level = document.createElement("div");
				div_level.setAttribute("class","col-sm-1");
				div_level.setAttribute("style","width:5.33%");
				var p=document.createElement("p");
				p.setAttribute("style","text-align: center;font-size:15px")
				p.textContent=("Buffer")
				div_level.appendChild(p);

				var div_col2 = document.createElement("div");
				div_col2.setAttribute("class","col-sm-6");
				div_col2.setAttribute("style","width:51%");
				var div_col2_row = document.createElement("div");
				div_col2_row.setAttribute("class","row");
				div_col2.appendChild(div_col2_row);

				var div_col_tmp=document.createElement("div");
				div_col_tmp.setAttribute("class","col-sm-2")
				div_col_tmp.setAttribute("style","width:20%")
				var p_tmp=document.createElement("p");
				p_tmp.setAttribute("style","text-align: center;font-size:15px")
				p_tmp.textContent=("-");
				div_col_tmp.appendChild(p_tmp);
				div_col2_row.appendChild(div_col_tmp);


				for(i=1;i<=5;i++){
					var div_col_tmp=document.createElement("div");
	        div_col_tmp.setAttribute("class",class_array[i])
					div_col_tmp.setAttribute("style",style_array[i])
	        var p_tmp=document.createElement("p");
	        p_tmp.setAttribute("style","text-align: center;font-size:15px")
	        p_tmp.textContent=("-");
	        div_col_tmp.appendChild(p_tmp);
					div_col2_row.appendChild(div_col_tmp);
				}


        var div_col3=document.createElement("div");
        div_col3.setAttribute("class","col-sm-5")
        div_col3.setAttribute("style","text-align: center;")
        var button=document.createElement("button");
        button.setAttribute("class","lsm_button lsm_button_buffer");
        var text=document.createTextNode(numberWithCommas(Math.floor(mbuffer/E)))
        button.appendChild(text);
        div_col3.appendChild(button);


				div_new_row.appendChild(div_level);
        div_new_row.appendChild(div_col2);
        div_new_row.appendChild(div_col3);
    result_div.appendChild(div_new_row);


    var isMobile;
    if (window.matchMedia)
    {
        isMobile = window.matchMedia('(max-device-width: 768px)').matches;
    }
    else
    {
        isMobile = screen.width <= 768;
    }

    if (window.matchMedia)
    {
        isMobile = window.matchMedia('(max-device-width: 768px)').matches;
    }
    else
    {
        isMobile = screen.width <= 768;
    }

    console.log(screen.width)

    var max_button_size=457;
    if (screen.width<=1200)
    {
        max_button_size=Math.max(screen.width-700,350);
    }



    // if (isMobile)
    //     max_button_size=350;

    var lsm_button_size_ratio=(max_button_size-100)/filters.length;
    var cur_length=100;
    cur_length+=lsm_button_size_ratio;

    if (N<=(mbuffer/E))
    {
    	//nothing in LSM tree
//adding the buffer row
    var div_new_row=document.createElement("div");
    div_new_row.setAttribute("class","row");;
		var div_level = document.createElement("div");
		div_level.setAttribute("class","col-sm-1");
		div_level.setAttribute("style","width:5.33%");
		var p=document.createElement("p");
		p.setAttribute("style","text-align: center;font-size:15px")
		p.textContent=("Buffer")
		div_level.appendChild(p);
		var div_col2 = document.createElement("div");
		div_col2.setAttribute("class", "col-sm-6");
		div_col2.setAttribute("style","width:51%");
		var div_col2_row = document.createElement("div");
		div_col2_row.setAttribute("class","row");
		div_col2.appendChild(div_col2_row);

		for(i=0;i<=5;i++){
			var div_col_tmp=document.createElement("div");
			div_col_tmp.setAttribute("class",class_array[i])
			div_col_tmp.setAttribute("style",style_array[i])
			var p_tmp=document.createElement("p");
			p_tmp.setAttribute("style","text-align: center;font-size:15px")
			p_tmp.textContent=("-");
			div_col_tmp.appendChild(p_tmp);
			div_col2_row.appendChild(div_col_tmp);
		}

        var div_col3=document.createElement("div");
        div_col3.setAttribute("class","col-sm-5")
        div_col3.setAttribute("style","text-align: center;")
        var p4=document.createElement("p");
        p4.setAttribute("style","text-align: center;")
        p4.textContent=("All data entries fit into the buffer")
        var p4b=document.createElement("p");
        p4b.setAttribute("style","text-align: center;")
        p4b.textContent=("Add more entries to see a tree!")

        div_col3.appendChild(p4);
        div_col3.appendChild(p4b);

				div_new_row.appendChild(div_level);
				div_new_row.appendChild(div_col2);
        div_new_row.appendChild(div_col3);

    result_div.appendChild(div_new_row);


	    var hr=document.createElement("hr");
	    result_div.appendChild(hr)

    }
    else
    {

		var last_is_smaller=false;
		var full_runs_in_last_level=0;
		var L = filters.length;

	    for (var i=0;i<L;i++)
	    {
		    var div_new_row=document.createElement("div");
		    div_new_row.setAttribute("class","row");

				var div_level = document.createElement("div");
				div_level.setAttribute("class","col-sm-1");
				div_level.setAttribute("style","width:5.33%");
				var p=document.createElement("p");
				p.setAttribute("style","text-align: center;font-size:15px")
				var span=document.createElement("span");
				var MSD=getMostSignificantDigit(filters[i].fp);
	        	if (MSD>10)
	        		MSD=10;
				message  = "The false positive rate for Bloom filters at Level "+(i+1)+
				" is "+(filters[i].fp*100).toFixed(MSD+1)+"%. This entails "+(filters[i].mem/filters[i].nokeys).toFixed(2)+
				" bits-per-element, resulting in a total of "+formatBytes(filters[i].mem/8,1)+" for Bloom filters at this level out of "+formatBytes(mfilter_bits/8)+" for Bloom filters across all levels."
				span.setAttribute("data-tooltip",message);
				span.setAttribute("data-tooltip-position","bottom")
				span.textContent=(i+1+"");
				p.appendChild(span);
				div_level.appendChild(p);

				var div_col2 = document.createElement("div");
				div_col2.setAttribute("class", "col-sm-6");
				div_col2.setAttribute("style","width:51%");
				var div_col2_row = document.createElement("div");
				div_col2_row.setAttribute("class","row");
				div_col2.appendChild(div_col2_row);

				for(j=0;j <= 5;j++){
					var div_col2_tmp=document.createElement("div");
	        div_col2_tmp.setAttribute("class",class_array[j]);
					div_col2_tmp.setAttribute("style",style_array[j])
					var p2_tmp=document.createElement("p");
					var span2_tmp=document.createElement("span");

					var cost = leveled_function_array[j](i+1, L, filters, N, T, B, Y, K, Z, s, Mu, isOptimalFPR)
					var threshold_flag=false;
					var message;
					var msg_cost = cost;
					if(cost*1000%1 != 0){
						msg_cost=cost.toExponential(5);
					}
					if(j != 5){
						message = text_array[j] + " at this level has " + msg_cost + " I/O cost(s)."
					}else{
						message = text_array[j] + " at this level is " + msg_cost + "."
					}
					if(cost <= THRESHOLD){
						if(cost != 0){
							threshold_flag=true;
						}
						cost = 0.0;
					}else if(typeof cost == 'number'  && cost*1000 < 1){
						cost = myFloor(cost, 1).toExponential(1)
					}else if(cost*1000%1 != 0){
						cost = (Math.floor(cost*1000)/1000).toFixed(3)
					}
					if(threshold_flag){
						message += "Because the value here is too small (less than 1e-9), it is noted as 0 in breakdown table. "
					}

					span2_tmp.setAttribute("data-tooltip",message);
	        span2_tmp.setAttribute("data-tooltip-position","bottom")
					p2_tmp.setAttribute("style","text-align: center;font-size:15px")
					p2_tmp.textContent=(cost+"")
					span2_tmp.appendChild(p2_tmp);
					div_col2_tmp.appendChild(span2_tmp);
					div_col2_row.appendChild(div_col2_tmp);
				}





		   	  var div_col3=document.createElement("div");
			    div_col3.setAttribute("class","col-sm-5")
			    div_col3.setAttribute("style","text-align: center;")

			    var levelcss=i+1;
			    if (L<5)
			    	levelcss=5-L+1+i;
                // console.log(i+":"+levelcss)
								var n;
                if (i >= filters.length-Y-1) {
									maxRuns = Z;
									n = Math.min(Z, 7);
                } else {
									maxRuns = K;
                  n = Math.min(K, 7);
								}

                    for (var j = 0; j < n; j++) {
                        if (maxRuns > 6 && j == 5) {
                            var span =document.createElement("span");
                            var message="This level contains "+maxRuns+" runs";
                            span.setAttribute("data-tooltip", message);
                            span.setAttribute("data-tooltip-position", "left");
                            span.setAttribute("style", "width:19.27px; font-size: 20px; color: #a51c30; padding: 0px 2px");
														span.id = i + "span";
                            span.textContent=" ...";
                            div_col3.appendChild(span);
                        } else {
                            var button=document.createElement("button");
			                button.setAttribute("class","lsm_button lsm_button"+(levelcss));
			                if (last_is_smaller && j >= full_runs_in_last_level)
			    	            button.setAttribute("class","lsm_button_not_solid");
			                else
			    	            button.setAttribute("class","lsm_button");
												if(maxRuns >= 7){
													button.setAttribute("style","width: "+(cur_length- 19.27)/6+"px; height: 36px");
												}else{
													button.setAttribute("style","width: "+cur_length/n+"px; height: 36px");
												}


                            var message = "At this level, each run contains "+numberWithCommas(Math.floor(filters[i].nokeys/maxRuns)+" entries"+((last_is_smaller)? " (level is not full)" : ""));
                            button.setAttribute("data-tooltip", message);
                            button.setAttribute("data-tooltip-position", "left");
                            div_col3.appendChild(button);
                        }
                    }
                    cur_length+=lsm_button_size_ratio;

				        if (i==(L-2)) {
				            if (T*filters[i].nokeys>filters[i+1].nokeys)
				    	            last_is_smaller=true;
													full_runs_in_last_level = Math.floor(filters[i+1].nokeys/(T*filters[i].nokeys/Z))
                }

								div_new_row.appendChild(div_level);
								div_new_row.appendChild(div_col2);
				        div_new_row.appendChild(div_col3);
	    	  result_div.appendChild(div_new_row);


	    }



	    var hr=document.createElement("hr");
	    result_div.appendChild(hr)


	    //total cost line
	    var div_new_row=document.createElement("div");
	    div_new_row.setAttribute("class","row")

/*
	    	var Rbaseline=eval_R(filters_baseline, leveltier, T);
	    	var Rmonkey=eval_R(filters_monkey, leveltier, T);
				var RfluidLSMTree;
				if(leveltier == 3){
					RfluidLSMTree = eval_R(filters_fluid_lsm_tree, leveltier, T, K, Z)
				}
*/
				var div_col1 = document.createElement("div");
				div_col1.setAttribute("class","col-sm-1")
				div_col1.setAttribute("style","width:5.33%");
			var p1=document.createElement("p");
			p1.setAttribute("style","text-align: center;font-size:15px;font-weight:bold")
			p1.textContent=("Total:")
			div_col1.appendChild(p1);


	        var div_col2=document.createElement("div");
		    div_col2.setAttribute("class","col-sm-6");
				div_col2.setAttribute("style","width:51%");
				var div_col2_row=document.createElement("div");
				div_col2_row.setAttribute("class","row");
				div_col2.appendChild(div_col2_row);

				var total_cost=0;

				for(j=0;j<=5;j++){
					var div_col2_tmp=document.createElement("div");
	        div_col2_tmp.setAttribute("class",class_array[j]);
					div_col2_tmp.setAttribute("style",style_array[j]);
					var p2_tmp=document.createElement("p");
					var span2_tmp=document.createElement("span");

					var cost = total_function_array[j](i+1, L, filters, N, T, B, Y, K, Z, s, Mu, isOptimalFPR)
					total_cost += coefficient_array[j]*cost;
					var msg_cost = cost;
					if(cost*1000%1 != 0){
						msg_cost=cost.toExponential(5);
					}
					var message;
					if(j != 5){
						message = "The total cost of " +text_array[j]+ " is " + 	msg_cost + " I/O(s)."
					}else{
						message = "The total "+text_array[j] +" is " + 	msg_cost + "."
					}
					var threshold_flag=false;
					if(cost <= THRESHOLD){
						if(cost != 0){
							threshold_flag=true;
						}
						cost = 0.0;
					}else if(typeof cost == 'number'  && cost*1000 < 1){
						cost = myCeil(cost, 1).toExponential(1)
					}else if(cost*1000%1 != 0){
						cost = (Math.ceil(cost*1000)/1000).toFixed(3)
					}
					if(threshold_flag){
						message += "Because the value here is too small (less than 1e-9), it is noted as 0 in breakdown table. "
					}
					p2_tmp.textContent=(cost+"")
					span2_tmp.setAttribute("data-tooltip",message);
	        span2_tmp.setAttribute("data-tooltip-position","bottom")
					p2_tmp.setAttribute("style","text-align: center;font-size:15px;font-weight:bold")
					span2_tmp.appendChild(p2_tmp);
					div_col2_tmp.appendChild(span2_tmp);
					div_col2_row.appendChild(div_col2_tmp);
				}


					var div_col3=document.createElement("div");
		    div_col3.setAttribute("class","col-sm-5");
				//div_col3.setAttribute("style","width:40%");
				var p3_tmp=document.createElement("p");
				var span3_tmp=document.createElement("span");
				if(total_cost*1000%1 != 0){
					message="Total cost under the specified workload is " + total_cost.toExponential(5) + " I/Os"
				}else{
					message="Total cost under the specified workload is " + total_cost + " I/Os"
				}
				p3_tmp.textContent=message;
				span3_tmp.setAttribute("data-tooltip",message);
				span3_tmp.setAttribute("data-tooltip-position","bottom")
				p3_tmp.setAttribute("style","text-align: center;font-size:15px;font-weight:bold")
				span3_tmp.appendChild(p3_tmp);
				div_col3.appendChild(span3_tmp)
				div_new_row.appendChild(div_col1);
				div_new_row.appendChild(div_col2);
				div_new_row.appendChild(div_col3);


	        /*if (leveltier==0)
	            message+=" TIERING: When a run is flushed to the next level it is NOT merged with the already existing runs (if any). Hence, every time a run is pushed to the next level it is written only once.";
	        else if (leveltier==1)
	            message+=" LEVELING: When a run is flushed to the next level it is ALWAYS merged with the already existing runs (if any). Hence, every flush causes (up to) size_ratio-1 merges in the next level.";
            */


		    result_div.appendChild(div_new_row);

    // var hr=document.createElement("hr");
    // result_div.appendChild(hr)
	}

}

function getTotalCost(input_conf, L_decimal_flag=false){
	var N=input_conf.N;
	var M=input_conf.mfilter_per_entry*N/8;
	var T=input_conf.T;
	var E=input_conf.E;
	var mbuffer=input_conf.mbuffer;
	var B=input_conf.P/E;
	var P=mbuffer/input_conf.P;
	var K = input_conf.fluidK;
	var Z = input_conf.fluidZ;
	var s = input_conf.s;
	var Mu = input_conf.Mu;
	var leveltier=input_conf.leveltier;
	if(leveltier==0){
		K = T - 1;
		Z = T - 1;
	}else if(leveltier == 1){
		K = 1;
		Z = 1;
	}else if(leveltier == 2){
		K = T - 1;
		Z = 1;
	}

	var w=input_conf.w;
	var v=input_conf.v;
	var r=input_conf.r;
	var qL=input_conf.qL;
	var qS=input_conf.qS;
	var sum=w+v+r+qL+qS;
	var L;
	if(L_decimal_flag){
		L = Math.log(N*E*(T - 1)/(mbuffer*T)+ 1/T)/Math.log(T);
	}else{
		L = Math.ceil(Math.log(N*E*(T - 1)/(mbuffer*T)+ 1/T)/Math.log(T));
	}

	var Y = calc_Y(M/N, K, Z, T, L);
	var R = get_accurate_R(M, T, N, K, Z, B, P, input_conf.leveltier);
	var fpr_Z = (R - Y*Z)*(T - 1)/(T - Math.pow(T, Y - L + 1))/Z;
	var V;
	if(Y > 0){
		V = (R - Z) + Math.floor(Z/2) + 1;
	}else{
		V = (R - fpr_Z*Z) + fpr_Z*Math.floor(Z/2) + 1;
	}
	var W = ((T - 1)*((L - Y - 1)/(K + 1) + (Y + 1.0)/(Z + 1)))/Mu/B;
	var Q = qS*((Y+1)*Z + K*(L - Y - 1));
	for(j = 1; j < L-Y; j++){
		Q += qL*Math.ceil(s/B/Math.pow(T, L - j)/Mu);
	}
	Q += qL*Math.ceil(Z*s/B/Mu)
	var total = (w*W + r*R + v*V + Q)/sum;
	return total;


}


function AutoTune(tuneFlag, min_value, max_value, input_conf) {

	var delta_value = Math.ceil((max_value - min_value)/2);
	var tmp, tmp1, tmp2;
	var cost, cost1, cost2;
	var optimal_conf, optimal_conf1, optimal_conf2;
	var unchanged_flag = 4;
	do{
		if(unchanged_flag != 0){
			tmp = Math.ceil(min_value + delta_value);
			if(tuneFlag == 0){
				input_conf.T = tmp;
				optimal_conf = AutoTune(1, 1, tmp-1, input_conf);

			}else if(tuneFlag == 1){

				input_conf.K = tmp;
				optimal_conf = AutoTune(2, 1, input_conf.T-1, input_conf);

			}else{
				input_conf.Z = tmp;
				optimal_conf = [];
				optimal_conf.push(getTotalCost(input_conf));
				optimal_conf.push([tmp]);
			}

			cost = optimal_conf[0];
		}

		if(unchanged_flag != 1 && unchanged_flag != 2){

			if(tuneFlag == 0){
				tmp1 = Math.max(Math.floor(tmp - delta_value), min_value);
				input_conf.T = tmp1;
				optimal_conf1 = AutoTune(1, 1, tmp1-1, input_conf);
				cost1 = optimal_conf1[0];

				tmp2 = Math.min(Math.ceil(tmp + delta_value), max_value);
				input_conf.T = tmp2;
				optimal_conf2 = AutoTune(1, 1, tmp2-1, input_conf);
				cost2 = optimal_conf2[0];
			}else if(tuneFlag == 1){
				tmp1 = Math.max(Math.floor(tmp - delta_value), min_value);
				input_conf.K = tmp1;
				optimal_conf1 = AutoTune(2, 1, input_conf.T-1, input_conf);
				cost1 = optimal_conf1[0];

				tmp2 = Math.min(Math.ceil(tmp + delta_value), max_value);
				input_conf.K = tmp2;
				optimal_conf2 = AutoTune(2, 1, input_conf.T-1, input_conf);

				cost2 = optimal_conf2[0];
			}else{
				tmp1 = Math.max(Math.floor(tmp - delta_value), min_value);
				input_conf.Z = tmp1;
				optimal_conf1 = []
				optimal_conf1.push(getTotalCost(input_conf));
				optimal_conf1.push([tmp1]);
				cost1 = optimal_conf1[0];

				tmp2 = Math.min(Math.ceil(tmp + delta_value), max_value);
				input_conf.Z = tmp2;
				optimal_conf2 = []
				optimal_conf2.push(getTotalCost(input_conf));
				optimal_conf2.push([tmp2]);
				cost2 = optimal_conf2[0];
			}

		}
		if(cost1 < cost2 && cost1 < cost){
			tmp = tmp1;
			cost = cost1;
			optimal_conf = optimal_conf1;
			unchanged_flag = 1;
		}else if(cost2 < cost){
			tmp = tmp2;
			cost = cost2;
			optimal_conf = optimal_conf2;
			unchanged_flag = 2;
		}else{
			unchanged_flag = 0;
		}
		delta_value = delta_value/2;

	}while(delta_value >= 1)
	if(tuneFlag != 2){
		optimal_conf[1].push(tmp)
	}
	console.log(optimal_conf);
	return optimal_conf
}

function AutoTuneT_test(min_value, max_value, input_conf, start_value=NaN, steps=20){
	var tmp, tmp1, tmp2;
	tmp_min = min_value;
	tmp_max = max_value;
	var cost, cost1, cost;
	var optimal_conf;
	if(isNaN(start_value)){
			tmp = Math.ceil((tmp_max - tmp_min)/2) + min_value;
	}else{
		tmp = start_value;
	}

	do{
		input_conf.T = tmp;
		cost = getTotalCost(input_conf, true);


		if(tmp <= min_value){
			cost1 = Number.MAX_VALUE;
		}else{
			input_conf.T = tmp - 1;
			cost1 = getTotalCost(input_conf, true);
		}


		if(tmp >= max_value){
			cost2 = Number.MAX_VALUE;
		}else{
			input_conf.T = tmp + 1;
			cost2 = getTotalCost(input_conf, true);
		}

		if(cost < cost1 && cost < cost2){
			break;
		}
		var delta = (cost2 - cost1)/2*steps;
		if(delta > 0){
			tmp = tmp - Math.ceil(delta);
		}else{
			tmp = tmp - Math.floor(delta)
		}
	}while(tmp < tmp_max && tmp > tmp_min);

	optimal_conf = [];
	optimal_conf.push(cost);
	optimal_conf.push([tmp]);

	return optimal_conf;
}

function AutoTuneT(min_value, max_value, input_conf) {

	var delta_value = Math.ceil((max_value - min_value)/2);
	var tmp, tmp1, tmp2;
	var cost, cost1, cost2;
	var optimal_conf, optimal_conf1, optimal_conf2;
	var unchanged_flag = 4;
	do{
		if(unchanged_flag != 0){
			tmp = min_value + delta_value;
			input_conf.T = tmp;
			optimal_conf = [];
			optimal_conf.push(getTotalCost(input_conf));
			optimal_conf.push([tmp]);
			cost = optimal_conf[0];
		}

		if(unchanged_flag != 1 && unchanged_flag != 2){

			tmp1 = Math.max(Math.floor(tmp - delta_value), min_value);
			input_conf.T = tmp1;
			optimal_conf1 = []
			optimal_conf1.push(getTotalCost(input_conf));
			optimal_conf1.push([tmp1]);
			cost1 = optimal_conf1[0];

			tmp2 = Math.min(Math.ceil(tmp + delta_value), max_value);
			input_conf.T = tmp2;
			optimal_conf2 = []
			optimal_conf2.push(getTotalCost(input_conf));
			optimal_conf2.push([tmp2]);
			cost2 = optimal_conf2[0];

	}
		if(cost1 < cost2 && cost1 < cost){
			max_value = tmp;
			tmp = tmp1;
			cost = cost1;
			optimal_conf = optimal_conf1;
			unchanged_flag = 1;
		}else if(cost2 < cost){
			min_value = tmp;
			tmp = tmp2;
			cost = cost2;
			optimal_conf = optimal_conf2;
			unchanged_flag = 2;
		}else{
			max_value = tmp2;
			min_value = tmp1;
			unchanged_flag = 0;
		}
		delta_value = delta_value/2;

	}while(delta_value >= 1)
	console.log(optimal_conf)
	return optimal_conf
}

function calc_T(N, mbuffer, E, L){
	var tmp = N/(mbuffer/E);
	var test_T = Math.ceil(Math.pow(tmp/(L+1), 1.0/L));
	var test_T1 = test_T;
	var test_T2 = test_T;
	while((Math.pow(test_T1, L+1) - 1)/(test_T1-1) < tmp){
		test_T1 -= 1
	}
	while((Math.pow(test_T2, L+1) - 1)/(test_T2-1) < tmp){
		test_T2 += 1
	}
	if(test_T1 >= 2 && test_T1 != test_T2){
		return [test_T1, test_T2];
	}else{
		return [test_T2];
	}
}

function AutoTuneL(input_conf, tree_types_array){
	var optimal_conf = []
	var cost = Number.MAX_VALUE;
	var Lmax = Math.ceil(Math.log(input_conf.N/(2*input_conf.mbuffer/input_conf.E)+ 1/2)/Math.log(2));
	for(var L=1;L<=Lmax;L++){
		T_array = calc_T(input_conf.N, input_conf.mbuffer, input_conf.E, L);
		for(var i=0;i<T_array.length;i++){
			input_conf.T = T_array[i];
			for(var j=0;j<tree_types_array.length;j++){
				input_conf.leveltier=tree_types_array[j];
				var tmp_cost=getTotalCost(input_conf);
				if(tmp_cost < cost){
					cost = tmp_cost;
					optimal_conf = []
					optimal_conf.push(cost);
					optimal_conf.push([tree_types_array[j], L, T_array[i]])
				}
			}


		}
	}
	console.log(optimal_conf);
	return optimal_conf;
}
//////////////////////////////////////
//AUX HOLISTIC TUNING FUNCTIONS

function LSM_config() {
    var P;
    var T;
    var L;
    var N;
    var R;
    var W;
    var M;
    var B;
    var E;
		var K;
		var Z;
		var Y;
		var Mu;
		var mfilter;
		var leveltier;
    var valid;
    var throughput;
    var num_levels;
}

// returns write cost
function get_W(M, N, T, K, Z, B, P, Mu, leveltier) {
	if(leveltier==0){
		K = T - 1;
		Z = T - 1;
	}else if(leveltier == 1){
		K = 1;
		Z = 1;
	}else if(leveltier == 2){
		K = T - 1;
		Z = 1;
	}

    if (P * B >= N) {
        return 0;
    }
		M = M*8;
		var L = Math.log(N*(T - 1)/(B*P*T)+ 1/T)/Math.log(T);
		var Y = calc_Y(M/N, K, Z, T, L);
    var result = (T - 1)*((L - Y - 1)/(K + 1) + (Y + 1.0)/(Z + 1))/(Mu*B);
    return result;
}

// read cost under state-of-the-art approach
function get_R_uniform_strategy(M, T, N, K, Z, B, P, leveltier) {
	  if(leveltier==0){
			K = T - 1;
			Z = T - 1;
		}else if(leveltier == 1){
			K = 1;
			Z = 1;
		}else if(leveltier == 2){
			K = T - 1;
			Z = 1;
		}
    M = M * 8;

		var L = Math.log(N*(T - 1) / (B * P * T)) / Math.log(T);

    if (P * B >= N) {
        return 0;
    }

    if (M <= 0) {
        return K*(L - 1) + Z;
    }

    //double L = !leveled ? floor(1 +  log(N/(B * P)) / log(T)) : ceil(1 +  log(N/(B * P)) / log(T));
    //double L = !leveled ? floor(1 +  log(N/(B * P)) / log(T)) : ceil(1 +  log(N/(B * P)) / log(T));

    var T_part = 1/( 1 - Math.pow(T, -L) );

    var exponent = (M / N) * Math.log(2) * Math.log(2) * T_part;

    var EULER = 2.71828182845904523536;
    var bottom = Math.pow(EULER, exponent);
    var R = 1.0 / bottom;

    return R*((L - 1)*K + Z);
}


// read cost under MonKey
function get_accurate_R(M, T, N, K, Z, B, P, leveltier) {

	if(leveltier==0){
		K = T - 1;
		Z = T - 1;
	}else if(leveltier == 1){
		K = 1;
		Z = 1;
	}else if(leveltier == 2){
		K = T - 1;
		Z = 1;
	}

	M = M*8;


	var L = Math.log(N*(T - 1)/(B*P*T)+ 1/T)/Math.log(T);
	var Y = calc_Y(M/N, K, Z, T, L);

	if (P * B >= N) {
			return 0;
	}

	if (M <= 0) {
			return K*(L - 1) + Z;
	}

	if(Y <= 0){
		Y = 0;
	}

	var tmp = 1/(Math.pow(T, 1-Y) - Math.pow(T, 1-L));
	var exponent = ((Math.pow(T, 1-Y) - Math.pow(T, 2-L))/(T-1) - (L-Y-1)/Math.pow(T, L-1))*tmp;
	var result = ((T - Math.pow(T, Y-L+1))/(T-1))*K*Math.pow(T, exponent)*Math.pow(Z*1.0/K, tmp*(T - 1.0)/Math.pow(T, Y))*Math.exp(-(M/N)*Math.log(2)*Math.log(2)*tmp*(T - Math.pow(T, -L)))+Y*Z


  return result;
}

//////////////////////////////////////




function smooth(to_print) {
    var num_elimintated = 0;
    for (var i = 2; i < to_print.length; i++) {
        var c1 = to_print[i-2];
        var c2 = to_print[i-1];
        var c3 = to_print[i];
        var rise1 = c2.R - c1.R;
        var run1 = c2.W - c1.W;
        var gradient1 = Math.abs(rise1 / run1);

        var rise2 = c3.R - c1.R;
        var run2 = c3.W - c1.W;
        var gradient2 = Math.abs(rise2 / run2);

        if (gradient1 < gradient2) {
            // to_print.erase(to_print.begin() + i - 1);
            to_print.splice(i - 1,1);
            num_elimintated++;
        }
    }
    return num_elimintated;
}

function print_csv_line(c, num_commas, print_details, differentiate_tiered_leveled) {
    // printf("%.7f, ", c.W);
    var message="";
    message+=((c.W).toFixed(4));

    message+=(", ");



    // if (c.L == 1 && differentiate_tiered_leveled) {
    //  // printf(", ");
    //  message+=(", ");
    // }

    // for (var j = 0; j < num_commas; j++) {
    //  // printf(", ");
    //  message+=(", ");
    //  if (differentiate_tiered_leveled) {
    //      message+=(", ");
    //      // printf(", ");
    //  }
    // }

    // printf("%.4f", c.R);
    message+=" ";
    message+=(  pad((c.R).toFixed(4),7," ")  );

    if (print_details) {
        // printf(",\t%g,\t%g,\t%g,\t%g", c.T, c.P, c.L, c.num_levels);
        var TL=((c.L==0)?"tier":"level");
        var memory_for_level_0=(c.P*c.B*c.E);
        var memory_for_BF=c.M-memory_for_level_0;
        var level_0_size_MB=memory_for_level_0/1024/1024;
        var BF_size_MB=memory_for_BF/1024/1024;
        message+=(",\t"+pad(c.T,4," ")+",\t"+pad(level_0_size_MB.toFixed(2),7," ")+",\t"+pad(TL,5," ")+",\t"+c.num_levels+",\t\t"+(BF_size_MB).toFixed(2));
    }
    console.log(message);
    //document.getElementById("explanationHolistic").value+=(message+"\n");
    // printf("\n");
}




// Iterate through leveling/tiering, size ratio T, and buffer size P to find optimal parameters
function find_optimal_R(input_conf, constant_buffer_size = -1, isOptimalFPR = true) {
    var conf = new LSM_config();
    conf.P=input_conf.P;
    conf.T=input_conf.T;
    conf.L=input_conf.L;
    conf.N=input_conf.N;
    conf.R=input_conf.R;
    conf.W=input_conf.W;
    conf.M=input_conf.M;
    conf.B=input_conf.B;
    conf.E=input_conf.E;
		conf.K=input_conf.K;
		conf.Z=input_conf.Z;
		conf.Y=input_conf.Y;
		conf.mfilter=input_conf.mfilter;
		conf.leveltier=input_conf.leveltier;
		conf.Mu=input_conf.Mu;
    conf.valid=input_conf.valid;
    conf.throughput=input_conf.throughput;
    conf.num_levels=input_conf.num_levels;

    var N = conf.N;
    var B = conf.B;
    var E = conf.E;
    var M = conf.M;
    var W = conf.W;
		var K = conf.K;
		var Z = conf.Z;
		var Y = conf.Y;
		var Mu = conf.Mu;
    conf.valid = false;

    var max_elements_in_buffer = N;

    var starting_buffer_size = 1;
    if (constant_buffer_size > -1) {
        starting_buffer_size = constant_buffer_size / (E * B);
        max_elements_in_buffer = starting_buffer_size * B;
    }

    var min_W = 1;
    var min_R = (1 + Math.log(N/2) / Math.log(2)) * B, best_T = 2, best_P = 1;
    var is_leveled = 1;
    for (var T = 2; T < B * 4; T++) {
        for (var P = starting_buffer_size;  P * B <= max_elements_in_buffer; P *= T) {
            for (var L = 0; L <= 1; L++) {
                var leveled = L;
                var num_levels =  Math.ceil(Math.log(N*(T - 1)/ (B * P * T)) / Math.log(T));
                var mem_for_bloom = M - P * B * E;

                var current_R = isOptimalFPR ? get_accurate_R(mem_for_bloom, T, N, K, Z, B, P, L) :
                        get_R_uniform_strategy(mem_for_bloom, T, N, K, Z, B, P, L);

                var current_W = get_W(mem_for_bloom, N, T, K, Z, B, P, Mu, L);
                if (mem_for_bloom >= 0 && current_R < min_R && current_W <= W) {
                    min_R = current_R;
                    best_T = T;
                    best_P = P;
                    min_W = current_W;
                    is_leveled = leveled;
                    conf.W = current_W;
                    conf.valid = true;
                    conf.num_levels = 1 +  Math.log(N*(T - 1)/(B * P * T)) / Math.log(T);
                }
            }
        }
    }
    conf.T = best_T;
    conf.P = best_P;
    conf.leveltier = is_leveled;
    conf.R = min_R;
    return conf;
}




function getPoint(leveltier, T, mfilter, conf, isOptimalFPR) {
                var meetingR;
                if (isOptimalFPR) {
                    meetingR = get_accurate_R(mfilter, T, conf.N, conf.K, conf.Z, conf.B, conf.P, leveltier);
                }
                else {
                    meetingR = get_R_uniform_strategy(mfilter, T, conf.N, conf.K, conf.Z, conf.B, conf.P, leveltier);
                }
                var meetingW = get_W(conf.M-conf.B*conf.P*conf.E, conf.N, T, conf.K, conf.Z, conf.B, conf.P, conf.Mu, leveltier);
                return {W: meetingW, R: meetingR};
}



// Find the corresponding optimal read cost for different upper bounds on write cost
function print_csv_experiment(input_conf, num_commas, print_details, fix_buffer_size = -1, isOptimalFPR = true, smoothing = false, differentiate_tiered_leveled = true) {

    var conf = new LSM_config();
    conf.P=input_conf.P;
    conf.T=input_conf.T;
    conf.L=input_conf.L;
    conf.N=input_conf.N;
    conf.R=input_conf.R;
    conf.W=input_conf.W;
    conf.M=input_conf.M;
    conf.B=input_conf.B;
    conf.E=input_conf.E;
		conf.K=input_conf.K;
		conf.Z=input_conf.Z;
		conf.Y=input_conf.Y;
		conf.mfilter=input_conf.mfilter;
		conf.Mu=input_conf.Mu;

    conf.valid=input_conf.valid;
    conf.throughput=input_conf.throughput;
    conf.num_levels=input_conf.num_levels;

    var prev_W = -1;
    var prev_L = 0;

    var array_to_print = new Array();

    var meeting_W = getPoint(0, 2, conf.mfilter, conf, isOptimalFPR).W;
    var min_W = getPoint(0, 4, conf.mfilter, conf, isOptimalFPR).W;
    var max_W = getPoint(1, 10, conf.mfilter, conf, isOptimalFPR).W;

    var starting_W = min_W * 0.1;
    var finishing_W = max_W * 10;
    var steps_for_leveling = (finishing_W - meeting_W) / 200;

    //print_write_optimized_extreme(conf, use_new_strategy);
    // console.log(input_conf)
    // console.log(conf)

    var last_c = LSM_config();
    for (var W = starting_W; W <= finishing_W; W += 0.001) {
        conf.W = W;

        var c = find_optimal_R(conf, fix_buffer_size, isOptimalFPR);
        // console.log(c)
        // print_csv_line(c,num_commas,print_details,differentiate_tiered_leveled);

        // document.getElementById("explanation").value += (c.W.toFixed(4) + "\t\t" + c.R.toFixed(2) + "\t\t" + c.L + "\t\t" + c.T + "\t\t" + c.valid + "\n");

        // return;

        if ((prev_W == c.W && prev_L == c.L) || !c.valid) {
            continue;
        }
        prev_W = c.W;
        prev_L = c.leveltier;

        /*if (c.L == 1 && last_c.L == 0 && last_c.T == 2 && differentiate_tiered_leveled) {
            last_c.L = 1;
            array_to_print.push(last_c);
            //print_csv_line(last_c, num_commas, print_details);
        }*/

        array_to_print.push(c);
        //print_csv_line(c, num_commas, print_details);

        last_c = c;

        if (W >= meeting_W) {
            W += steps_for_leveling;
        }
    }

    var num_eliminated = 1;
    while (smoothing && num_eliminated > 0) {
        num_eliminated = smooth(array_to_print);
        //printf("%d   \n", num_eliminated);
    }

    // console.log(array_to_print)
    return (array_to_print);

    // console.log("printing "+array_to_print.length+" rows")

    // document.getElementById("explanation").value+=("Printing the "+array_to_print.length+ " configurations on the performance skyline.\n");
    // document.getElementById("explanation").value+=("write,\tread,\tT,\tL0 (MB),  merge,  levels,\tBF (MB)\n");

    // document.getElementById("explanation").value+="array size " + array_to_print.length;

    // for (var i = 0; i < array_to_print.length; i++) {
    //     print_csv_line(array_to_print[i], num_commas, print_details, differentiate_tiered_leveled);
    // }

    // printf("\n");
}

function AutoTune1(e){
	var inputParameters = parseInputTextBoxes();
	var optimal_conf = AutoTuneL(inputParameters, [inputParameters.leveltier]);
	var startT = optimal_conf[1][optimal_conf[1].length-1];
	var Tmax = 4* inputParameters.P/inputParameters.E;
	var tmp_optimal_conf=AutoTuneT_test(2, Tmax, inputParameters, startT);
	console.log(tmp_optimal_conf);
	if(tmp_optimal_conf[0] < optimal_conf[0]){
		optimal_conf = tmp_optimal_conf;
	}
	document.getElementById("T").value=optimal_conf[1][optimal_conf[1].length-1];
	e.target.id="T";
	re_run(e);
	clickbloomTuningButton(true)
}

function AutoTune2(e){
	var inputParameters = parseInputTextBoxes();
	var optimal_conf = AutoTuneL(inputParameters, [0, 1, 2]);
	//var Tmax = 4* inputParameters.P/inputParameters.E;
	//AutoTuneT(2, Tmax, inputParameters);
	var tree_type=optimal_conf[1][0];
	var startT = optimal_conf[1][optimal_conf[1].length-1];
	var Tmax = 4* inputParameters.P/inputParameters.E;
	var tmp_optimal_conf=AutoTuneT_test(2, Tmax, inputParameters, startT);
	console.log(tmp_optimal_conf);
	if(tmp_optimal_conf[0] < optimal_conf[0]){
		optimal_conf = tmp_optimal_conf;
	}
	document.getElementById("T").value=optimal_conf[1][optimal_conf[1].length-1];
	var radios = document.getElementsByName("ltradio");
	for(var i = 0; i < radios.length; i++){
		if(radios[i].value==tree_type){
			radios[i].checked=true;
		}else{
			radios[i].checked=false;
		}

	}
	e.target.id="T";
	re_run(e);
	clickbloomTuningButton(true)
}
function AutoTune3(e){
	var inputParameters = parseInputTextBoxes();
	var Tmax=4* inputParameters.P/inputParameters.E;
	var optimal_conf = AutoTune(0, 2, Tmax, inputParameters);
	document.getElementById("T").value=optimal_conf[1][optimal_conf[1].length-1];
	document.getElementById("Fluid LSM-Tree K").value=optimal_conf[1][1];
	document.getElementById("Fluid LSM-Tree Z").value=optimal_conf[1][0];
	e.target.id="T";
	re_run(e);
	clickbloomTuningButton(true)
}

function MergeByFliudLSMTree(){

	document.getElementById("input-group-Fluid-K").style.display = '';
	document.getElementById("input-group-Fluid-Z").style.display = '';
	var T = document.getElementById("T").value;
	if(lastTreeType == 1){
		document.getElementById("Fluid LSM-Tree K").value = 1;
		document.getElementById("Fluid LSM-Tree Z").value = 1;
	}else if(lastTreeType == 0){
		document.getElementById("Fluid LSM-Tree K").value = T - 1;
		document.getElementById("Fluid LSM-Tree Z").value = T - 1;
	}else if(lastTreeType == 2){
		document.getElementById("Fluid LSM-Tree K").value = T - 1;
		document.getElementById("Fluid LSM-Tree Z").value = 1;
	}
	re_run(event, 'input7');
}

function MergeNotyFliudLSMTree(){
	document.getElementById("input-group-Fluid-K").style.display = 'none';
	document.getElementById("input-group-Fluid-Z").style.display = 'none';
	re_run(event, 'input7');
}

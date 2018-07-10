function InputTextBoxes()
{
	var N;
	var E;
    var mbuffer;
    var T;
    var mfilter;
    var P;
    var leveltier;  // tiered is 0, leveled is 1, lazy_level is 2, fluid_LSM-Tree is 3
    var isLeveled;
		var fluidK;
		var fluidZ;
}


function parseInputTextBoxes()
{
	var parsedBoxes = new InputTextBoxes();
    parsedBoxes.N = parseInt(document.getElementById("N").value.replace(/\D/g,''),10);
    parsedBoxes.E = parseInt(document.getElementById("E").value.replace(/\D/g,''),10);
    parsedBoxes.mbuffer = parseFloat(document.getElementById("mbuffer").value.replace(/\D/g,''))*1048576;
    parsedBoxes.T = parseInt(document.getElementById("T").value.replace(/\D/g,''), 10);
    parsedBoxes.mfilter = parseFloat(document.getElementById("mfilter").value.replace(/\D/g,''))*1048576;
    parsedBoxes.P = parseInt(document.getElementById("P").value.replace(/\D/g,''), 10);
		parsedBoxes.L = parseInt(document.getElementById("L").value.replace(/\D/g,''), 10);
    parsedBoxes.isLeveled = isRadioLeveled("ltradio");  // tiered is 0, leveled is 1
    parsedBoxes.leveltier = getRadioValueByName("ltradio");
		parsedBoxes.fluidK = parseInt(document.getElementById("Fluid LSM-Tree K").value.replace(/\D/g,''), 10);
		parsedBoxes.fluidZ = parseInt(document.getElementById("Fluid LSM-Tree Z").value.replace(/\D/g,''), 10);
    return parsedBoxes;
}

function Filter() {
    var nokeys;
    var fp;
    var mem;
}

var log10 = Math.log(10);
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
    var value = Math.exp(-(f.mem / f.nokeys) * denom);
    return value;
}


// TODO for tiering, we now assume there are (T-1) runs in the last level are all equal to each other in size,
// but if the last level is not full to capacity, then there may in reality be less than (T-1) runs, and they
// may have different sizes. To fix this problem, we can insert multiple runs per filter to the filters_array
// until the number of remaining keys is 0.
function initFilters(N,E,mbuffer,T,mfilter,P,leveltier) {
     mfilter_bits=8*mfilter;

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

    //Initialize the memory per level to be equal
    for (var i=0;i<filter_array.length;i++)
    {
        filter_array[i].mem=mfilter_bits/filter_array.length;
    }
    return filter_array;
}


function getBaselineFPassigment(N,E,mbuffer,T,mfilter,P,leveltier)
{

    var filter_array = initFilters(N,E,mbuffer,T,mfilter,P,leveltier);

    // console.log(filter_array);

    var limit_on_M=mbuffer*8; //amount of memory for filters in bits
    var diff = limit_on_M;
    var change = true;
    var iteration = 0;

    while (diff > 1) {
        change = false;
        for (var i = 0; i < filter_array.length - 1; i++) {
            for (var j = i + 1; j < filter_array.length ; j++) {

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


                if (diff_new < diff_orig && filter_array[j].mem > 0 && filter_array[i].mem > 0) {
                    change = true;
                    continue;
                }
                filter_array[i].mem -= diff * 2;
                filter_array[j].mem += diff * 2;

                f1_new = calc_R(filter_array[i]);
                f2_new = calc_R(filter_array[j]);
                diff_new = Math.abs(f1_new - f2_new);

                if (diff_new < diff_orig && filter_array[j].mem > 0 && filter_array[i].mem > 0) {
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


function getMonkeyFPassigment(N, E, mbuffer, T, K, Z, mfilter, P, leveltier)
{

    var filter_array = initFilters(N,E,mbuffer,T,mfilter,P,leveltier);

    // console.log(filter_array);
    var limit_on_M=mbuffer*8; //amount of memory for filters in bits
    var diff = limit_on_M;
    var change = true;
    var iteration = 0;
    var current_R = eval_R(filter_array, leveltier, T, K, Z);
    var original = current_R;
    var value = 0;
    while (diff > 1) {
        change = false;
        for (var i = 0; i < filter_array.length - 1; i++) {
            for (var j = i + 1; j < filter_array.length ; j++) {
                filter_array[i].mem += diff;
                filter_array[j].mem -= diff;
                value = eval_R(filter_array, leveltier, T, K, Z);
                if (value < current_R && value > 0 && filter_array[j].mem > 0 ) {
                    current_R = value;
                    change = true;
                    continue;
                }
                filter_array[i].mem -= diff * 2;
                filter_array[j].mem += diff * 2;

                value = eval_R(filter_array, leveltier, T, K, Z);

                if (value < current_R && value > 0 && filter_array[i].mem > 0 ) {
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



function eval_R(filters, leveltier, T, K = 9, Z = 1)
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
    n -= 1;
    for (var i = 0; i < n ; i++)
    {
        var val = calc_R(filters[i]);
				total += val * K;
    }

    var val = calc_R(filters[filters.length - 1]);
    total += val * Z;

    return total;
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
    document.getElementById("mfilter").value=numberWithCommas(81920); //in MB (5 bits per element)
    document.getElementById("P").value=4096; //in B
    document.getElementsByName("ltradio")[0].checked=true;
    document.getElementsByName("ltradio")[1].checked=false;
		document.getElementById("Fluid LSM-Tree K").value = 1;
		document.getElementById("Fluid LSM-Tree Z").value = 1;

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
    document.getElementById("mfilter").value=numberWithCommas(81920); //in MB (10 bits per element)
    document.getElementById("P").value=4096; //in B
    document.getElementsByName("ltradio")[0].checked=false;
    document.getElementsByName("ltradio")[1].checked=true;
		document.getElementById("Fluid LSM-Tree K").value = 3;
		document.getElementById("Fluid LSM-Tree Z").value = 3;

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
    document.getElementById("mfilter").value=numberWithCommas(0); //in MB (10 bits per element)
    document.getElementById("P").value=4096; //in B
    document.getElementsByName("ltradio")[0].checked=true;
    document.getElementsByName("ltradio")[1].checked=false;
		document.getElementById("Fluid LSM-Tree K").value = 1;
		document.getElementById("Fluid LSM-Tree Z").value = 1;

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
    var mfilter=inputParameters.mfilter;
    var P=inputParameters.P;
    var leveltier=inputParameters.leveltier;
		var K = inputParameters.fluidK;
		var Z = inputParameters.fluidZ;


    if (document.getElementById("N").value=="" || document.getElementById("E").value=="" || document.getElementById("T").value=="" || document.getElementById("P").value==""
        || document.getElementById("mbuffer").value=="" || document.getElementById("mfilter").value=="" || isNaN(N)
        || isNaN(E) || isNaN(mbuffer) || isNaN(T) || isNaN(mfilter) || isNaN(P) || isNaN(leveltier) || isNaN(K) || isNaN(Z))
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

    //get BF allocation
    var filters_baseline=getBaselineFPassigment(N, E, mbuffer, T, mfilter, P, leveltier);
    var filters_monkey=getMonkeyFPassigment(N, E, mbuffer, T, T - 1, 1, mfilter, P, leveltier);
		var filters_fluid_lsm_tree;
		if(leveltier == 3){
			filters_fluid_lsm_tree=getMonkeyFPassigment(N, E, mbuffer, T, K, Z, mfilter, P, leveltier);
		}

    google.charts.setOnLoadCallback(drawChart);


	//present it nicely in the html!
    var result_div=document.getElementById("result_tuning")
	while (result_div.firstChild) {
	    result_div.removeChild(result_div.firstChild);
	}

    // var br=document.createElement("br");
    // result_div.appendChild(br);

    var div_row1=document.createElement("div");
    div_row1.setAttribute("class","row")
	    var div_title_col=document.createElement("div");
	    div_title_col.setAttribute("class","col-sm-12");
	    var lsm_header=document.createElement("h4");
			lsm_header.textContent="";
			if(leveltier != 3){
				lsm_header.textContent="Monkey vs. State-Of-The-Art";
			} else {
				lsm_header.textContent="Fluid LSM-Tree vs. Lazy Leveling vs. State-Of-The-Art";
			}
	    lsm_header.setAttribute("style","text-align: center;")
		div_title_col.appendChild(lsm_header);
		div_row1.appendChild(div_title_col);
	result_div.appendChild(div_row1); // added main row title

    var div_row2=document.createElement("div");
    div_row2.setAttribute("class","row")
	    var div_title_col1a=document.createElement("div");
	    div_title_col1a.setAttribute("class","col-sm-2")
		div_row2.appendChild(div_title_col1a);

	    var div_title_col1=document.createElement("div");
	    div_title_col1.setAttribute("class","col-sm-5")
	    var lsm_header1=document.createElement("h5");
	    lsm_header1.textContent="False Positive Rates";
	    lsm_header1.setAttribute("style","text-align: center;")
		div_title_col1.appendChild(lsm_header1);
		div_row2.appendChild(div_title_col1); // added secondary row title

	    var div_title_col2=document.createElement("div");
	    div_title_col2.setAttribute("class","col-sm-5")
	    var lsm_header2=document.createElement("h5");
	    lsm_header2.textContent="Entries per Level";
	    lsm_header2.setAttribute("style","text-align: center;")
		div_title_col2.appendChild(lsm_header2);
		div_row2.appendChild(div_title_col2); // added secondary row title
	result_div.appendChild(div_row2);


    var first_col_alignment="center";
	//titles
    var div_new_row=document.createElement("div");
    div_new_row.setAttribute("class","row");
		var div_empty = document.createElement("div");
		div_empty.setAttribute("class","col-sm-1");
		var div_schema = document.createElement("div");
		div_schema.setAttribute("class","col-sm-6");
		var div_schema_row = document.createElement("div");
		div_schema_row.setAttribute("class","row");
		div_new_row.appendChild(div_empty);
		div_new_row.appendChild(div_schema);
		div_schema.appendChild(div_schema_row);

   	    var div_col1=document.createElement("div");
	    div_col1.setAttribute("class","col-sm-2")
	    var p1=document.createElement("p");
	    p1.setAttribute("style","text-align: "+first_col_alignment+";font-size:20px")
	    p1.textContent=("Level")
	    div_col1.appendChild(p1);


   	    var div_col2=document.createElement("div");
	    div_col2.setAttribute("class","col-sm-6")
	    var p2=document.createElement("p");
	    p2.setAttribute("style","text-align: center;font-size:20px")
	    p2.textContent=("State-Of-The-Art")
	    div_col2.appendChild(p2);

   	    var div_col3=document.createElement("div");
	    div_col3.setAttribute("class","col-sm-4")
	    var p3=document.createElement("p");
	    p3.setAttribute("style","text-align: center;font-size:20px")
	    p3.textContent=("Monkey")
			if(leveltier == 3){
				p3.textContent=("Lazy Leveling")
			}
	    div_col3.appendChild(p3);




    	div_schema_row.appendChild(div_col1);
    	div_schema_row.appendChild(div_col2);
    	div_schema_row.appendChild(div_col3);

			if(leveltier == 3){
				//div_col1.setAttribute("class","col-sm-3");
				div_col2.setAttribute("class","col-sm-3");
				div_col3.setAttribute("class","col-sm-3");
				//p1.setAttribute("style","text-align: center;font-size:15px");
				//p2.setAttribute("style","text-align: center;font-size:15px");
				//p3.setAttribute("style","text-align: center;font-size:15px");
				var div_col_tmp=document.createElement("div");
	    div_col_tmp.setAttribute("class","col-sm-3");
	    var p4=document.createElement("p");
	    p4.setAttribute("style","text-align: center;");
	    p4.textContent=("Fluid LSM-Tree");
	    div_col_tmp.appendChild(p4);
			div_schema_row.appendChild(div_col_tmp);
			}






	result_div.appendChild(div_new_row);

//adding the buffer row
    var div_new_row=document.createElement("div");
    div_new_row.setAttribute("class","row")

        var div_col1=document.createElement("div");
        div_col1.setAttribute("class","col-sm-1")
        var p1=document.createElement("p");
        p1.setAttribute("style","text-align: right;")
        p1.textContent=("Buffer:")
				if(leveltier==3){
					p1.setAttribute("style","text-align: right;")
				}
        div_col1.appendChild(p1);

				var div_col2 = document.createElement("div");
				div_col2.setAttribute("class","col-sm-6");
				var div_col2_row = document.createElement("div");
				div_col2_row.setAttribute("class","row");
				div_col2.appendChild(div_col2_row);

        var div_col2a=document.createElement("div");
        div_col2a.setAttribute("class","col-sm-2")
        var p2a=document.createElement("p");
        p2a.setAttribute("style","text-align: center;")
        p2a.textContent=("0");
        div_col2a.appendChild(p2a);


        var div_col2b=document.createElement("div");
        div_col2b.setAttribute("class","col-sm-6")
        var p2b=document.createElement("p");
        p2b.setAttribute("style","text-align: center;")
        p2b.textContent=("--");
        div_col2b.appendChild(p2b);

        var div_col2c=document.createElement("div");
        div_col2c.setAttribute("class","col-sm-4")
        var p2c=document.createElement("p");
        p2c.setAttribute("style","text-align: center;")
        p2c.textContent=("--");

        div_col2c.appendChild(p2c);


				div_col2_row.appendChild(div_col2a);
				div_col2_row.appendChild(div_col2b);
				div_col2_row.appendChild(div_col2c);

				if(leveltier==3){
					//p2a.setAttribute("style","text-align: center;font-size:15px")
					//p2b.setAttribute("style","text-align: center;font-size:15px")
					//div_col2a.setAttribute("class","col-sm-3")
					div_col2b.setAttribute("class","col-sm-3")
					//p2c.setAttribute("style","text-align: center;font-size:15px")
					div_col2c.setAttribute("class","col-sm-3")

					var div_col2d=document.createElement("div");
	        div_col2d.setAttribute("class","col-sm-3")
					var p2d=document.createElement("p");
	        p2d.setAttribute("style","text-align: center")
	        p2d.textContent=("--");
	        div_col2d.appendChild(p2d);
					div_col2_row.appendChild(div_col2d);
				}


        var div_col3=document.createElement("div");
        div_col3.setAttribute("class","col-sm-5")
        div_col3.setAttribute("style","text-align: center;")
        var button=document.createElement("button");
        button.setAttribute("class","lsm_button lsm_button_buffer");
        var text=document.createTextNode(numberWithCommas(Math.floor(mbuffer/E)))
        button.appendChild(text);
        div_col3.appendChild(button);


        div_new_row.appendChild(div_col1);
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

    var lsm_button_size_ratio=(max_button_size-100)/filters_monkey.length;
    var cur_length=100;
    cur_length+=lsm_button_size_ratio;

    if (N<=(mbuffer/E))
    {
    	//nothing in LSM tree
//adding the buffer row
    var div_new_row=document.createElement("div");
    div_new_row.setAttribute("class","row");
		var div_empty = document.createElement("div");
		div_empty.setAttribute("class","col-sm-1");
		var div_col2 = document.createElement("div");
		div_col2.setAttribute("class", "col-sm-6");
		var div_col2_row = document.createElement("div");
		div_col2_row.setAttribute("class","row");
		div_col2.appendChild(div_col2_row);




        var div_col2a=document.createElement("div");
        div_col2a.setAttribute("class","col-sm-2")


        var div_col2b=document.createElement("div");
        div_col2b.setAttribute("class","col-sm-6")

        var div_col2c=document.createElement("div");
        div_col2c.setAttribute("class","col-sm-4");

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


        div_col2_row.appendChild(div_col2a);
        div_col2_row.appendChild(div_col2b);
        div_col2_row.appendChild(div_col2c);

				if(leveltier==3){
					div_col2a.setAttribute("class","col-sm-3")
					div_col2_row.removeChild(div_empty);
					div_col2b.setAttribute("class","col-sm-3")
					div_col2c.setAttribute("class","col-sm-3")
					var div_col2d=document.createElement("div");
	        div_col2d.setAttribute("class","col-sm-3")
	        div_col2d.appendChild(p2d);
					div_col2_row.appendChild(div_col2d);
				}

				div_new_row.appendChild(div_empty);
				div_new_row.appendChild(div_col2);
        div_new_row.appendChild(div_col3);

    result_div.appendChild(div_new_row);


	    var hr=document.createElement("hr");
	    result_div.appendChild(hr)

    }
    else
    {

		var last_is_smaller=false;
		var L = filters_monkey.length;
		if(leveltier == 3){
			L = filters_fluid_lsm_tree.length;
		}
	    for (var i=0;i<L;i++)
	    {
		    var div_new_row=document.createElement("div");
		    div_new_row.setAttribute("class","row");

				var div_empty = document.createElement("div");
				div_empty.setAttribute("class","col-sm-1");
				var div_col2 = document.createElement("div");
				div_col2.setAttribute("class", "col-sm-6");
				var div_col2_row = document.createElement("div");
				div_col2_row.setAttribute("class","row");
				div_col2.appendChild(div_col2_row);

				var div_col2a=document.createElement("div");
        div_col2a.setAttribute("class","col-sm-2")
				var p2a=document.createElement("p");
				p2a.setAttribute("style","text-align: center;")
				p2a.textContent=((i+1)+"")
				div_col2a.appendChild(p2a);


        var div_col2b=document.createElement("div");
        div_col2b.setAttribute("class","col-sm-6")
				var MSD2b=getMostSignificantDigit(filters_baseline[i].fp);
				if (MSD2b>10)
					MSD2b=10;
				var message2b="The false positive rate for Bloom filters at Level "+(i+1)+" is "+(filters_baseline[i].fp*100).toFixed(MSD2b+1)+"%. This entails "+(filters_baseline[i].mem/filters_baseline[i].nokeys).toFixed(2)+" bits-per-element, resulting in a total of "+formatBytes(filters_baseline[i].mem/8,1)+" for Bloom filters at this level out of "+formatBytes(mfilter)+" for Bloom filters across all levels."
				var span2b=document.createElement("span");
				span2b.setAttribute("data-tooltip",message2b);
				span2b.setAttribute("data-tooltip-position","left")
			var p2b=document.createElement("p");
			p2b.setAttribute("style","text-align: center;")
				if (MSD2b<1)
					MSD2b=1;
			p2b.textContent=(filters_baseline[i].fp*100).toFixed(MSD2b-1)+"%"
				span2b.appendChild(p2b);
			div_col2b.appendChild(span2b);

        var div_col2c=document.createElement("div");
        div_col2c.setAttribute("class","col-sm-4");

	        	var MSD2c=getMostSignificantDigit(filters_monkey[i].fp);
	        	if (MSD2c>10)
	        		MSD2c=10;
	        	var message2c="The false positive rate for Bloom filters at Level "+(i+1)+" is "+(filters_monkey[i].fp*100).toFixed(MSD2c+1)+"%. This entails "+(filters_monkey[i].mem/filters_monkey[i].nokeys).toFixed(2)+" bits-per-element, resulting in a total of "+formatBytes(filters_monkey[i].mem/8,1)+" for Bloom filters at this level out of "+formatBytes(mfilter)+" for Bloom filters across all levels."
	        	var span2c=document.createElement("span");
	        	span2c.setAttribute("data-tooltip",message2c);
	        	span2c.setAttribute("data-tooltip-position","right")
			    var p2c=document.createElement("p");
			    p2c.setAttribute("style","text-align: center;")
	        	if (MSD2c<1)
	        		MSD2c=1;
		    	p2c.textContent=(filters_monkey[i].fp*100).toFixed(MSD2c-1)+"%"
	        	span2c.appendChild(p2c);
			    div_col2c.appendChild(span2c);


		   	    var div_col3=document.createElement("div");
			    div_col3.setAttribute("class","col-sm-5")
			    div_col3.setAttribute("style","text-align: center;")

			    var levelcss=i+1;
			    if (L<5)
			    	levelcss=5-L+1+i;
                // console.log(i+":"+levelcss)
								var n;
                if (i == filters_monkey.length-1) {
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
			                if (last_is_smaller)
			    	            button.setAttribute("class","lsm_button_not_solid");
			                else
			    	            button.setAttribute("class","lsm_button");
												if(maxRuns >= 7){
													button.setAttribute("style","width: "+(cur_length- 19.27)/6+"px; height: 36px");
												}else{
													button.setAttribute("style","width: "+cur_length/n+"px; height: 36px");
												}


                            var message = "At this level, each run contains "+numberWithCommas(Math.floor(filters_monkey[i].nokeys/maxRuns)+" entries"+((last_is_smaller)? " (level is not full)" : ""));
                            button.setAttribute("data-tooltip", message);
                            button.setAttribute("data-tooltip-position", "left");
                            div_col3.appendChild(button);
                        }
                    }
                    cur_length+=lsm_button_size_ratio;

				if (i==(L-2)) {
				    if ((leveltier != 3 && T*filters_monkey[i].nokeys>filters_monkey[i+1].nokeys) || (leveltier == 3 && T*filters_fluid_lsm_tree[i].nokeys>filters_fluid_lsm_tree[i+1].nokeys))
				    	last_is_smaller=true;
                }


								div_col2_row.appendChild(div_col2a);
				        div_col2_row.appendChild(div_col2b);
				        div_col2_row.appendChild(div_col2c);

								if(leveltier==3){
									//p2a.setAttribute("style","text-align: center;font-size:15px")
									//p2b.setAttribute("style","text-align: center;font-size:15px")
									div_col2b.setAttribute("class","col-sm-3")
									//p2c.setAttribute("style","text-align: center;font-size:15px")
									div_col2c.setAttribute("class","col-sm-3")

									var div_col2d=document.createElement("div");
					        div_col2d.setAttribute("class","col-sm-3")



									var MSD2d=getMostSignificantDigit(filters_fluid_lsm_tree[i].fp);
				        	if (MSD2d>10)
				        		MSD2d=10;
				        	var message2d="The false positive rate for Bloom filters at Level "+(i+1)+" is "+(filters_fluid_lsm_tree[i].fp*100).toFixed(MSD2d+1)+"%. This entails "+(filters_fluid_lsm_tree[i].mem/filters_fluid_lsm_tree[i].nokeys).toFixed(2)+" bits-per-element, resulting in a total of "+formatBytes(filters_fluid_lsm_tree[i].mem/8,1)+" for Bloom filters at this level out of "+formatBytes(mfilter)+" for Bloom filters across all levels."
				        	var span2d=document.createElement("span");
				        	span2d.setAttribute("data-tooltip",message2d);
				        	span2d.setAttribute("data-tooltip-position","right")
									var p2d=document.createElement("p");
									p2d.setAttribute("style","text-align: center;")
				        	if (MSD2d<1)
				        		MSD2d=1;
					    	p2d.textContent=(filters_fluid_lsm_tree[i].fp*100).toFixed(MSD2d-1)+"%"
				        	span2d.appendChild(p2d);
						    div_col2d.appendChild(span2d);

									div_col2_row.appendChild(div_col2d);
								}

								div_new_row.appendChild(div_empty);
								div_new_row.appendChild(div_col2);
				        div_new_row.appendChild(div_col3);
	    	  result_div.appendChild(div_new_row);


	    }



	    var hr=document.createElement("hr");
	    result_div.appendChild(hr)


	    //total read cost line
	    var div_new_row=document.createElement("div");
	    div_new_row.setAttribute("class","row")

	    	var Rbaseline=eval_R(filters_baseline, leveltier, T);
	    	var Rmonkey=eval_R(filters_monkey, leveltier, T);
				var RfluidLSMTree;
				if(leveltier == 3){
					RfluidLSMTree = eval_R(filters_fluid_lsm_tree, leveltier, T, K, Z)
				}

				var div_col1 = document.createElement("div");
				div_col1.setAttribute("class","col-sm-2")
			var p1=document.createElement("p");
			p1.setAttribute("style","text-align: right;")
			p1.textContent=("Lookup cost:")
			div_col1.appendChild(p1);



	        var message;
	        if (leveltier==0) {
                message="The average I/O cost for a zero-result lookup is the sum of the false positives rates across all runs. With tiering, there are (T-1) runs per level, where T is the size ratio. Therefore, lookup cost is the sum of false positive rate prescriptions across all levels multiplied by (T-1).";
            } else if (leveltier==1) {
                message="The average I/O cost for a zero-result lookup is the sum of the false positives rates across all levels. ";
            } else if (leveltier==2) {
                message="The average I/O cost for a zero-result lookup is the sum of the false positives rates across all runs. With lazy leveling, there are (T-1) runs per level on levels except the last one, where T is the size ratio. Therefore, lookup cost is the sum of false positive rate prescriptions across all levels except the last one multiplied by (T-1) plus the false positive rate of the last level.";
            } else if (leveltier==3) {
							message="The average I/O cost for a zero-result lookup is the sum of the false positives rates across all runs. With fluid LSM-Tree, there are K runs per level on levels except the last one, while there are Z runs in the last level. Therefore, lookup cost is the sum of false positive rate prescriptions across all levels except the last one multiplied by K plus the multiplication between Z and the false positive rate of the last level.";
						}



	        var div_col2=document.createElement("div");
		    div_col2.setAttribute("class","col-sm-3")
	        var span2=document.createElement("span");
	        span2.setAttribute("data-tooltip",message);
	        span2.setAttribute("data-tooltip-position","bottom")
		    var p2=document.createElement("p");
		    p2.setAttribute("style","text-align: center;")
	        var em2=document.createElement("em");
	        em2.textContent=(Rbaseline.toFixed(3)+" I/Os")
	        p2.appendChild(em2);
	        span2.appendChild(p2);
	        div_col2.appendChild(span2);

	        var div_col3=document.createElement("div");
	        div_col3.setAttribute("class","col-sm-2")
	        var span3=document.createElement("span");
	        span3.setAttribute("data-tooltip",message);
	        span3.setAttribute("data-tooltip-position","bottom")
	        var p3=document.createElement("p");
	        p3.setAttribute("style","text-align: center;")
	        var em3=document.createElement("em");
	        em3.textContent=(Rmonkey.toFixed(3)+" I/Os")
	        p3.appendChild(em3);
	        span3.appendChild(p3);
	        div_col3.appendChild(span3);

					div_new_row.appendChild(div_col1);
		    	div_new_row.appendChild(div_col2);
		    	div_new_row.appendChild(div_col3);

            var speedup=(Rbaseline/Rmonkey).toFixed(2);

            if (isNaN(speedup))
                speedup="1.0"


					if(leveltier != 3){
						message=" Lookups are "+speedup+"x faster with Monkey. Update cost is the same as the state-of-the-art."
					}else{
						var fluidSpeedUp=(Rbaseline/RfluidLSMTree).toFixed(2);
						message=" Lazy leveling/Fluid LSM-Tree are"+speedup+"x/"+fluidSpeedUp+"x faster."

						div_col2.setAttribute("style","width:12.5%")
						div_col3.setAttribute("style","width:12.5%")

						var div_col3a=document.createElement("div");
		        div_col3a.setAttribute("class","col-sm-2")
						div_col3a.setAttribute("style","width:12.5%")
		        var span3a=document.createElement("span");
		        span3a.setAttribute("data-tooltip",message);
		        span3a.setAttribute("data-tooltip-position","bottom")
		        var p3a=document.createElement("p");
		        p3a.setAttribute("style","text-align: center;")
		        var em3a=document.createElement("em");
		        em3a.textContent=(RfluidLSMTree.toFixed(3)+" I/Os")
		        p3a.appendChild(em3a);
		        span3a.appendChild(p3a);
		        div_col3a.appendChild(span3a);

						div_new_row.appendChild(div_col3a);
					}

	        /*if (leveltier==0)
	            message+=" TIERING: When a run is flushed to the next level it is NOT merged with the already existing runs (if any). Hence, every time a run is pushed to the next level it is written only once.";
	        else if (leveltier==1)
	            message+=" LEVELING: When a run is flushed to the next level it is ALWAYS merged with the already existing runs (if any). Hence, every flush causes (up to) size_ratio-1 merges in the next level.";
            */
						var div_col4=document.createElement("div");
		        div_col4.setAttribute("class","col-sm-5")
						div_col4.setAttribute("style","overflow:  visible;height: 3vh;float:right")
	        var span4=document.createElement("span");
	        span4.setAttribute("data-tooltip",message);
	        span4.setAttribute("data-tooltip-position","bottom")
	        var p4=document.createElement("p");
	        p4.setAttribute("style","text-align: center;")
	        var em4=document.createElement("em");
	        em4.textContent=("Monkey lookups are "+speedup+"x faster!")
					if(leveltier != 3){
						em4.textContent=("Monkey lookups are "+speedup+"x faster!")
					}else{
						em4.textContent=("Lazy leveling/Fluid LSM-Tree lookups are "+speedup+"x/" + fluidSpeedUp + "x faster!")
					}
	        p4.appendChild(em4);
	        span4.appendChild(p4);
	        div_col4.appendChild(span4);

					div_new_row.appendChild(div_col4);


		    result_div.appendChild(div_new_row);

	       //total write cost line
	       var div_new_row=document.createElement("div");
	       div_new_row.setAttribute("class","row")

	        var W;
					var W_lazy;
	        var entries_per_page=Math.floor(P/E);
					W = (1/entries_per_page)*((T-1)*(filters_monkey.length - 1)/(K + 1) + (T-1)/(Z + 1));
					W_lazy = (1/entries_per_page)*((T-1)*(filters_monkey.length - 1)/T + (T-1)/2);



	        var div_col1=document.createElement("div");
	        div_col1.setAttribute("class","col-sm-2")
	        var p1=document.createElement("p");
	        p1.setAttribute("style","text-align: right;")
	        p1.textContent=("Update cost:")
	        div_col1.appendChild(p1);


	        var div_col2=document.createElement("div");
	        div_col2.setAttribute("class","col-sm-3")
	        var p2=document.createElement("p");
	        p2.setAttribute("style","text-align: center;")
					var em2=document.createElement("em");
					em2.textContent=(W.toFixed(3)+" I/Os")
	        p2.appendChild(em2);
	        div_col2.appendChild(p2);

	        var div_col3=document.createElement("div");
	        div_col3.setAttribute("class","col-sm-2")
	        var p3=document.createElement("p");
	        p3.setAttribute("style","text-align: center;")
					var em3=document.createElement("em");
					em3.textContent=(W.toFixed(3)+" I/Os")
	        p3.appendChild(em3);
	        div_col3.appendChild(p3);

					div_new_row.appendChild(div_col1);
	        div_new_row.appendChild(div_col2);
	        div_new_row.appendChild(div_col3);


					if(leveltier == 3){

						div_col2.setAttribute("style","width:12.5%")
						div_col3.setAttribute("style","width:12.5%")
						em3.textContent=(W_lazy.toFixed(3)+" I/Os")

						var div_col3a=document.createElement("div");
		        div_col3a.setAttribute("class","col-sm-2")
						div_col3a.setAttribute("style","width:12.5%")
		        var span3a=document.createElement("span");
		        span3a.setAttribute("data-tooltip",message);
		        span3a.setAttribute("data-tooltip-position","bottom")
		        var p3a=document.createElement("p");
		        p3a.setAttribute("style","text-align: center;")
		        var em3a=document.createElement("em");
		        em3a.textContent=(W.toFixed(3)+" I/Os")
		        p3a.appendChild(em3a);
		        span3a.appendChild(p3a);
		        div_col3a.appendChild(span3a);

						div_new_row.appendChild(div_col3a);
					}

					var div_col4=document.createElement("div");
					div_col4.setAttribute("class","col-sm-5")
					var p4=document.createElement("p");
					p4.setAttribute("style","text-align: center;")
					p4.textContent=("")
					div_col4.appendChild(p4);


	        div_new_row.appendChild(div_col4);
	        result_div.appendChild(div_new_row);

            //total write cost line
           var div_new_row=document.createElement("div");
           div_new_row.setAttribute("class","row")

            var total_mem = mfilter + mbuffer;
            //alert("mfilter " + mfilter + "   mbuffer" + mbuffer + "  total_mem " + total_mem);

            var div_col1=document.createElement("div");
            div_col1.setAttribute("class","col-sm-2")
            var p1=document.createElement("p");
            p1.setAttribute("style","text-align: right;")
            p1.textContent=("Main memory:")
            div_col1.appendChild(p1);

            var div_col2=document.createElement("div");
            div_col2.setAttribute("class","col-sm-3")
            var p2=document.createElement("p");
            p2.setAttribute("style","text-align: center;")
						var em2=document.createElement("em");
		        em2.textContent=(formatBytes(total_mem,2))
		        p2.appendChild(em2);
            div_col2.appendChild(p2);

            var div_col3=document.createElement("div");
            div_col3.setAttribute("class","col-sm-2")
            var p3=document.createElement("p");
            p3.setAttribute("style","text-align: center;")
						var em3=document.createElement("em");
		        em3.textContent=(formatBytes(total_mem,2))
		        p3.appendChild(em3);
            div_col3.appendChild(p3);

						div_new_row.appendChild(div_col1);
            div_new_row.appendChild(div_col2);
            div_new_row.appendChild(div_col3);

						if(leveltier == 3){

							div_col2.setAttribute("style","width:12.5%")
							div_col3.setAttribute("style","width:12.5%")

							var div_col3a=document.createElement("div");
			        div_col3a.setAttribute("class","col-sm-2")
							div_col3a.setAttribute("style","width:12.5%")
			        var p3a=document.createElement("p");
			        p3a.setAttribute("style","text-align: center;")
			        var em3a=document.createElement("em");
			        em3a.textContent=(formatBytes(total_mem,2))
			        p3a.appendChild(em3a);
			        div_col3a.appendChild(p3a);

							div_new_row.appendChild(div_col3a);
						}

            var div_col4=document.createElement("div");
            div_col4.setAttribute("class","col-sm-5")
            var p4=document.createElement("p");
            p4.setAttribute("style","text-align: center;")
            p4.textContent=("")
            div_col4.appendChild(p4);



            div_new_row.appendChild(div_col4);
            result_div.appendChild(div_new_row);

	}
    // var hr=document.createElement("hr");
    // result_div.appendChild(hr)

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
    var valid;
    var throughput;
    var num_levels;
}

// returns write cost
function get_W(N, T, B, P, leveled) {
    if (P * B >= N) {
        return 0;
    }

    var loga = Math.ceil(Math.log(N / (B * P)) / Math.log(T));
    var denom = leveled ? 2.0 : T;
    var all =  ((T - 1.0) / denom);
    all = loga * all ;
  //  all = all + 1;
    all = all / B;
    return all;
}

// read cost under state-of-the-art approach
function get_R_uniform_strategy(M, T, N, B, P, leveled) {
    M = M * (8 * 1024);

    if (P * B >= N) {
        //return 0;
    }

    if (M <= 0) {
        var c = Math.log(N / (B * P));
        var  b = c / Math.log(T);
        var  a = 1.0 + b ;
        a *= leveled ? 1 : (T - 1);
        return a;
    }

    //double L = !leveled ? floor(1 +  log(N/(B * P)) / log(T)) : ceil(1 +  log(N/(B * P)) / log(T));
    //double L = !leveled ? floor(1 +  log(N/(B * P)) / log(T)) : ceil(1 +  log(N/(B * P)) / log(T));
    var L = 1 +  Math.log(N/ (B * P)) / Math.log(T);
    if (!leveled) {
        //L--;
    }
    var T_part = (T - 1) / (T * ( 1 - Math.pow(T, -L) ));

    var exponent = (M / N) * Math.log(2) * Math.log(2) * T_part;

    var EULER = 2.71828182845904523536;
    var bottom = Math.pow(EULER, exponent);
    var R = L / bottom;
    if (!leveled) {
        R *= T - 1;
    }
    return R;
}


// read cost under MonKey
function get_accurate_R(M, T, N, B, P, leveled) {

    if (P * B >= N) {
        //return 0;
    }

    if (M <= 0) {
        var c = Math.log(N / (B * P));
        var b = c / Math.log(T);
        var a = 1.0 +  b ;
        a *= leveled ? 1 : (T - 1);
        return a;
    }



    var F = N / Math.pow(Math.log(2), 2);
    F *= T / Math.pow((T - 1), 2);
    F *= Math.log(T);
    F /= 8 * 1024;

    var j = Math.max(Math.log(F / M) / Math.log(T), 0.0);
    var L = 1 + Math.log(N / (B * P)) / Math.log(T) - j;
    //printf("%f\n", L);
    //double L = !leveled ? floor(1 +  log(N/(B * P)) / log(T)) : ceil(1 +  log(N/(B * P)) / log(T));

    N /= Math.pow(T, j);

    var Y = 1 - L / Math.pow(T, L);
    Y += (1 - Math.pow(1/T, L - 1.0)) / (T - 1.0);
    Y = Math.pow(T, Y);

    var S = leveled ? (T - 1) : 1.0;
    var X = Math.pow(T, L) * (S / (Math.pow(T, L) - 1));;
    X = Math.pow(X, 1 - Math.pow(1/T, L));

    var E = ((M * 8 * 1024) / N) * ((T - 1) / T) * Math.pow(Math.log(2), 2);

    var EULER = 2.71828182845904523536;
    E = Math.pow(EULER, E);

    var product = Y / (X * E);

    var raised = Math.pow(product, Math.pow(T, L) / (Math.pow(T, L) - 1));

    return raised + j * (leveled ? 1 : T - 1.0);
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
function find_optimal_R(input_conf, constant_buffer_size = -1, use_new_bloom_tuning = true) {
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
    conf.valid=input_conf.valid;
    conf.throughput=input_conf.throughput;
    conf.num_levels=input_conf.num_levels;

    var N = conf.N;
    var B = conf.B;
    var E = conf.E;
    var M = conf.M;
    var W = conf.W;
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
                var N_used = N * (T-1) / (T);
                var num_levels = 1 + Math.log(N_used / (B * P)) / Math.log(T);
                var mem_for_bloom = M - P * B * E;

                var current_R = use_new_bloom_tuning ? get_accurate_R(mem_for_bloom / 1024, T, N_used, B, P, leveled) :
                        get_R_uniform_strategy(mem_for_bloom / 1024, T, N_used, B, P, leveled);

                var current_W = get_W(N_used, T, B, P, leveled);
                if (mem_for_bloom >= 0 && current_R < min_R && current_W <= W) {
                    min_R = current_R;
                    best_T = T;
                    best_P = P;
                    min_W = current_W;
                    is_leveled = leveled;
                    conf.W = current_W;
                    conf.valid = true;
                    conf.num_levels = 1 +  Math.log(N_used/(B * P)) / Math.log(T);
                }
            }
        }
    }
    conf.T = best_T;
    conf.P = best_P;
    conf.L = is_leveled;
    conf.R = min_R;
    return conf;
}




function getPoint(tieringVsLeveling, T, mfilter, conf, use_monkey) {
                var meetingR;
                var N_used = conf.N * (T-1) / T;
                if (use_monkey) {
                    meetingR = get_accurate_R(mfilter / 1024, T, N_used, conf.B, conf.P, tieringVsLeveling);
                }
                else {
                    meetingR = get_R_uniform_strategy(mfilter  / 1024, T, N_used, conf.B, conf.P, tieringVsLeveling);
                }
                var meetingW = get_W(N_used, T, conf.B, conf.P, tieringVsLeveling);
                return {W: meetingW, R: meetingR};
}



// Find the corresponding optimal read cost for different upper bounds on write cost
function print_csv_experiment(input_conf, num_commas, print_details, fix_buffer_size = -1, use_monkey = true, smoothing = false, differentiate_tiered_leveled = true) {

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
    conf.valid=input_conf.valid;
    conf.throughput=input_conf.throughput;
    conf.num_levels=input_conf.num_levels;

    var prev_W = -1;
    var prev_L = 0;

    var array_to_print = new Array();

    var meeting_W = getPoint(0, 2, conf.M, conf, true).W;
    var min_W = getPoint(0, 4, conf.M, conf, true).W;
    var max_W = getPoint(1, 10, conf.M, conf, true).W;

    var starting_W = min_W * 0.01;
    var finishing_W = max_W * 10;
    var steps_for_leveling = (finishing_W - meeting_W) / 130;

    //print_write_optimized_extreme(conf, use_new_strategy);
    // console.log(input_conf)
    // console.log(conf)

    var last_c = LSM_config();
    for (var W = starting_W; W <= finishing_W; W += 0.001) {
        conf.W = W;

        var c = find_optimal_R(conf, fix_buffer_size, use_monkey);
        // console.log(c)
        // print_csv_line(c,num_commas,print_details,differentiate_tiered_leveled);

        // document.getElementById("explanation").value += (c.W.toFixed(4) + "\t\t" + c.R.toFixed(2) + "\t\t" + c.L + "\t\t" + c.T + "\t\t" + c.valid + "\n");

        // return;

        if ((prev_W == c.W && prev_L == c.L) || !c.valid) {
            continue;
        }
        prev_W = c.W;
        prev_L = c.L;

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



function MergeByFliudLSMTree(){
	document.getElementById("Fluid LSM-Tree K").value = 9;
	document.getElementById("input-group-Fluid-K").style.display = '';
	document.getElementById("Fluid LSM-Tree Z").value = 1;
	document.getElementById("input-group-Fluid-Z").style.display = '';
	re_run(event, 'input7');
}

function MergeNotyFliudLSMTree(){
	document.getElementById("input-group-Fluid-K").style.display = 'none';
	document.getElementById("input-group-Fluid-Z").style.display = 'none';
	re_run(event, 'input7');
}

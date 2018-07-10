var timer=null;

function re_run(e) {

    if(timer){
        clearTimeout(timer);
        timer = null;
    }

    var event = e || event;
    // console.log(event)

    var x = event.which || event.keyCode;  // Use either which or keyCode, depending on browser support
    // console.log(event.timeStamp)
    // console.log(x)

    if (!((x==38 || x==40) && (event.target.id=="E" || event.target.id=="mbuffer" || event.target.id=="P" || event.target.id=="T")))
    {
	    if (x>=37 && x<=40 || x==9 || x==16 || x==17 || x==91 || x==18 || x==65)
	    {
	        // if (event.target.id=="N" || event.target.id=="mfilter")
	        console.log("User is clicking ["+x+"] on arrows or tab (9) or shift (16) or ctrl (17) or cmd (91) or alt (18). No rerun.")
	        return;
	    }
	}

    if(event.target.id=="mfilter")
    {
        var mfilter = parseInt(document.getElementById("mfilter").value.replace(/\D/g,''),10);
        document.getElementById("mfilter").value=numberWithCommas(mfilter)
        // console.log(numberWithCommas(N))
    }

    if(event.target.id=="mbuffer")
    {
        var mbuffer = parseFloat(document.getElementById("mbuffer").value.replace(/\D/g,''))*1048576;
        document.getElementById("mbuffer").value=mbuffer/1048576;
        // console.log(numberWithCommas(N))
    }

    if(event.target.id=="E")
    {
        var mfilter = parseInt(document.getElementById("E").value.replace(/\D/g,''),10);
        document.getElementById("E").value=E
        // console.log(numberWithCommas(N))
    }

    if(event.target.id=="T")
    {
        var T = parseInt(document.getElementById("T").value.replace(/\D/g,''),10);
        if (T<2)
            document.getElementById("T").value=2;
        // console.log(numberWithCommas(T))
    }

    if(event.target.id=="N")
    {
        var N = parseInt(document.getElementById("N").value.replace(/\D/g,''),10);
        document.getElementById("N").value=numberWithCommas(N)
        var mbuffer = parseFloat(document.getElementById("mbuffer").value.replace(/\D/g,''))*1048576;
        var E = parseInt(document.getElementById("E").value.replace(/\D/g,''),10);
        var T = parseInt(document.getElementById("T").value.replace(/\D/g,''),10);
        var L = Math.ceil(Math.log(N*E*(T - 1)/mbuffer + 1)/Math.log(T)) + 1;
        document.getElementById("L").value=L;
        // console.log(numberWithCommas(N))
    }

    if(event.target.id=="L")
    {
        var L = parseInt(document.getElementById("L").value.replace(/\D/g,''),10);
        document.getElementById("L").value=L
        var mbuffer = parseFloat(document.getElementById("mbuffer").value.replace(/\D/g,''))*1048576;
        var E = parseInt(document.getElementById("E").value.replace(/\D/g,''),10);
        var T = parseInt(document.getElementById("T").value.replace(/\D/g,''),10);
        var N = Math.floor(mbuffer*(Math.pow(T, L + 1) - 1)/(E*(T - 1)));
        document.getElementById("N").value=numberWithCommas(N)
        // console.log(numberWithCommas(N))
    }


    timer=setTimeout(re_run_now,250);

}


function re_run_now() {

    var inputParameters = parseInputTextBoxes();

    var N=inputParameters.N;
    var E=inputParameters.E;
    var mbuffer=inputParameters.mbuffer;
    var T=inputParameters.T;
    var mfilter=inputParameters.mfilter;
    var P=inputParameters.P;
    var leveltier=inputParameters.leveltier;
    var fluidK=inputParameters.fluidK;
    var fluidZ=inputParameters.fluidZ;

    if (!isNaN(N))
        document.getElementById("N").value=numberWithCommas(N);
    if (!isNaN(E))
        document.getElementById("E").value=E;
    if (!isNaN(mbuffer))
        document.getElementById("mbuffer").value=mbuffer/1048576;
    if (!isNaN(T))
        document.getElementById("T").value=T;
    if (!isNaN(mfilter))
        document.getElementById("mfilter").value=numberWithCommas(mfilter/1048576)
    if (!isNaN(P))
        document.getElementById("P").value=P;

    if(leveltier != 3){
      document.getElementById("Fluid LSM-Tree K").value=1;
      fluidK = 1;
      document.getElementById("Fluid LSM-Tree Z").value=1;
      fluidZ = 1;
    }

    if (!isNaN(fluidK)){
      if(fluidK <= 0){
        alert("K="+fluidK+" is too small")
        console.log("K is too small: "+fluidK)
        document.getElementById("Fluid LSM-Tree K").value=1;
        fluidK = 1;
      }else if(fluidK >= T){
        alert("K="+fluidK+" is too large")
        console.log("K is too large: "+fluidK)
        document.getElementById("Fluid LSM-Tree K").value=T-1;
        fluidK = T-1;
      }else{
        document.getElementById("Fluid LSM-Tree K").value=fluidK;
      }
    }

    if (!isNaN(fluidZ)){
      if(fluidZ <= 0){
        alert("Z="+fluidZ+" is too small")
        console.log("Z is too small: "+fluidZ)
        document.getElementById("Fluid LSM-Tree Z").value=1;
        fluidZ=1;
      }else if(fluidZ >= T){
        alert("Z="+fluidZ+" is too large")
        console.log("Z is too large: "+fluidZ)
        document.getElementById("Fluid LSM-Tree Z").value=T-1;
        fluidZ=T-1;
      }else{
        document.getElementById("Fluid LSM-Tree Z").value=fluidZ;
      }
    }

    if(leveltier != 3){
      console.log("Running with: N="+N+", E="+E+", buffer="+(mbuffer/1048576)+", T="+T+", P="+P+", filter="+(mfilter/1048576)+", isLeveled="+leveltier)
    }else{
      console.log("Running with: N="+N+", E="+E+", buffer="+(mbuffer/1048576)+", T="+T+", P="+P+", filter="+(mfilter/1048576)+", K="+fluidK+", Z="+fluidZ+" in fluid LSM-Tree.")
    }

    // console.log("")
    // console.log("N="+N)
    // console.log("E="+E)
    // console.log("mbuffer="+mbuffer)
    // console.log("T="+T)
    // console.log("P="+P)
    // console.log("mfilter="+mfilter)
    // console.log("leveltier="+leveltier)

    var color='#777';
    document.getElementById("scenario1").style.background=color;
    document.getElementById("scenario2").style.background=color;
    document.getElementById("scenario3").style.background=color;

    if (document.getElementById("N").value=="" || document.getElementById("E").value=="" || document.getElementById("T").value=="" || document.getElementById("P").value==""
        || document.getElementById("mbuffer").value=="" || document.getElementById("mfilter").value=="" || isNaN(N)
        || isNaN(E) || isNaN(mbuffer) || isNaN(T) || isNaN(mfilter) || isNaN(P) || isNaN(leveltier))
    {

        return;
    }

    if(leveltier == 3 && (document.getElementById("Fluid LSM-Tree K").value=="" || document.getElementById("Fluid LSM-Tree Z").value=="" || isNaN(fluidK) || isNaN(fluidZ))){
      return;
    }


    if (T<2)
    {
        alert("T="+T+" is too small. It should be 2 or higher.")
        console.log("T is too small: "+T)
    }
    else if (mbuffer==0)
    {
        alert("Buffer="+mbuffer+" is too small.")
        console.log("T is too small: "+T)
    }
    else
    {
        clickbloomTuningButton(false)
    }
}

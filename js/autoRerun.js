var timer=null;

function re_run() {
    // console.log(event)
    // console.log(event.target)
    var x = event.which || event.keyCode;  // Use either which or keyCode, depending on browser support
    // console.log(event.timeStamp)
    // console.log(x)
    // console.log(String.fromCharCode(x))
    // console.log("running in a few")
    
    if (x>=37 && x<=40 || x==9 || x==16 || x==17 || x==91 || x==18 || x==65)
    {
        console.log("User is clicking on arrows or tab (9) or shift (16) or ctrl (17) or cmd (91) or alt (18). No rerun.")
        return;
    }

    if(event.target.id=="N")
    {
        var N = parseInt(document.getElementById("N").value.replace(/\D/g,''),10);
        document.getElementById("N").value=numberWithCommas(N)
        // console.log(numberWithCommas(N))
    }


    // console.log(timer)
    if(timer){
        clearTimeout(timer);
        timer = null;
    }
    timer=setTimeout(re_run_now,200);

    // re_run_now();        
}


function re_run_now() {

    var N = parseInt(document.getElementById("N").value.replace(/\D/g,''),10);
    var E = parseInt(document.getElementById("E").value.replace(/\D/g,''),10);
    var mbuffer = parseFloat(document.getElementById("mbuffer").value.replace(/\D/g,''))*1048576;
    var T = parseInt(document.getElementById("T").value.replace(/\D/g,''), 10);
    var mfilter = parseFloat(document.getElementById("mfilter").value.replace(/\D/g,''))*1048576;
    var P = parseInt(document.getElementById("P").value.replace(/\D/g,''), 10);
    var leveltier = getRadioValueByName("ltradio");

    if (!isNaN(N))
        document.getElementById("N").value=numberWithCommas(N);
    if (!isNaN(E))
        document.getElementById("E").value=E;
    if (!isNaN(mbuffer))
        document.getElementById("mbuffer").value=mbuffer/1048576;
    if (!isNaN(T))
        document.getElementById("T").value=T;
    if (!isNaN(mfilter))
        document.getElementById("mfilter").value=mfilter/1048576;
    if (!isNaN(P))
        document.getElementById("P").value=P;

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


    if (T<2)
    {
        console.log("T is too small: "+T)
    }
    else
    {
        clickbloomTuningButton(false)
    }
}
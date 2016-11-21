function re_run() {

    var N = parseInt(document.getElementById("N").value,10);
    var E = parseInt(document.getElementById("E").value,10);
    var mbuffer = parseFloat(document.getElementById("mbuffer").value)*1048576;
    var T = parseInt(document.getElementById("T").value, 10);
    var mfilter = parseFloat(document.getElementById("mfilter").value)*1048576;
    var P = parseInt(document.getElementById("P").value, 10);
    var leveltier = getRadioValueByName("ltradio");


    console.log("")
    console.log("N="+N)
    console.log("E="+E)
    console.log("mbuffer="+mbuffer)
    console.log("T="+T)
    console.log("P="+P)
    console.log("mfilter="+mfilter)
    console.log("leveltier="+leveltier)

    document.getElementById("scenario1").style.background='#c51c30';   
    document.getElementById("scenario2").style.background='#c51c30';   
    document.getElementById("scenario3").style.background='#c51c30';   

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
        clickbloomTuningButton()
}
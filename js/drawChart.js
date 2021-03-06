
setTimeout(scenario1, 200);
// setTimeout(function(){
//   $(document).ready(function(){
//     $('[data-toggle="tooltip"]').tooltip({
//      'delay': { show: 0, hide: 0 }
//     });
//   });
// }, 50);

// $('a').tooltip({
     // 'delay': { show: 5000, hide: 3000 }
// });

function initChart() {

    var monkeyTiering = {
        y: [-1],
        x: [-1],
        mode: 'lines+markers',
        type: 'scatter',
        name: 'Monkey - Tiering',
        text: [ "" ],
        marker: { size: 7, symbol: 'circle' },
        line: {width: 1}
    };

    var monkeyLeveling = {
        y: [-1],
        x: [-1],
        mode: 'lines+markers',
        type: 'scatter',
        name: 'Monkey - Leveling',
        text: [ "" ],
        marker: { size: 7, symbol: 'circle'},
        line: {width: 1}
    };

    var monkeyLazyLeveling = {
        y: [-1],
        x: [-1],
        mode: 'lines+markers',
        type: 'scatter',
        name: 'Monkey - Lazy Leveling',
        text: [ "" ],
        marker: { size: 7, symbol: 'circle'},
        line: {width: 1}
    };

    var transitionOne = {
        y: [-1],
        x: [-1],
        mode: 'lines+markers',
        type: 'scatter',
        name: 'Transition 1',
        // text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
        text: [ "" ],
        marker: { size: 7, symbol: 'circle'},
        line: {dash: 'dashdot', width: 1}
    };

    var state_of_artTiering = {
        y: [-1],
        x: [-1],
        mode: 'lines+markers',
        type: 'scatter',
        name: 'State-of-the-art - Tiering',
        text: [ "" ],
        marker: { size: 7, symbol: 'square'  },
        line: {width: 1}
    };

    var state_of_artLeveling = {
        y: [-1],
        x: [-1],
        mode: 'lines+markers',
        type: 'scatter',
        name: 'State-of-the-art - Leveling',
        text: [ "" ],
        marker: { size: 7, symbol: 'square' },
        line: {width: 1}
    };

    var selected_point = {
        y: [-1],
        x: [-1],
        mode: 'markers',
        type: 'scatter',
        name: 'Selected point',
        text: [ "" ],
        marker: { size: 13, symbol: 'star', color: 'black' },
        line: {width: 1}
    };

  var monkeyPoint = {
        y: [-1],
        x: [-1],
        mode: 'markers',
        type: 'scatter',
        name: 'Selected Configuration',
        text: [ "" ],
        marker: { size: 13, symbol: 'star', color: 'black' },
        line: {width: 1}
    };

    var stateOfArtPoint = {
        y: [-1],
        x: [-1],
        mode: 'markers',
        type: 'scatter',
        name: 'State-of-the-art - Selected Configuration',
        text: [ "" ],
        marker: { size: 11, symbol: 'cross', color: 'black' },
        line: {width: 1}
    };

    var data = [ state_of_artTiering, state_of_artLeveling, stateOfArtPoint, monkeyTiering , monkeyLeveling, monkeyPoint];

    var xmin=0;
    var ymin=0;
    var xmax=2;
    var ymax=2;

    var layout =
    {
        xaxis: {
            title: 'Update cost (I/Os)',
            range: [ xmin, xmax ]
        },
        yaxis: {
            title: 'Lookup cost (I/Os)',
            range: [ymin, ymax]
        },
        //title:'Pareto frontiers for State-of-the-art and Monkey Tuning'
        title:''
    };

    Plotly.newPlot('myDiv', data, layout);


}

function addPoint(leveltier, T, mfilter, conf, monkeyTW, monkeyTR, monkeyT_Ratio, isOptimalFPR) {
                var meetingR;
                if (isOptimalFPR) {
                    meetingR = get_accurate_R(mfilter, T, conf.N, conf.K, conf.Z, conf.B, conf.P, leveltier);
                }
                else {
                    meetingR = get_R_uniform_strategy(mfilter, T, conf.N, conf.K, conf.Z, conf.B, conf.P, leveltier);
                }
                var meetingW = get_W(mfilter, conf.N, T, conf.K, conf.Z, conf.B, conf.P, conf.Mu, leveltier);
                    monkeyTW.push(meetingW.toExponential(4));
                    monkeyTR.push(meetingR.toExponential(3));
                    monkeyT_Ratio.push("Ratio: "+T);
                return {W: meetingW, R: meetingR};
}




function drawChart() {

    var inputParameters = parseInputTextBoxes();

    var N=inputParameters.N;
    var E=inputParameters.E;
    var mbuffer=inputParameters.mbuffer;
    var T=inputParameters.T;
    var K=inputParameters.fluidK;
    var Z=inputParameters.fluidZ;
    var Mu=inputParameters.Mu;
    var mfilter_per_entry=inputParameters.mfilter_per_entry;
    var mfilter = mfilter_per_entry*N/8;
    var P=inputParameters.P;
    var leveltier=inputParameters.leveltier;


    var conf = new LSM_config();
    conf.T=-1;
    conf.L=Math.ceil(Math.log(N*E*(T - 1)/mbuffer/T+ 1/T)/Math.log(T));
    conf.P=mbuffer / P;
    conf.N=N;
    conf.M=mbuffer+mfilter;
    conf.mfilter=mfilter;
    conf.E=E;
    conf.B=P/E;
    conf.K=K;
    conf.Z=Z;
    conf.Mu=Mu;
    conf.leveltier=leveltier;
    conf.Y=calc_Y(mfilter/N, K, Z, T, L);

    var smoothing=false;
//print_csv_experiment(input_conf, num_commas, print_details, fix_buffer_size = -1, use_monkey = true, smoothing = false, differentiate_tiered_leveled = true)
    var monkey_pareto = print_csv_experiment(conf, 0, false, mbuffer, true, smoothing)
    var state_of_art_pareto = print_csv_experiment(conf, 0, false, mbuffer, false, smoothing)

    var monkeyTW = new Array();
    var monkeyTR = new Array();
    var monkeyT_Ratio= new Array();
    var monkeyLW = new Array();
    var monkeyLR = new Array();
    var monkeyL_Ratio= new Array();
    //lazy Leveling
    var monkeyLazyLW = new Array();
    var monkeyLazyLR = new Array();
    var monkeyLazyL_Ratio= new Array();

    // console.log("isLeveled  " + isLeveled + " \n");
    var part1_monkey_point = getPoint(leveltier, T, mfilter, conf, inputParameters.isOptimalFPR);
    //var part1_state_of_the_art_point = getPoint(leveltier, T, mfilter, conf, 0);
    //console.log("leveltier  " + leveltier);
    // console.log(part1_monkey_point);
    // console.log(part1_state_of_the_art_point);

    /*var p1 = getPoint(leveltier, T, mfilter, conf, 0);
    var p2 = getPoint(!leveltier, T, mfilter, conf, 0);
    console.log("test");
    console.log(p1);
    console.log(p2);*/


    for (var i=0;i<monkey_pareto.length;i++)
    {
        if (monkey_pareto[i].leveltier==0)
        {
            if (monkey_pareto[i].T > 3) {
                monkeyTW.push(monkey_pareto[i].W.toExponential(4));
                monkeyTR.push(monkey_pareto[i].R.toExponential(3));
                monkeyT_Ratio.push("Ratio: " + monkey_pareto[i].T);
            }
        }
        else
        {
            if (monkeyLW.length==0)
            {
                //addPoint(0, 4, mfilter, conf, monkeyTW, monkeyTR, monkeyT_Ratio, 1);
                addPoint(0, 3, mfilter, conf, monkeyTW, monkeyTR, monkeyT_Ratio, 1);
                addPoint(0, 2, mfilter, conf, monkeyTW, monkeyTR, monkeyT_Ratio, 1);
                addPoint(1, 2, mfilter, conf, monkeyLW, monkeyLR, monkeyL_Ratio, 1);
                addPoint(1, 3, mfilter, conf, monkeyLW, monkeyLR, monkeyL_Ratio, 1);
                //addPoint(1, 4, mfilter, conf, monkeyLW, monkeyLR, monkeyL_Ratio, 1);
            }
            if( monkeyLazyLW.length==0){
              addPoint(2, 2, mfilter, conf, monkeyLazyLW, monkeyLazyLR, monkeyLazyL_Ratio, 1);
              addPoint(2, 3, mfilter, conf, monkeyLazyLW, monkeyLazyLR, monkeyLazyL_Ratio, 1);
            }
            if (monkey_pareto[i].T > 3) {
                monkeyLW.push(monkey_pareto[i].W.toExponential(4));
                monkeyLR.push(monkey_pareto[i].R.toExponential(3));
                monkeyL_Ratio.push("Ratio: " + monkey_pareto[i].T);
                addPoint(2, monkey_pareto[i].T, mfilter, conf, monkeyLazyLW, monkeyLazyLR, monkeyLazyL_Ratio, 1);
            }
        }
    }

    //transition
    var transitionOneW = new Array();
    var transitionOneR = new Array();
    var transitionOne_Ratio= new Array();
    //var transitionTwoW = new Array();
    //var transitionTwoR = new Array();
    //var transitionTwo_Ratio= new Array();

    var w_min = Number.MAX_VALUE;
    var w_min_index = -1;
    for(var i=0; i < monkeyLazyLW.length-1; i++){
      if(parseFloat(monkeyLazyLW[i]) < parseFloat(w_min)){
        w_min = monkeyLazyLW[i];
        w_min_index = i;
      }
    }

    var transitionT = parseInt(monkeyLazyL_Ratio[w_min_index].substr(6));
    var transitionK = transitionT - 1;
    conf.K = transitionK;
    //transition 1
    addPoint(2, transitionT, mfilter, conf, transitionOneW, transitionOneR, transitionOne_Ratio, 1);
    for(var tmpZ=2; tmpZ <= transitionT-1; tmpZ++){
      conf.Z = tmpZ;
      addPoint(3, transitionT, mfilter, conf, transitionOneW, transitionOneR, transitionOne_Ratio, 1);
    }
    //transition 2
    /*
    conf.Z = 1;
    addPoint(2, transitionT, mfilter, conf, transitionTwoW, transitionTwoR, transitionTwo_Ratio, 1);
    for(var tmpT=transitionT+1; tmpT <= monkey_pareto[monkey_pareto.length-1].T; tmpT++){
      addPoint(3, tmpT, mfilter, conf, transitionTwoW, transitionTwoR, transitionTwo_Ratio, 1);
    }*/

    var state_of_artTW = new Array();
    var state_of_artTR = new Array();
    var state_of_artT_Ratio= new Array();
    var state_of_artLW = new Array();
    var state_of_artLR = new Array();
    var state_of_artL_Ratio= new Array();
    var meetSoA_R,meetSoA_W;
    var highestSoA_R;
    for (var i=0;i<state_of_art_pareto.length;i++)
    {
        if (state_of_art_pareto[i].leveltier==0)
        {
            if (state_of_artTW.length == 0) {
                highestSoA_R = state_of_art_pareto[i].R;
                highestSoA_W = state_of_art_pareto[i].W;
            }
            if (state_of_art_pareto[i].T > 3) {
                state_of_artTW.push(state_of_art_pareto[i].W.toExponential(4));
                state_of_artTR.push(state_of_art_pareto[i].R.toExponential(3));
                state_of_artT_Ratio.push("Ratio: "+state_of_art_pareto[i].T);
            }
        }
        else
        {
            if (state_of_artLW.length==0)
            {
                //addPoint(0, 4, mfilter, conf, state_of_artTW, state_of_artTR, state_of_artT_Ratio, 0);
                addPoint(0, 3, mfilter, conf, state_of_artTW, state_of_artTR, state_of_artT_Ratio, 0);
                addPoint(0, 2, mfilter, conf, state_of_artTW, state_of_artTR, state_of_artT_Ratio, 0);
                var meeting_point = addPoint(1, 2, mfilter, conf, state_of_artLW, state_of_artLR, state_of_artL_Ratio, 0);
                addPoint(1, 3, mfilter, conf, state_of_artLW, state_of_artLR, state_of_artL_Ratio, 0);
                //addPoint(1, 4, mfilter, conf, state_of_artLW, state_of_artLR, state_of_artL_Ratio, 0);

                meetSoA_R = meeting_point.R;
                meetSoA_W = meeting_point.W;
            }
            if (state_of_art_pareto[i].T > 3) {
              state_of_artLW.push(state_of_art_pareto[i].W.toExponential(4));
               state_of_artLR.push(state_of_art_pareto[i].R.toExponential(3));
               state_of_artL_Ratio.push("Ratio: "+state_of_art_pareto[i].T);
            }
        }
    }





    var monkeyPoint = {
        y: [part1_monkey_point.R],
        x: [part1_monkey_point.W],
        mode: 'markers',
        type: 'scatter',
        name: 'Selected Configuration',
        text: [ "Ratio: "+T ],
        marker: { size: 13, symbol: 'star', color: 'black' },
        line: {width: 1}
    };
/*
    var stateOfArtPoint = {
        y: [part1_state_of_the_art_point.R],
        x: [part1_state_of_the_art_point.W],
        mode: 'markers',
        type: 'scatter',
        name: 'State-of-the-art - Selected Configuration',
        text: [ "Ratio: "+T ],
        marker: { size: 11, symbol: 'cross', color: 'black' },
        line: {width: 1}
    };

*/
    var monkeyTiering = {
        y: monkeyTR,
        x: monkeyTW,
        mode: 'lines+markers',
        type: 'scatter',
        name: 'Monkey - Tiering',
        text: monkeyT_Ratio,
        marker: { size: 7, symbol: 'circle' },
        line: {width: 1}
    };

    var monkeyLeveling = {
        y: monkeyLR,
        x: monkeyLW,
        mode: 'lines+markers',
        type: 'scatter',
        name: 'Monkey - Leveling',
        // text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
        text: monkeyL_Ratio,
        marker: { size: 7, symbol: 'circle'},
        line: {width: 1}
    };

    var monkeyLazyLeveling = {
        y: monkeyLazyLR,
        x: monkeyLazyLW,
        mode: 'lines+markers',
        type: 'scatter',
        name: 'Monkey - Lazy Leveling',
        // text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
        text: monkeyLazyL_Ratio,
        marker: { size: 7, symbol: 'circle'},
        line: {width: 1}
    };

    var transitionOne = {
        y: transitionOneR,
        x: transitionOneW,
        mode: 'lines+markers',
        type: 'scatter',
        name: 'Transition 1',
        // text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
        text: transitionOne_Ratio,
        marker: { size: 7, symbol: 'circle'},
        line: {dash: 'dashdot', width: 1}
    };
/*
    var transitionTwo = {
        y: transitionTwoR,
        x: transitionTwoW,
        mode: 'lines+markers',
        type: 'scatter',
        name: 'Transition 2',
        // text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
        text: transitionOne_Ratio,
        marker: { size: 7, symbol: 'circle'},
        line: {dash: 'dashdot', width: 1}
    };
*/
    var state_of_artTiering = {
        y: state_of_artTR,
        x: state_of_artTW,
        mode: 'lines+markers',
        type: 'scatter',
        name: 'State-of-the-art - Tiering',
        text: state_of_artT_Ratio,
        marker: { size: 7, symbol: 'square' },
        line: {width: 1}
    };

    var state_of_artLeveling = {
        y: state_of_artLR,
        x: state_of_artLW,
        mode: 'lines+markers',
        type: 'scatter',
        name: 'State-of-the-art - Leveling',
        // text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
        text: state_of_artL_Ratio,
        marker: { size: 7, symbol: 'square' },
        line: {width: 1}
    };



    var data = [ state_of_artTiering, state_of_artLeveling, monkeyTiering , transitionOne, monkeyLeveling, monkeyPoint, monkeyLazyLeveling];


    var xmin=0;
    var ymin=0;
    var xmax=meetSoA_W*3;
    var ymax=Math.min(highestSoA_R, meetSoA_R*3);


    var layout =
    {
        xaxis: {
            title: 'Update cost (I/Os)',
            range: [ xmin, xmax ]
        },
        yaxis: {
            title: 'Lookup cost (I/Os)',
            range: [ymin, ymax]
        },
        //title:'Pareto frontiers for State-of-the-art and Monkey Tuning'
        title:''
    };

    Plotly.newPlot('myDiv', data, layout);

}
